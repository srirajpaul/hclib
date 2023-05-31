
#ifndef SELECTOR_H
#define SELECTOR_H

#include "safe_buffer.h"
#include "hclib_bale_actor.h"
extern "C" {
#include "convey.h"
}

#ifndef NO_USE_BUFFER
#define USE_BUFFER
#endif

#define DONE_MARK -1
#define BUFFER_SIZE 1024
#ifndef ELASTIC_BUFFER_SIZE
#define ELASTIC_BUFFER_SIZE 128
#endif

namespace hclib {

template<typename T>
struct BufferPacket {
    T data;
    int64_t rank;

    BufferPacket() : rank(0) {}
    BufferPacket(T data, int64_t rank) : data(data), rank(rank) {}
};

template<typename T, int SIZE>
class Mailbox {

    hclib::conveyor::safe_buffer<BufferPacket<T>> *buff=nullptr;
    convey_t* conv=nullptr;
    BufferPacket<T> done_mark;
    hclib::promise_t<int> worker_loop_end;
    bool is_early_exit = false, is_done = false;
    Mailbox* dep_mb = nullptr;

  public:

    Mailbox() {
        buff = new hclib::conveyor::safe_buffer<BufferPacket<T>>(SIZE);
    }

    ~Mailbox() {
    }

    std::function<void (T, int)> process;

    hclib::future_t<int>* get_worker_loop_finish() {
        return worker_loop_end.get_future();
    }

    hclib::conveyor::safe_buffer<BufferPacket<T>>* get_buffer() {
        return buff;
    }

    void set_is_early_exit(bool val) {
        is_early_exit = val;
    }

    void set_dep_mb(Mailbox* val) {
        dep_mb = val;
    }

    Mailbox* get_dep_mb() {
        return dep_mb;
    }

    void start() {
        conv = convey_new(SIZE_MAX, 0, NULL, 0);
        assert( conv != nullptr );
        convey_begin(conv, sizeof(T), 0);
        done_mark.rank = DONE_MARK;
    }

    void end() {
        delete buff;
        convey_free(conv);
    }

    bool send(T pkt, int rank) {
        if(buff->full()) {
            if(is_early_exit)
                return false;
            else {
                while(buff->full()) { hclib::yield_at(nic);
                //fprintf(stderr, "buff is marked !full\n");
                }
            }
        }
        assert(!buff->full());
        buff->push_back(BufferPacket<T>(pkt, rank));
        return true;
    }

    void done() {
        is_done = true;
        while(buff->full()) hclib::yield_at(nic);
        assert(!buff->full());
        buff->push_back(done_mark);
    }

    int start_worker_loop(int status=0) {

        assert(status == 0);
        hclib::async_at([=]{
        while(true) {
            size_t buff_size = buff->size();
            if(buff_size > 0) break;
            hclib::yield_at(nic);
        }

        BufferPacket<T> bp;
        if(buff->size() > 0)
        bp = buff->at(0);

        //Assumes once 'advance' is called with done=true, the conveyor
        //enters endgame and subsequent value of 'done' is ignored
        while(convey_advance(conv, bp.rank == DONE_MARK)) {
            int i;
            size_t buff_size = buff->size();
            for(i=0;i<buff_size; i++){
                bp = buff->operator[](i);
                if( bp.rank == DONE_MARK) break;
                if( !convey_push(conv, &(bp.data), bp.rank)) {
                    fprintf(stderr, "no-push\n");
                    break;
                }
            }
            if(i>0)
            {
#ifdef USE_LOCK
            std::lock_guard<std::mutex> lg(buff->get_mutex());
#endif
                buff->erase_begin(i);
              }
            int64_t from;
            T pop;
              //while(!get_dep_mb()->get_buffer()->full() &&  convey_pull(conv, &pop, &from) == convey_OK) {
            
            while( convey_pull(conv, &pop, &from) == convey_OK) {
                process(pop, from);
            }
            hclib::yield_at(nic);
          }
          worker_loop_end.put(1);
        }, nic);
        return 0;
    }
};

template<int N, typename T=int64_t, int SIZE=BUFFER_SIZE>
class Selector {

  private:
    void start_worker_loop() {
        for(int i=0;i<N;i++) {
            mb[i].start_worker_loop();
        }
    }

    hclib::promise_t<int> end_prom;
    int num_work_loop_end = 0;

  protected:

  public:

    Mailbox<T, SIZE> mb[N];

    Selector(bool is_start = false) {
        if(is_start) {
            start();
        }
    }

    ~Selector() {
        for(int i=0; i<N; i++) {
            mb[i].end();
        }
    }

    void start() {
        for(int i=0; i<N; i++) {
            //mb[i].set_dep_mb(&mb[(i+1)%N]);
            mb[i].start();
        }
        start_worker_loop();
    }

    bool send(int mb_id, T pkt, int rank) {
        return mb[mb_id].send(pkt, rank);
    }

    bool send(T pkt, int rank) {
        assert(N==1);
        return send(0, pkt, rank);
    }

    void done(int mb_id) {
        mb[mb_id].done();
        hclib::async_await_at([=]() {
            num_work_loop_end++;
            if(num_work_loop_end < N) {
                done((mb_id+1)%N);
            }
            else {
                assert(num_work_loop_end == N);
                end_prom.put(1);
            }
        }, mb[mb_id].get_worker_loop_finish(), nic);
    }

    void done() {
        assert(N==1);
        done(0);
    }

    hclib::future_t<int>* get_future() {
        return end_prom.get_future();
    }
};

template<typename T=int64_t>
using Actor = Selector<1,T>;

}; // namespace hclib

#endif
