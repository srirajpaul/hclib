
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

#ifdef USE_LAMBDA
class BaseLambdaPacket {
  public:
    virtual void invoke() = 0;
    virtual size_t get_bytes() = 0;
    virtual ~BaseLambdaPacket() {}
};

template<typename L>
class LambdaPacket : public BaseLambdaPacket {
    L lambda;

  public:
    LambdaPacket(L lambda) : lambda(lambda) {}

    void invoke() {
        lambda();
    }

    size_t get_bytes() {
        return sizeof(*this);
    }
};

template<typename T>
struct BufferPacket {
    T data;
    int64_t rank;
    BaseLambdaPacket* lambda_pkt;

    BufferPacket() : rank(0) {}
    BufferPacket(T data, int64_t rank) : data(data), rank(rank) {}
    BufferPacket(int64_t rank, BaseLambdaPacket* lambda_pkt) : rank(rank), lambda_pkt(lambda_pkt) {}
};

#else // USE_LAMBDA

template<typename T>
struct BufferPacket {
    T data;
    int64_t rank;

    BufferPacket() : rank(0) {}
    BufferPacket(T data, int64_t rank) : data(data), rank(rank) {}
};

#endif // USE_LAMBDA

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
#ifdef SELECTOR_DEBUG
        printf("Creating Mailbox\n");
#endif
    }

    ~Mailbox() {
#ifdef SELECTOR_DEBUG
        printf("Deleting Mailbox\n");
#endif
        //delete buff;
        //convey_free(conv);
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
#ifdef USE_LAMBDA
        conv = convey_new_elastic(ELASTIC_BUFFER_SIZE, SIZE_MAX, 0, NULL, convey_opt_PROGRESS);
#else
        //buff = new hclib::conveyor::safe_buffer<BufferPacket<T>>(SIZE);
        //conv = convey_new(SIZE_MAX, 0, NULL, 0);
        conv = convey_new(SIZE_MAX, 0, NULL, convey_opt_PROGRESS);
#endif
        assert( conv != nullptr );
        convey_begin(conv, sizeof(T), 0);
        done_mark.rank = DONE_MARK;
    }

    void end() {
        delete buff;
        convey_free(conv);
    }

#ifdef USE_LAMBDA
    template<typename L>
    bool send(int rank, L lambda) {
        //printf("size %d\n", sizeof(lambda));
        if(buff->full()) {
            if(is_early_exit)
                return false;
            else
                while(buff->full()) hclib::yield_at(nic);
        }
        assert(!buff->full());
        buff->push_back(BufferPacket<T>(rank, new LambdaPacket<L>(lambda)));
        return true;
    }
#else
    bool send(T pkt, int rank) {
#ifdef USE_BUFFER
        if(buff->full()) {
            if(is_early_exit)
                return false;
            else
                while(buff->full()) hclib::yield_at(nic);
        }
        assert(!buff->full());
        buff->push_back(BufferPacket<T>(pkt, rank));
        return true;
#else // USE_BUFFER
        int ret = convey_push(conv, &pkt, rank);
        if(is_early_exit)
            return ret == convey_OK;
        else if(ret != convey_OK)
            while(convey_push(conv, &pkt, rank) != convey_OK) hclib::yield_at(nic);
        return true;
#endif // USE_BUFFER
    }
#endif // USE_LAMBDA

    void done() {
        is_done = true;
        while(buff->full()) hclib::yield_at(nic);
        assert(!buff->full());
        buff->push_back(done_mark);
    }


#ifndef YIELD_LOOP
    int start_worker_loop(int status=0) {

        assert(status == 0);
        hclib::async_at([=]{
#ifdef USE_BUFFER
        
        /* PATCH 0.2023 (by Shubh S). Because of this case, we realise that buf.size() = 0 
        assumes the fact : if there is "nothing to push", there would be "nothing to pull".
        This assumption is violated specifically in the graph algorithms where, (for a PE, and a mailbox)
        
        let's say the node has incoming edges only, and
        the decision on arriving at that node by some algorithm A follows no push but only pulls from 
        all the neighbours in order to continue and finish. 
        
        In such case, the program would hang with the code below, as it would yield(), and
        since there is nothing to push, it would yield from any other operation as well.

        Types of Transactions:

        send() - if partial, return true T_0. Else full yield(in a while). T_1
        starter_worker_loop - if buf.size() == 0, yield T_2, else convey advance (push, pull) and then yield T_3

        T_0 would complete and continue working on some other tasks. T_1 would yield because of our assumption of no push, only pull and buf.size being 0.
        T_2 would yield, and T_3 would not enter ever provided T_2 exists!

        Problem: Because T_0 is the only task running in this case, program would hang!   
        Solution: Allow pulls even when buf.size() = 0, as after pulls, yield occurs!

        Note: This is not a performance optimization but a bug solved which resulted from Bottom BFS experimentation with 2 PE's and scale = 16.

         while(true) {
                size_t buff_size = buff->size();
                if(buff_size > 0) break;
                hclib::yield_at(nic);
         } 
        */

          BufferPacket<T> bp;
          if(buff->size() > 0)
            bp = buff->at(0);

          //Assumes once 'advance' is called with done=true, the conveyor
          //enters endgame and subsequent value of 'done' is ignored
          while(convey_advance(conv, bp.rank == DONE_MARK)) {
            // if(bp.rank == DONE_MARK) {
            //     fprintf(stderr, "endgame begin, my_pe: %d\n", shmem_my_pe());
            //   }
              int i;
              size_t buff_size = buff->size();
              for(i=0;i<buff_size; i++){
                  bp = buff->operator[](i);
                  if( bp.rank == DONE_MARK) break;
#ifdef USE_LAMBDA
                  //printf("size %d\n", bp.lambda_pkt->get_bytes());
                  if( !convey_epush(conv, bp.lambda_pkt->get_bytes(), bp.lambda_pkt, bp.rank)) break;
                  //delete bp.lambda_pkt;
#else
                  if( !convey_push(conv, &(bp.data), bp.rank)) break;
#endif //  USE_LAMBDA
              }

	          if(i>0)
              {
#ifdef USE_LOCK
                  std::lock_guard<std::mutex> lg(buff->get_mutex());
#endif
                  buff->erase_begin(i);
              }
#else // USE_BUFFER
          while(convey_advance(conv, is_done)) {
#endif // USE_BUFFER
              int64_t from;
#ifdef USE_LAMBDA
              convey_item_t item;
              while( convey_epull(conv, &item) == convey_OK) {
                  BaseLambdaPacket* data = (BaseLambdaPacket*)item.data;
                  data->invoke();
              }
#else

              T pop;
              //while(!get_dep_mb()->get_buffer()->full() &&  convey_pull(conv, &pop, &from) == convey_OK) {
              while( convey_pull(conv, &pop, &from) == convey_OK) {
                  //hclib::async([=]() { process(pop, from); });
                  process(pop, from);
              }
#endif // USE_LAMBDA

              hclib::yield_at(nic);
          }
          worker_loop_end.put(1);
        }, nic);
        return 0;
    }

#else // YIELD_LOOP

    int start_worker_loop(int status=0) {

          while(true) {
              size_t buff_size = buff->size();
              if(buff_size > 0) break;
              if(status == 1)
                  return 1;
              else {
                  assert(status == 2);
                  break;
              }
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
#ifdef USE_LAMBDA
                  //printf("size %d\n", bp.lambda_pkt->get_bytes());
                  if( !convey_epush(conv, bp.lambda_pkt->get_bytes(), bp.lambda_pkt, bp.rank)) break;
                  //delete bp.lambda_pkt;
#else
                  if( !convey_push(conv, &(bp.data), bp.rank)) break;
#endif //  USE_LAMBDA
              }

	          if(i>0)
              {
#ifdef USE_LOCK
                  std::lock_guard<std::mutex> lg(buff->get_mutex());
#endif
                  buff->erase_begin(i);
              }
              int64_t from;
#ifdef USE_LAMBDA
              convey_item_t item;
              while( convey_epull(conv, &item) == convey_OK) {
                  BaseLambdaPacket* data = (BaseLambdaPacket*)item.data;
                  data->invoke();
              }
#else

              T pop;
              //while(!get_dep_mb()->get_buffer()->full() &&  convey_pull(conv, &pop, &from) == convey_OK) {
              while( convey_pull(conv, &pop, &from) == convey_OK) {
                  //hclib::async([=]() { process(pop, from); });
                  process(pop, from);
              }
#endif // USE_LAMBDA

              return 2;
          }
          worker_loop_end.put(1);
          return 0;
    }

#endif // YIELD_LOOP

};

template<int N, typename T=int64_t, int SIZE=BUFFER_SIZE>
class Selector {

  private:
#ifndef YIELD_LOOP
    void start_worker_loop() {
        for(int i=0;i<N;i++) {
            mb[i].start_worker_loop();
        }
    }
#else
    void start_worker_loop() {
        hclib::async_at([=]{
            int loop_stat[N];
            std::fill_n(loop_stat, N, 1);
            int finish_count = 0;

            while(finish_count < N) {
              for(int i=0;i<N;i++) {
                if(loop_stat[i] != 0) {
                  loop_stat[i] = mb[i].start_worker_loop(loop_stat[i]);
                  if(loop_stat[i] == 0)
                    finish_count++;
                }
              }
              hclib::yield_at(nic);
            }
        }, nic);
    }
#endif

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

#ifdef USE_LAMBDA
    template<typename L>
    bool send(int mb_id, int rank, L lambda) {
        return mb[mb_id].send(rank, lambda);
    }

    template<typename L>
    bool send(int rank, L lambda) {
        assert(N==1);
        return send(0, rank, lambda);
    }
#else
    bool send(int mb_id, T pkt, int rank) {
        return mb[mb_id].send(pkt, rank);
    }

    bool send(T pkt, int rank) {
        assert(N==1);
        return send(0, pkt, rank);
    }
#endif // USE_LAMBDA

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