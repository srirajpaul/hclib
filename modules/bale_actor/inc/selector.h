
#ifndef SELECTOR_H
#define SELECTOR_H

#include<iostream>
#include<vector>
#include "safe_buffer.h"
#include "hclib_bale_actor.h"
extern "C" {
#include "convey.h"
}

#ifndef NO_USE_BUFFER
#define USE_BUFFER
#endif

// #define ENABLE_TRACE
#define DONE_MARK -1
#define BUFFER_SIZE 1024
#ifndef ELASTIC_BUFFER_SIZE
#define ELASTIC_BUFFER_SIZE 128
#endif

namespace hclib {
//PEtoNodeMap will be used for tracing purpose when ENABLE_TRACE is defined
#ifdef ENABLE_TRACE
#include<papi.h>

#ifndef USE_SHMEM
#error "Trace generation can be used only when USE_SHMEM=1"
#else
#include "shmem.h"
#endif

FILE *get_trace_fptr(bool is_new, const char name[] = "") {
    static FILE *fptr = NULL;
    if (fptr == NULL || is_new) {
        char fname[256];
	if (name != NULL && name[0] != '\0')
	  snprintf(fname, sizeof(fname), "PE%d_%s_send.csv", shmem_my_pe(), name);
	else
	  snprintf(fname, sizeof(fname), "PE%d_send.csv", shmem_my_pe());
        fptr = fopen(fname, "w");
    }
    return fptr;
}

FILE *get_papi_trace_fptr(bool is_new, const char name[] = "") {
    static FILE *fptr = NULL;
    if (fptr == NULL || is_new) {
        char fname[256];
	if (name != NULL && name[0] != '\0')
	  snprintf(fname, sizeof(fname), "PE%d_%s_PAPI.csv", shmem_my_pe(), name);
	else
	  snprintf(fname, sizeof(fname), "PE%d_PAPI.csv", shmem_my_pe());
        fptr = fopen(fname, "w");
    }
    return fptr;
}

/*
  Library function to create a new file for trace output, e.g.:
    char name[32];
    sprintf(name, "phase%d", 3);
    hclib::new_file_for_selector_trace(name);
    // Trace output is written to the file named "PE*_phase3_send.csv" and "PE*_phase3_PAPI.csv".
 */
void new_file_for_selector_trace(char name[]) {
  FILE *fptr1 = get_trace_fptr(true, name);
  FILE *fptr2 = get_papi_trace_fptr(true, name);
}

double get_clock_time() {
    struct timespec tv;
    clock_gettime(CLOCK_REALTIME, &tv);
    double stamp = (double)tv.tv_sec + (double)tv.tv_nsec / 1000000000L;
    return stamp;
}

int *PEtoNodeMap;
void trace_send(int64_t src, int64_t dst, size_t pkg_size) {
    FILE *fptr = get_trace_fptr(false);
    double stamp = get_clock_time();
    fprintf(fptr, "%d, %ld, %d, %ld, %ld, %lf\n", PEtoNodeMap[src], src, PEtoNodeMap[dst], dst, pkg_size, stamp);
}

/*
  Future work: Allows user to specify the PAPI events to be measured.
 */
const char* papi_events[3] = {"PAPI_TOT_CYC", "PAPI_TOT_INS", "PAPI_LST_INS"};
int papi_num_events = 3;
long long papi_tmp_counts[3];

int papi_eventset;
int HCLIB_PAPI_Init() {
    int retval;
    //fprintf(stderr, "[HCLIB PAPI] Initializing PAPI\n");
    retval = PAPI_library_init(PAPI_VER_CURRENT);
    if (retval != PAPI_VER_CURRENT) {
        fprintf(stderr,"Error initializing PAPI! %s\n",
                PAPI_strerror(retval));
        return 0;
    }
    papi_eventset = PAPI_NULL;

    retval = PAPI_create_eventset(&papi_eventset);
    if (retval != PAPI_OK) {
        fprintf(stderr,"Error creating eventset! %s\n",
                PAPI_strerror(retval));
    }
    for (int i = 0; i < papi_num_events; i++) {
        //fprintf(stderr, "[HCLIB PAPI] adding %s \n", papi_events[i]);
        retval = PAPI_add_named_event(papi_eventset, papi_events[i]);
        if (retval != PAPI_OK) {
            fprintf(stderr,"Error adding %s: %s\n",
                    papi_events[i],
                    PAPI_strerror(retval));
        }
    }
    return 0;
}

void inline HCLIB_PAPI_Start() {
    int retval;
    PAPI_reset(papi_eventset);
    retval = PAPI_start(papi_eventset);
    if (retval != PAPI_OK) {
        fprintf(stderr,"Error starting PAPI: %s\n",
                PAPI_strerror(retval));
    }
}

void inline HCLIB_PAPI_Stop(long long *counters, bool accum = true) {
    int retval;
    retval = PAPI_stop(papi_eventset, papi_tmp_counts);
    if (retval != PAPI_OK) {
        fprintf(stderr,"Error stopping:  %s\n",
                PAPI_strerror(retval));
    }
    if (accum) {
        for (int i = 0; i < papi_num_events; i++)
            counters[i] += papi_tmp_counts[i];
    } else {
        for (int i = 0; i < papi_num_events; i++)
            counters[i] = papi_tmp_counts[i];
    }
}

void HCLIB_PAPI_Show(int64_t src, size_t pkg_size, int mbox_id, long long *counters, int num_sends) {
    FILE *fp = get_papi_trace_fptr(false);
    int64_t dst = shmem_my_pe();
    int nd_s = PEtoNodeMap[src];
    int nd_d = PEtoNodeMap[dst];
    long long tot_ins = counters[1];
    long long lst_ins = counters[2];
    double stamp = get_clock_time();
    fprintf(fp, "%d, %ld, %d, %ld, %ld, %d, %d, %lld, %lld, %lf\n",
	    nd_s, src, nd_d, dst, pkg_size, mbox_id, num_sends, tot_ins, lst_ins, stamp);
}

class PapiCounter {
 public:
  int num_sends;
  long long *papi_counters;

  PapiCounter() {
    num_sends = 0;
    papi_counters = (long long *) calloc(papi_num_events, sizeof(long long));
  }

  void init() {
    num_sends = 0;
    for (int i = 0; i < papi_num_events; i++)
      papi_counters[i] = 0;
  }
};

// #define PAPI_TRACER_DEBUG1
// #define PAPI_TRACER_DEBUG2
class PapiTracer {
  std::vector<PapiCounter> counters;
  std::vector<int> reuse_idcs;
  PapiCounter *curr;
  int curr_idx;

 public:
  void init() {
    counters.clear();
    reuse_idcs.clear();
    curr = NULL;
    curr_idx = -1;
#ifdef PAPI_TRACER_DEBUG1
    printf("PapiTracer: PE%d init\n", shmem_my_pe());
#endif
    HCLIB_PAPI_Init();
  }

  void start() {
    assert(curr_idx == -1);
    bool to_add = reuse_idcs.empty();
    if (to_add) {
      counters.push_back(PapiCounter());
      curr_idx = counters.size() - 1;
    } else {
      curr_idx = reuse_idcs.back();
      reuse_idcs.pop_back();
      counters[curr_idx].init();
    }
    curr = &(counters[curr_idx]);
#ifdef PAPI_TRACER_DEBUG2
    const char *msg = to_add ? "added" : "reused";
    printf("PapiTracer: PE%d start %d (counter %s)\n", shmem_my_pe(), curr_idx, msg);
#endif
    HCLIB_PAPI_Start();
  }

  int pause() {
    HCLIB_PAPI_Stop(curr->papi_counters);
    assert(curr_idx >= 0);
#ifdef PAPI_TRACER_DEBUG2
    printf("PapiTracer: PE%d pause %d\n", shmem_my_pe(), curr_idx);
#endif
    int paused = curr_idx;
    curr_idx = -1;
    curr = NULL;
    return paused;
  }

  void resume(int paused) {
    assert(curr_idx == -1);
    curr_idx = paused;
    curr = &(counters[curr_idx]);
    (curr->num_sends)++;
#ifdef PAPI_TRACER_DEBUG2
    printf("PapiTracer: PE%d resume %d\n", shmem_my_pe(), curr_idx);
#endif
    HCLIB_PAPI_Start();
  }

  void end_and_dump(int64_t src, size_t pkg_size, int mb_id) {
    HCLIB_PAPI_Stop(curr->papi_counters);
    assert(curr_idx >= 0);
    HCLIB_PAPI_Show(src, pkg_size, mb_id, curr->papi_counters, curr->num_sends);
#ifdef PAPI_TRACER_DEBUG2
    printf("PapiTracer: PE%d end %d\n", shmem_my_pe(), curr_idx);
#endif
    reuse_idcs.push_back(curr_idx);
    curr_idx = -1;
    curr = NULL;
  }
};
#endif

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
    int mb_id;
#ifdef ENABLE_TRACE
    PapiTracer *papi_tracer = NULL;
#endif

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

#ifdef ENABLE_TRACE
    void start(int mid, PapiTracer *ptrc) {
      papi_tracer = ptrc;
#else
    void start(int mid) {
#endif
        mb_id = mid;
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

#ifdef ENABLE_TRACE
        trace_send(shmem_my_pe(), rank, sizeof(T));
#endif 

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
#ifdef ENABLE_TRACE
                  papi_tracer->start();
#endif
                  process(pop, from);
#ifdef ENABLE_TRACE
                  papi_tracer->end_and_dump(from, sizeof(T), mb_id);
#endif
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
#ifdef ENABLE_TRACE
                  papi_tracer->start();
#endif
                  process(pop, from);
#ifdef ENABLE_TRACE
                  papi_tracer->end_and_dump(from, sizeof(T), mb_id);
#endif
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
#ifdef ENABLE_TRACE
    void createPEtoNodeMap() {
        PEtoNodeMap = (int*)shmem_malloc(shmem_n_pes()*sizeof(int));
        if(PEtoNodeMap==NULL){
            std::cout << "ERROR: Unable to allocate space for PEtoNodeMap pointer\n" << std::endl;
            abort();
        }
        char *nid = getenv("SLURM_NODEID");
        if(nid==NULL){
            std::cout << "ERROR: Unable to retrieve NodeID\n" << std::endl;
            abort();
        }
        int myNodeID = atoi(nid);
        for (int i = 0; i < shmem_n_pes(); i++) {
            shmem_put32(PEtoNodeMap+shmem_my_pe(), &myNodeID, 1, i);
        }
        shmem_barrier_all();
        if(myNodeID==0 && shmem_my_pe()==0){
            printf("Logical actor message trace enabled\n");
            for (int i = 0; i < shmem_n_pes(); i++) {
                printf("PE:%d, Node %d\n", i, PEtoNodeMap[i]);
            }
        }
    }

  PapiTracer papi_tracer;
#endif
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
        #ifdef ENABLE_TRACE
        createPEtoNodeMap();
        #endif
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
#ifdef ENABLE_TRACE
        papi_tracer.init();
#endif
        for(int i=0; i<N; i++) {
            //mb[i].set_dep_mb(&mb[(i+1)%N]);
#ifdef ENABLE_TRACE
            mb[i].start(i, &papi_tracer);
#else
            mb[i].start(i);
#endif
        }
        start_worker_loop();
#ifdef ENABLE_TRACE
        papi_tracer.start();
#endif
    }

#ifdef USE_LAMBDA
    template<typename L>
    bool send(int mb_id, int rank, L lambda) {
#ifdef ENABLE_TRACE
        int idx = papi_tracer.pause();
        bool ret = mb[mb_id].send(rank, lambda);
        papi_tracer.resume(idx);
        return ret;
#else
        return mb[mb_id].send(rank, lambda);
#endif
    }

    template<typename L>
    bool send(int rank, L lambda) {
        assert(N==1);
        return send(0, rank, lambda);
    }
#else
    bool send(int mb_id, T pkt, int rank) {
#ifdef ENABLE_TRACE
        int idx = papi_tracer.pause();
        bool ret = mb[mb_id].send(pkt, rank);
        papi_tracer.resume(idx);
        return ret;
#else
        return mb[mb_id].send(pkt, rank);
#endif
    }

    void yield() {
#ifdef ENABLE_TRACE
        int idx = papi_tracer.pause();
        hclib::yield();
        papi_tracer.resume(idx);
#else
        hclib::yield();
#endif
    }

    bool send(T pkt, int rank) {
        assert(N==1);
        return send(0, pkt, rank);
    }
#endif // USE_LAMBDA

    void done(int mb_id) {
#ifdef ENABLE_TRACE
        if (mb_id == 0)
            papi_tracer.end_and_dump(shmem_my_pe(), 0, -1);
#endif
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
