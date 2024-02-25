
#ifndef YIELD_LOOP
#define YIELD_LOOP
#endif

#include <shmem.h>
#include <stdio.h>
#include "hclib_bale_actor.h"
#include "selector.h"

class TestSelector: public hclib::Selector<3, int64_t> {
    int64_t* counterPtr_;

    void process0(int64_t pkt, int sender_rank) {
        printf("Process0 In rank %d, val %d, from rank %d\n", shmem_my_pe(), pkt, sender_rank);
        send(1, pkt, sender_rank);
    }

    void process1(int64_t pkt, int sender_rank) {
        printf("Process1 In rank %d, val %d, from rank %d\n", shmem_my_pe(), pkt, sender_rank);
        send(2, pkt, sender_rank);
    }

    void process2(int64_t pkt, int sender_rank) {
        printf("Process2 In rank %d, val %d, from rank %d\n", shmem_my_pe(), pkt, sender_rank);
        (*counterPtr_)++;
    }

  public:

    TestSelector(int64_t* counterPtr) : counterPtr_(counterPtr) {
        mb[0].process = [this](int64_t pkt, int sender_rank) { this->process0(pkt, sender_rank); };
        mb[1].process = [this](int64_t pkt, int sender_rank) { this->process1(pkt, sender_rank); };
        mb[2].process = [this](int64_t pkt, int sender_rank) { this->process2(pkt, sender_rank); };
    }
};

int main() {

  const char *deps[] = { "system", "bale_actor" };
  hclib::launch(deps, 2, [=] {
    int64_t counter = 0;
    int64_t* counterPtr = &counter;

    TestSelector *ts_ptr = new TestSelector(counterPtr);

    printf("Finished inititialization\n");
    hclib::finish([=]() {
      ts_ptr->start();
      //while (*counterPtr < 10) {
        printf("DEBUG: counterPtr has value %lu \n", *counterPtr);
        int num = 10, dest_rank = (shmem_my_pe() + 1)%shmem_n_pes();
        for(int i=0;i<num;i++) {
            int64_t val = shmem_my_pe() * 1000 + i;
            ts_ptr->send(0, val, dest_rank);
        }

        //hclib::yield();
        //printf("DEBUG: after yield counterPtr has value %lu \n", *counterPtr);
      //}
      printf("DEBUG: ready to teardown\n");

      ts_ptr->done(0); // Indicate that we are done with sending messages to the REQUEST mailbox
    });
    printf("Outside Finish\n");
  });
  printf("Outside Launch\n");
  
  return 0;
}
