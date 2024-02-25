
#include <shmem.h>
#include <stdio.h>
#include "selector.h"

class TestSelector: public hclib::Selector<1, int64_t> {

    void process(int64_t pkt, int sender_rank) {
        printf("In rank %d, val %d, from rank %d\n", shmem_my_pe(), pkt, sender_rank);
    }

  public:

    TestSelector(): hclib::Selector<1, int64_t>(true) {
        mb[0].process = [this](int64_t pkt, int sender_rank) { this->process(pkt, sender_rank); };
    }
};

int main() {

  const char *deps[] = { "system", "bale_actor" };
  hclib::launch(deps, 2, [=] {

    TestSelector *ts_ptr = new TestSelector();

    int num = 10, dest_rank = (shmem_my_pe() + 1)%shmem_n_pes();
    for(int i=0;i<num;i++) {
        int64_t val = shmem_my_pe() * 1000 + i;
        ts_ptr->send(0, val, dest_rank);
    }
    ts_ptr->done(0); // Indicate that we are done with sending messages to the REQUEST mailbox

    hclib::async_await([=]() {
        printf("Done\n");
    }, ts_ptr->get_future());

  });
  printf("Outside Launch\n");
  
  return 0;
}
