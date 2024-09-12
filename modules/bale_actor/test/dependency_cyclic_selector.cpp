/******************************************************************
//
// 
 *****************************************************************/ 
/*! \file dependancy_cylic_selector.cpp
 * \brief Demo application for graph termination dependancies (with complex cyclic/acyclic dependencies and cycles). 
 *
 */

#include <math.h>
#include <shmem.h>
extern "C" {
#include "spmat.h"
}
#include "selector.h"

#define THREADS shmem_n_pes()
#define MYTHREAD shmem_my_pe()

enum MailBoxType {A, B, C, D, E};

typedef struct DepPkt {
    int64_t i;
} DepPkt;

class DepSelector: public hclib::Selector<5, DepPkt> {
public:
    DepSelector(int64_t* sum_arr) : sum_arr_(sum_arr) {
        mb[A].process = [this] (DepPkt pkt, int sender_rank) { 
            this->a_process(pkt, sender_rank);
        }; 
        mb[A].add_dep_mbs({B,D}); // successors = {B,D}

        mb[B].process = [this] (DepPkt pkt, int sender_rank) { 
            this->b_process(pkt, sender_rank);
        };
        mb[B].add_dep_mbs({}); // successors = {}
        
        mb[C].process = [this] (DepPkt pkt, int sender_rank) { 
            this->c_process(pkt, sender_rank);
        };
        mb[C].add_dep_mbs({D,E}); // successors = {D,E}

        mb[D].process = [this] (DepPkt pkt, int sender_rank) { 
            this->d_process(pkt, sender_rank);
        };
        mb[D].add_dep_mbs({}); // successors = {}

        mb[E].process = [this] (DepPkt pkt, int sender_rank) { 
            this->e_process(pkt, sender_rank);
        };
        mb[E].add_dep_mbs({}); // successors = {}
    }

private:
    //shared variables
    int64_t* sum_arr_;
    int counter = 0;

    void a_process(DepPkt pkg, int sender_rank) {
        sum_arr_[0] = counter;
        counter++;
        send(B, pkg, sender_rank);
        printf("sending message to mailbox B ...\n");
        send(D, pkg, sender_rank);
        printf("sending message to mailbox D ...\n");
        send(A, pkg, sender_rank);
        printf("sending message to mailbox A ...\n"); //cyclic
        if (counter == 5) {
            done_extended(A);
        }
    }

    void b_process(DepPkt pkg, int sender_rank) {
        sum_arr_[1] = MYTHREAD;
    }
    
    void c_process(DepPkt pkg, int sender_rank) {
        sum_arr_[2] = MYTHREAD;
        send(D, pkg, sender_rank);
        printf("sending message to mailbox D ...\n");
        send(E, pkg, sender_rank);
        printf("sending message to mailbox E ...\n");
    }

    void d_process(DepPkt pkg, int sender_rank) {
        sum_arr_[3] = MYTHREAD;
    }

    void e_process(DepPkt pkg, int sender_rank) {
        sum_arr_[4] = MYTHREAD;
    }
    
};

int main(int argc, char* argv[]) {
    const char *deps[] = { "system", "bale_actor" };

    hclib::launch(deps, 2, [=] {

        int64_t* sum_arr = (int64_t*)lgp_all_alloc(5, sizeof(int64_t));
        for (int64_t x = 0; x < 5; x++) sum_arr[x] = 0;

        printf("\n\n");
        printf("Initial array contents: ");
        for (int64_t x = 0; x < 5; x++) printf("%d ", sum_arr[x]);
        printf("\n                        A B C D E");
        printf("\n\n");

        DepSelector* depSelector = new DepSelector(sum_arr);

        hclib::finish([=]() {
            depSelector->start();
            DepPkt pkg;
            pkg.i = 0;

            /* ---------------------------- */
            /*      DEPENDANCY GRAPH:       */
            //           ----
            //           ˅  |
            //          A   -->   B
            //              \
            //               \
            //                >
            //          C   -->   D
            //              \
            //               \
            //                 >  E
            /* ---------------------------- */

            // send to mailbox A        (A -> B)
            printf("sending message to mailbox A ...\n");
            depSelector->send(A, pkg, MYTHREAD);

            // send to mailbox C        (C -> D       C -> E)
            printf("sending message to mailbox C ...\n");
            depSelector->send(C, pkg, MYTHREAD);

            // invoke done on mailbox A and C
            depSelector->done_extended(C);
        });

        printf("\n\n");
        printf("Final array contents: ");
        for (int64_t x = 0; x < 5; x++) printf("%d ", sum_arr[x]);
        printf("\n                      A B C D E");
        printf("\n\n");

        lgp_barrier();

        lgp_finalize();

    });

    return 0;
}
