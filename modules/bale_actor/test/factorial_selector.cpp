/******************************************************************
//
// 
 *****************************************************************/ 
/*! \file dependancy_*_selector.cpp
 * \brief Demo application for graph termination dependancies. 
 *          Calculate the factorial of a randomly generated number from 1 to 10
 *
 */

#include <math.h>
#include <shmem.h>
extern "C" {
#include "spmat.h"
}
#include "selector.h"
#include <cstdlib>
#include <iostream>
using namespace std;

#define THREADS shmem_n_pes()
#define MYTHREAD shmem_my_pe()

enum MailBoxType {A};

typedef struct DepPkt {
    int64_t n;
    int64_t index;
} DepPkt;

class DepSelector: public hclib::Selector<1, DepPkt> {
public:
    DepSelector(int64_t* factorial) : factorial_(factorial) {
        mb[A].process = [this] (DepPkt pkt, int sender_rank) { 
            this->a_process(pkt, sender_rank);
        };
        //mb[A].add_dep_mbs({A}); // successors = {A}
    }

private:
    //shared variables
    int64_t* factorial_;

    void a_process(DepPkt pkg, int sender_rank) {
        if ((pkg.n == 0) || (pkg.n == 1)) {
            *factorial_ *= 1;
            done(A); //done_extended(A);
        } else {*factorial_ *= pkg.n;}
        DepPkt next_pkg;
        next_pkg.n = pkg.n - 1;
        send(A, next_pkg, sender_rank);
    }
};

int main(int argc, char* argv[]) {
    const char *deps[] = { "system", "bale_actor" };

    hclib::launch(deps, 2, [=] {

        //srand(time(0)); 

        T0_fprintf(stderr, "\n\nFactorial of a randomly generated number between 1 and 10 (selector)...\n");

        int64_t num_input = 1+ (rand() % 10);

        int64_t factorial = 1;

        DepSelector* depSelector = new DepSelector(&factorial);

        hclib::finish([=]() {
            depSelector->start();
            DepPkt pkg;
            pkg.n = num_input;
            depSelector->send(A, pkg, MYTHREAD);
        });

        lgp_barrier();

        printf("\n[PE %d] FACTORIAL OF %d IS (%d!): %d", MYTHREAD, num_input, num_input, factorial);
        printf("\n");

        lgp_finalize();

    });

    return 0;
}
