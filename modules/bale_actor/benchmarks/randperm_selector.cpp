/******************************************************************
//
//
//  Copyright(C) 2019, Institute for Defense Analyses
//  4850 Mark Center Drive, Alexandria, VA; 703-845-2500
//  This material may be reproduced by or for the US Government
//  pursuant to the copyright license under the clauses at DFARS
//  252.227-7013 and 252.227-7014.
// 
//
//  All rights reserved.
//  
//  Redistribution and use in source and binary forms, with or without
//  modification, are permitted provided that the following conditions are met:
//    * Redistributions of source code must retain the above copyright
//      notice, this list of conditions and the following disclaimer.
//    * Redistributions in binary form must reproduce the above copyright
//      notice, this list of conditions and the following disclaimer in the
//      documentation and/or other materials provided with the distribution.
//    * Neither the name of the copyright holder nor the
//      names of its contributors may be used to endorse or promote products
//      derived from this software without specific prior written permission.
// 
//  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
//  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
//  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
//  FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
//  COPYRIGHT HOLDER NOR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
//  INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
//  (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
//  SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
//  HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
//  STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
//  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
//  OF THE POSSIBILITY OF SUCH DAMAGE.
// 
 *****************************************************************/ 

#include <shmem.h>
extern "C" {
#include <spmat.h>
}
#include "selector.h"

#define THREADS shmem_n_pes()
#define MYTHREAD shmem_my_pe()

/*! \file randperm_selector.cpp
 * \brief Demo program that runs the variants of randperm kernel. This program
 * generates a random permutation in parallel.
 */

/*! 
\page randperm_page Random Permutation

Demo program that runs the variants of randperm kernel. This program
generates a random permutation in parallel. The algorithm used is the 
"dart throwing algorithm" found in 
 P.B.Gibbon, Y.Matias, and V.L.Ramachandran. Efficient low-contention Parallel algorithms.
 J. of Computer and System Sciences, 53:417-442, Dec 1992.

Interestingly, we discovered what looks to be a superior and simpler algorithm. This is implemented
in alternates/randperm_agi_opt.upc.

See files spmat_agi.upc, spmat_exstack.upc, spmat_exstack2.upc, and spmat_conveyor.upc
for the source for the kernels.

Usage:
randperm [-h][-n num][-M mask][-s seed]
- -h Print usage banner only
- -b count is the number of packages in an exstack(2) buffer
- -n=num Set the permutation entries per PE to n (default = 1000).
- -M=mask Set the models mask (1,2,4,8,16,32 for AGI,exstack,exstack2,conveyor,alternate)
- -s=seed Set a seed for the random number generation.
 */

/*
    Regular Benchmark stuff starts here
*/
typedef struct pkg_t {
    int64_t idx; 
    int64_t val;
} RPpkg;

/*
    In the original code the reply packet is just int64_t,
    since Hclib Selector doesn't support different mailboxes of different
    packet type, the reply packge will have idx field = data, val should 
    always 0
*/

enum PhaseOneMailBoxType {THROW, REPLY};
enum PhaseTwoMailBoxType {MSG};

class RandPermSelectorPhaseOne: public hclib::Selector<2, RPpkg> {
public:
    RandPermSelectorPhaseOne(
        int64_t* ltarget,
        int64_t* lperm,
        int64_t* iend,
        int64_t* hits
    ) : ltarget_(ltarget), lperm_(lperm), iend_(iend), hits_(hits) {
        mb[THROW].process = [this](RPpkg pkg, int senderRank) { this->throwProcess(pkg, senderRank); };
        mb[REPLY].process = [this](RPpkg pkg, int senderRank) { this->replyProcess(pkg, senderRank); };
    }

private:
    int64_t* ltarget_;
    int64_t* lperm_;
    int64_t* iend_;
    int64_t* hits_;

    void throwProcess(RPpkg pkg, int senderRank) {
        if (ltarget_[pkg.idx] == -1L) {
            int64_t val = pkg.val;

            send(REPLY, RPpkg{val, 0}, senderRank);
            ltarget_[pkg.idx] = pkg.val;
        } else {
            int64_t val = -(pkg.val + 1);

            send(REPLY, RPpkg{val, 0}, senderRank);
        }
    }

    void replyProcess(RPpkg pkg, int senderRank) {
        if (pkg.idx < 0L) {
            lperm_[--(*iend_)] = -(pkg.idx) - 1;
        } else {
            (*hits_)++;
        }
    }
};

class RandPermSelectorPhaseTwo: public hclib::Selector<1, int64_t> {
public:
    RandPermSelectorPhaseTwo(int64_t* lperm, int64_t& pos) : lperm_(lperm), pos_(pos) {
        mb[MSG].process = [this](int64_t pkg, int senderRank) { this->processMsg(pkg, senderRank); };
    }
private:
    int64_t* lperm_;
    int64_t& pos_;

    void processMsg(int64_t val, int senderRank) {
        lperm_[pos_++] = val;
    }
};


int64_t* copied_rand_permp_selector(int64_t N, int seed) {
    int64_t lN = (N + THREADS - MYTHREAD - 1) / THREADS;
    int64_t M = N * 2;
    int64_t lM = (M + THREADS - MYTHREAD - 1) / THREADS;

    int64_t* perm = (int64_t*)lgp_all_alloc(N, sizeof(int64_t));
    if (!perm) return nullptr;
    int64_t* lperm = lgp_local_part(int64_t, perm);

    int64_t* target = (int64_t*)lgp_all_alloc(M, sizeof(int64_t));
    if (!target) return nullptr;
    int64_t* ltarget = lgp_local_part(int64_t, target);;
  
    /* initialize perm[i] = i,  the darts*/
    for (int64_t i = 0; i < lN; i++)
        lperm[i] = i * THREADS + MYTHREAD;

    /* initialize target[i] = -1 */
    for (int64_t i = 0; i < lM; i++)
        ltarget[i] = -1L;

    if (seed != 0) srand48(seed);

    int64_t hits = 0;
    int64_t iend = 0;
    int64_t* hitsPtr = &hits;
    int64_t* iendPtr = &iend;

    RandPermSelectorPhaseOne* phaseOneSelector = new RandPermSelectorPhaseOne(ltarget, lperm, iendPtr, hitsPtr);
    
    // setup finish, start algorithm
    lgp_barrier();
    double t1 = wall_seconds();

    hclib::finish([lN, M, lperm, hitsPtr, iendPtr, phaseOneSelector]() {
        phaseOneSelector->mb[THROW].set_is_early_exit(true);
        phaseOneSelector->start();
        pkg_t pkg;
        //int64_t i = 0;

        // Since a throw a fail, we need to keep track of hits
        // instead of believing our lN request will all result in hits
        while (*hitsPtr != lN) {
            //i = *iendPtr;

            while (*iendPtr < lN) {
                int64_t r = lrand48() % M;
                int64_t pe = r % THREADS;
                pkg.idx = r / THREADS;
                pkg.val = lperm[*iendPtr];

                bool ret = phaseOneSelector->send(THROW, pkg, pe);
                //i++;

                if(ret)
                    (*iendPtr)++;
                else
                    hclib::yield();

            }

            //*iendPtr = i;

            // let the mailbox process in order for hits to update
            hclib::yield();

            // If enough hits are processed, break and teardown
            if (*hitsPtr >= lN) { break; }
        }
        phaseOneSelector->done(THROW);
    });

    lgp_barrier();
    t1 = wall_seconds() - t1;
    delete phaseOneSelector;

    // T0_printf("phase 1 t1 = %8.3lf\n", t1);
    // T0_printf("DEBUG: phase 1 finished, no task starving\n");

    /* now locally pack the values you have in target */
    int64_t cnt = 0;
    for (int64_t i = 0; i < lM; i++) {
        if (ltarget[i] != -1L) {
            ltarget[cnt++] = ltarget[i];
        }
    }
    lgp_barrier();
  
    /* sanity check */
    int64_t total = lgp_reduce_add_l(cnt);
    if (total != N) {
        T0_printf("ERROR: rand_permp_selector: total = %ld should be %ld\n", total, N);

        return nullptr;
    }

    int64_t offset = lgp_prior_add_l(cnt);
    int64_t pos = 0;
    RandPermSelectorPhaseTwo* phaseTwoSelector = new RandPermSelectorPhaseTwo(lperm, pos);

    hclib::finish([phaseTwoSelector, cnt, offset, ltarget]() {
        phaseTwoSelector->start();
        int64_t i = 0;
        int64_t pe = offset % THREADS;

        while (i < cnt) {
            phaseTwoSelector->send(MSG, ltarget[i], pe);

            i++;
            pe++;

            if (pe == THREADS) pe = 0;
        }

        phaseTwoSelector->done(MSG);
    });

    pos = lgp_reduce_add_l(pos);
    if (pos != N) {
        printf("ERROR! in rand_permp_selector! sum of pos = %ld lN = %ld\n", pos, N);

        return nullptr;
    }

    lgp_barrier();
    delete phaseTwoSelector;

    return perm;
}

int main(int argc, char* argv[]) {
    const char *deps[] = { "system", "bale_actor" };
    hclib::launch(deps, 2, [=] {
        int64_t i;
        int64_t models_mask = 0xF;
        int printhelp = 0;
        int64_t l_numrows = 1000000;
        int64_t numrows;
        int64_t buf_cnt = 1024;
        int64_t seed = 101892 + MYTHREAD;
        int64_t cores_per_node = 1;

        int opt; 
        while ((opt = getopt(argc, argv, "c:hn:M:s:")) != -1 ){
            switch (opt) {
                case 'h': printhelp = 1; break;
                case 'b': sscanf(optarg,"%ld",&buf_cnt);  break;
                case 'c': sscanf(optarg,"%ld" ,&cores_per_node); break;
                case 'n': sscanf(optarg,"%ld",&l_numrows);   break;
                case 'M': sscanf(optarg,"%ld",&models_mask);  break;
                case 's': sscanf(optarg,"%ld", &seed); break;      
                default:  break;
            }
        }

        T0_fprintf(stderr, "Running randperm on %d PEs\n", THREADS);
        T0_fprintf(stderr, "This is a demo program that runs various implementations of the randperm kernel.\n");
        T0_fprintf(stderr, "Usage:\n");
        T0_fprintf(stderr, "Permutation size per PE (-n) = %ld\n", l_numrows);
        T0_fprintf(stderr, "models_mask (-M)                 = %ld or of 1,2,4,8 for atomics,classic,exstack2,conveyor\n", models_mask);
        T0_fprintf(stderr, "seed (-s)                        = %ld\n", seed);

        numrows = l_numrows * THREADS;

        double t1;
        minavgmaxD_t stat[1];
        int64_t error = 0;
        int64_t* out;

        int64_t use_model;
  
        T0_fprintf(stderr, "Running rand_permp_selector\n");
        t1 = wall_seconds();
        out = copied_rand_permp_selector(numrows, seed);
        t1 = wall_seconds() - t1;
        
        T0_fprintf(stderr, "rand_permp_selector:           \n");
        lgp_min_avg_max_d(stat, t1, THREADS);
        T0_fprintf(stderr, " %8.3lf seconds\n", stat->avg);

#if USE_ERROR_CHECK
        if (!is_perm(out, numrows)) {
            error++;
            T0_printf("\nERROR: rand_permp_selector failed!\n\n");
        }
        lgp_all_free(out);
  
        if (error) {
            T0_fprintf(stderr,"BALE FAIL!!!!\n"); 
        }
#endif // USE_ERROR_CHECK

    });

    return 0;
}

