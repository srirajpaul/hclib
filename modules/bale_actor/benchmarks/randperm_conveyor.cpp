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
#include "shmem.h"
extern "C" {
#include <spmat.h>
}

/*! \file randperm_conveyor.cpp
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

int64_t* copied_rand_permp_conveyor(int64_t N, int seed) {
    int ret;
    int64_t i, j, cnt, pe, pos, fromth, istart, iend;
    int64_t val;
  
    int64_t lN = (N + THREADS - MYTHREAD - 1) / THREADS;
    int64_t M = N * 2;
    int64_t lM = (M + THREADS - MYTHREAD - 1) / THREADS;

    typedef struct pkg_t {
        int64_t idx; 
        int64_t val;
    } pkg_t;

    int64_t* perm = (int64_t*)lgp_all_alloc(N, sizeof(int64_t));
    if (!perm) return nullptr;
    int64_t* lperm = lgp_local_part(int64_t, perm);

    int64_t* target = (int64_t*)lgp_all_alloc(M, sizeof(int64_t));
    if (!target) return nullptr;
    int64_t* ltarget = lgp_local_part(int64_t, target);
  
    /* initialize perm[i] = i,  the darts*/
    for (i = 0; i < lN; i++)
        lperm[i] = i * THREADS + MYTHREAD;

    /* initialize target[i] = -1 */
    for (i = 0; i < lM; i++)
        ltarget[i] = -1L;

    if (seed != 0) srand48(seed);

    lgp_barrier();
  
    double t1 = wall_seconds();
    int64_t rejects = 0;
    pkg_t pkg;
    convey_t* conv_throw = convey_new(SIZE_MAX, 0, NULL, convey_opt_PROGRESS);
    convey_t* conv_reply = convey_new(SIZE_MAX, 0, NULL, convey_opt_PROGRESS);
    if (!conv_throw) { return nullptr; }
    if (!conv_reply) { return nullptr; }

    convey_begin(conv_throw, sizeof(pkg_t), 0);
    convey_begin(conv_reply, sizeof(int64_t), 0);
  
    bool more;
    int64_t hits = 0;
    iend = 0;
    while (more = convey_advance(conv_throw, (hits == lN)), more | convey_advance(conv_reply, !more)) {
        i = iend;
        while (i < lN) {
            int64_t r = lrand48() % M;
            pe = r % THREADS;
            pkg.idx = r / THREADS;
            pkg.val = lperm[i];
            
            if (!convey_push(conv_throw, &pkg, pe)) { break; }
            i++;
        }

        iend = i;
    
        while (convey_pull(conv_throw, &pkg, &fromth) == convey_OK) {
            if (ltarget[pkg.idx] == -1L) {
                val = pkg.val;
                if (!convey_push(conv_reply, &val, fromth)) {
                    convey_unpull(conv_throw);
                    break;
                }

                ltarget[pkg.idx] = pkg.val;
            } else {
                val = -(pkg.val + 1);
                if (!convey_push(conv_reply, &val, fromth)) {
                    convey_unpull(conv_throw);
                    rejects++;
                    break;
                }
            }
        }

        while (convey_pull(conv_reply, &val, NULL) == convey_OK) {
            if (val < 0L) {
                lperm[--iend] = -val - 1;
            } else {
                hits++;
            }
        }
    }
  
    lgp_barrier();
    t1 = wall_seconds() - t1;
    //T0_printf("phase 1 t1 = %8.3lf\n", t1);
    //rejects = lgp_reduce_add_l(rejects);
    //T0_printf("rejects = %ld\n", rejects);

    convey_free(conv_throw);
    convey_reset(conv_reply);
    convey_begin(conv_reply, sizeof(int64_t), 0);
  
    /* now locally pack the values you have in target */
    cnt = 0;
    for (i = 0; i < lM; i++) {
        if (ltarget[i] != -1L) {
            ltarget[cnt++] = ltarget[i];
        }
    }
    lgp_barrier();
  
    /* sanity check */
    int64_t total = lgp_reduce_add_l(cnt);
    if (total != N) {
        T0_printf("ERROR: rand_permp_convey: total = %ld should be %ld\n", total, N);
        return nullptr;
    }

    int64_t offset = lgp_prior_add_l(cnt);
    pe = offset % THREADS;
    i = pos = 0;
    while (convey_advance(conv_reply, (i == cnt))) {
        while (i < cnt) {
            if (!convey_push(conv_reply, &ltarget[i], pe)) break;
        
            i++;
            pe++;

            if (pe == THREADS) pe = 0;
        }

        while (convey_pull(conv_reply, &val, &fromth)) {
            lperm[pos++] = val;
        }
    }

    pos = lgp_reduce_add_l(pos);
    if (pos != N) {
        printf("ERROR! in rand_permp_convey! sum of pos = %ld lN = %ld\n", pos, N);
        return nullptr;
    }
    lgp_barrier();
  
    convey_free(conv_reply);

    return perm;
}

int main(int argc, char* argv[]) {
    lgp_init(argc, argv);
  
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
  
    t1 = wall_seconds();
    out = copied_rand_permp_conveyor(numrows, seed);
    t1 = wall_seconds() - t1;

    T0_fprintf(stderr, "rand_permp_conveyor:           \n");
    lgp_min_avg_max_d(stat, t1, THREADS);
    T0_fprintf(stderr, " %8.3lf seconds\n", stat->avg);

#if USE_ERROR_CHECK
    if (!is_perm(out, numrows)) {
        error++;
        T0_printf("\nERROR: rand_permp_%ld failed!\n\n", use_model & models_mask);
    }
    lgp_all_free(out);
  
    if (error) {
        T0_fprintf(stderr,"BALE FAIL!!!!\n"); 
    }
#endif // USE_ERROR_CHECK

    lgp_finalize();
  
    return error;
}
