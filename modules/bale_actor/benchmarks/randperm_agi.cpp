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

/*! \file randperm_agi.cpp
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

int64_t* copied_rand_permp_agi(int64_t N, int seed) {
    int64_t* ltarget, *lperm;
    int64_t r, i, j;
    int64_t pos, numdarts, numtargets, lnumtargets;

    if (seed != 0) srand48(seed);

    //T0_printf("Entering rand_permp_atomic...");fflush(0);

    int64_t* perm = (int64_t*)lgp_all_alloc(N, sizeof(int64_t));
    if (!perm) return nullptr;
    lperm = lgp_local_part(int64_t, perm);

    int64_t l_N = (N + THREADS - MYTHREAD - 1) / THREADS;
    int64_t M = 2 * N;
    int64_t l_M = (M + THREADS - MYTHREAD - 1) / THREADS;

    int64_t* target = (int64_t*)lgp_all_alloc(M, sizeof(int64_t));
    if (!target) return nullptr;
    ltarget = lgp_local_part(int64_t, target);
  
    for(i=0; i<l_M; i++)
        ltarget[i] = -1L;
  
    lgp_barrier();

    i = 0;
    while (i < l_N) {                // throw the darts until you get l_N hits
        r = lrand48() % M;
        if (lgp_cmp_and_swap(target, r, -1L, (i*THREADS + MYTHREAD)) == (-1L)) {
            i++;
        }
    }
  
    lgp_barrier();

    numdarts = 0;
    for (i = 0; i < l_M; i++)    // count how many darts I got
        numdarts += (ltarget[i] != -1L);

    pos = lgp_prior_add_l(numdarts);    // my first index in the perm array is the number 
                                        // of elements produce by the smaller threads
    for (i = 0; i < l_M; i++) {
        if (ltarget[i] != -1L) {
            lgp_put_int64(perm, pos, ltarget[i]);
            pos++;
        }
    }

    lgp_all_free(target);
    lgp_barrier();
    //T0_printf("done!\n");
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
    out = copied_rand_permp_agi(numrows, seed);
    t1 = wall_seconds() - t1;
    
    T0_fprintf(stderr, "rand_permp_AGI:           \n");
    lgp_min_avg_max_d(stat, t1, THREADS);
    T0_fprintf(stderr, " %8.3lf seconds\n", stat->avg);
    if (!is_perm(out, numrows)) {
        error++;
        T0_printf("\nERROR: rand_permp_%ld failed!\n\n", use_model & models_mask);
    }
    lgp_all_free(out);
  
    if (error) {
        T0_fprintf(stderr,"BALE FAIL!!!!\n"); 
    }
  
    lgp_finalize();
  
    return error;
}
