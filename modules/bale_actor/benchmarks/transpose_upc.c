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
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
//#include <upc.h>
#include <upc_relaxed.h>
#include <upc_atomic.h>
#include <upc_collective.h>
#include "myupc.h"

/*! \file transpose_matrix.upc
 * \brief Demo program that runs the variants of transpose_matrix kernel.
 */

/*! 
 * \page transpose_matrix_page Transpose Matrix
 *
 * Demo program that runs the variants of transpose matrix kernel. This program
 * generates a random square asymmetrix (via the Erdos-Renyi model) and then transposes
 * it in parallel.
 *
 * See files spmat_agi.upc, spmat_exstack.upc, spmat_exstack2.upc, and spmat_conveyor.upc
 * for the source for the kernels.
 * 
 * Usage:
 * transpose_matrix [-h][-b count][-M mask][-n num][-T tabsize][-c num]
 * - -h Print usage banner only
 * - -b count is the number of packages in an exstack(2) buffer
 * - -e=p Set the Erdos-Renyi probability to p.
 * - -M=m Set the models mask (1,2,4,8,16,32 for AGI, exstack, exstack2,conveyor,alternate)
 * - -n=n Set the number of rows per PE to n (default = 1000).
 * - -s=s Set a seed for the random number generation.
 * - -Z=z Set the avg number of nonzeros per row to z (default = 10, overrides Erdos-Renyi p).
 */
sparsemat_t * transpose_matrix_upc(sparsemat_t *A) {
  int64_t counted_nnz_At;
  int64_t lnnz, i, j, col, row, fromth, idx;
  int64_t pos;
  sparsemat_t * At;
  
  //T0_printf("UPC version of matrix transpose...");
  
  // find the number of nnz.s per thread

  shared int64_t * shtmp = upc_all_alloc( A->numcols + THREADS, sizeof(int64_t));
  if( shtmp == NULL ) return(NULL);
  int64_t * l_shtmp = (int64_t*)(shtmp+MYTHREAD);
  int64_t lnc = (A->numcols + THREADS - MYTHREAD - 1)/THREADS;
  for(i=0; i < lnc; i++)
    l_shtmp[i] = 0;
  upc_barrier;

  for( i=0; i< A->lnnz; i++) {                   // histogram the column counts of A
    assert( A->lnonzero[i] < A->numcols );
    assert( A->lnonzero[i] >= 0 ); 
    upc_atomic_relaxed(upc_atomic_domain, &pos, UPC_INC, &shtmp[A->lnonzero[i]], NULL, NULL);
  }
  upc_barrier;


  lnnz = 0;
  for( i = 0; i < lnc; i++) {
    lnnz += l_shtmp[i];
  }
  
  At = init_matrix(A->numcols, A->numrows, lnnz);
  if(!At){printf("ERROR: transpose_matrix_upc: init_matrix failed!\n");return(NULL);}

  int64_t sum = myupc_reduce_add_l(lnnz);      // check the histogram counted everything
  assert( A->nnz == sum ); 

  // compute the local offsets
  At->loffset[0] = 0;
  for(i = 1; i < At->lnumrows+1; i++)
    At->loffset[i] = At->loffset[i-1] + l_shtmp[i-1];

  // get the global indices of the start of each row of At
  for(i = 0; i < At->lnumrows; i++)
    l_shtmp[i] = MYTHREAD + THREADS * (At->loffset[i]);
    
  upc_barrier;

  //redistribute the nonzeros 
  for(row=0; row<A->lnumrows; row++) {
    for(j=A->loffset[row]; j<A->loffset[row+1]; j++){
      int64_t add_value = (int64_t) THREADS;
      upc_atomic_relaxed(upc_atomic_domain, &pos, UPC_ADD, &shtmp[A->lnonzero[j]], &add_value, NULL);
      (At->nonzero)[pos] = row*THREADS + MYTHREAD;
    }
  }
  upc_barrier;

  //upc_barrier;
  //if(!MYTHREAD)printf("done\n");
  upc_all_free(shtmp);

  return(At);
}

int main(int argc, char** argv) {

    int64_t i;
    int64_t check = 1;
    int64_t models_mask = 0xF;
    int printhelp = 0;
    double erdos_renyi_prob = 0.0;
    int64_t nz_per_row = -1;
    int64_t l_numrows = 10000;
    int64_t numrows;
    int64_t seed = 101892+MYTHREAD;
    sparsemat_t * inmat, * outmat;

    int opt; 
    while( (opt = getopt(argc, argv, "hCe:n:M:s:Z:")) != -1 ) {
        switch(opt) {
            case 'h': printhelp = 1; break;
            case 'C': check = 0; break;
            case 'e': sscanf(optarg,"%lf", &erdos_renyi_prob); break;
            case 'n': sscanf(optarg,"%ld", &l_numrows);   break;
            case 'M': sscanf(optarg,"%ld", &models_mask);  break;
            case 's': sscanf(optarg,"%ld", &seed); break;
            case 'Z': sscanf(optarg,"%ld", &nz_per_row);  break;
            default:  break;
        }
    }
  
    upc_atomic_domain = upc_all_atomicdomain_alloc(UPC_INT64, UPC_ADD|UPC_INC|UPC_CSWAP, 0);
    numrows = l_numrows * THREADS;

    /* set erdos_renyi_prob and nz_per_row to be consistent */
    if (nz_per_row == -1 && erdos_renyi_prob == 0) {
        nz_per_row = 10;
    } else if (nz_per_row == -1) {
        nz_per_row = erdos_renyi_prob*numrows;
    }
    
    erdos_renyi_prob = (2.0 * (nz_per_row - 1)) / numrows;
    if (erdos_renyi_prob > 1.0)
        erdos_renyi_prob = 1.0;

    T0_fprintf(stderr,"Running transpose_matrix on %d PEs\n", THREADS);
    T0_fprintf(stderr,"Erdos-Renyi edge probability(-e) = %lf\n", erdos_renyi_prob);
    T0_fprintf(stderr,"rows per PE (-n)                 = %ld\n", l_numrows);
    T0_fprintf(stderr,"models_mask (-M)                 = %ld or of 1,2,4,8,16 for gets,classic,exstack2,conveyor,alternate\n", models_mask);
    T0_fprintf(stderr,"seed (-s)                        = %ld\n", seed);
    T0_fprintf(stderr,"Avg # of nonzeros per row   (-Z) = %ld\n", nz_per_row);


    double t1;
    int64_t error = 0;
    
    inmat = gen_erdos_renyi_graph_dist(numrows, erdos_renyi_prob, 0, 3, seed + 2);
    if (!inmat) {
        T0_printf("ERROR: inmat is null!\n");

        return -1;
    }
  
    t1 = wall_seconds();
    outmat = transpose_matrix_upc(inmat);
    T0_fprintf(stderr, "upc:     \n");
    t1 = wall_seconds() - t1;

    double stat = myupc_reduce_add_d(t1)/ THREADS;
    T0_fprintf(stderr, " %8.3lf  seconds\n", stat);

    /* correctness check */
    if (check) {      
        sparsemat_t * outmatT = transpose_matrix(outmat);
        if (compare_matrix(outmatT, inmat)) {
            T0_fprintf(stderr,"ERROR: transpose of transpose does not match!\n");
        }
        clear_matrix(outmatT);
    }

    clear_matrix(outmat);
    clear_matrix(inmat);
    upc_barrier;
    return(error);
}
