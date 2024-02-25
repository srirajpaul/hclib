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

/*! \file permute_matrix.upc
 * \brief Demo program that runs the variants of permute_matrix kernel. This program generates
 * a random square matrix according to the Erdos-Renyi model, two random permutations and then
 * permutes the rows and columns of the matrix according to the two random permutations.
 */

/*! \brief apply row and column permutations to a sparse matrix using straight UPC
 * \param A pointer to the original matrix
 * \param rperminv pointer to the global array holding the inverse of the row permutation
 * \param cperminv pointer to the global array holding the inverse of the column permutation
 * rperminv[i] = j means that row i of A goes to row j in matrix Ap
 * cperminv[i] = j means that col i of A goes to col j in matrix Ap
 * \return a pointer to the matrix that has been produced or NULL if the model can't be used
 * \ingroup spmatgrp
 */
sparsemat_t * copied_permute_matrix_agi(sparsemat_t *A, int64_t *rperminv, int64_t *cperminv) {
  //T0_printf("Permuting matrix with single puts\n");
  int64_t i, j, col, row, pos;
  int64_t * lrperminv = lgp_local_part(int64_t, rperminv);
  int64_t * rperm = (int64_t*)lgp_all_alloc(A->numrows, sizeof(int64_t));
  if( rperm == NULL ) return(NULL);
  int64_t *lrperm = lgp_local_part(int64_t, rperm);

  //compute rperm from rperminv 
  for(i=0; i < A->lnumrows; i++){
    lgp_put_int64(rperm, lrperminv[i], i*THREADS + MYTHREAD);
  }

  lgp_barrier();
  
  int64_t cnt = 0, off, nxtoff;
  for(i = 0; i < A->lnumrows; i++){
    row = lrperm[i];
    off    = lgp_get_int64(A->offset, row);
    nxtoff = lgp_get_int64(A->offset, row + THREADS);
    cnt += nxtoff - off;
  }
  lgp_barrier();

  sparsemat_t * Ap = init_matrix(A->numrows, A->numcols, cnt);
  
  // fill in permuted rows
  Ap->loffset[0] = pos = 0;
  for(i = 0; i < Ap->lnumrows; i++){
    row = lrperm[i];
    off    = lgp_get_int64(A->offset, row);
    nxtoff = lgp_get_int64(A->offset, row + THREADS);
    for(j = off; j < nxtoff; j++){
      Ap->lnonzero[pos++] = lgp_get_int64(A->nonzero, j*THREADS + row%THREADS);
    }
    Ap->loffset[i+1] = pos;
  }
  
  assert(pos == cnt);

  lgp_barrier();
  
  // finally permute column indices
  for(i = 0; i < Ap->lnumrows; i++){
    for(j = Ap->loffset[i]; j < Ap->loffset[i+1]; j++){
      Ap->lnonzero[j] = lgp_get_int64(cperminv, Ap->lnonzero[j]);      
    }
  }
  lgp_barrier();

  lgp_all_free(rperm);

  T0_printf("done\n");
  return(Ap);
}
  
/*!
 * 
 * \page permute_matrix_page Permute Matrix
 *
 * This is a demo program that runs the variants of permute_matrix kernel. This program generates
 * a random square matrix according to the Erdos-Renyi model, two random permutations and then
 * permutes the rows and columns of the matrix according to the two random permutations.
 *
 * See files spmat_agi.upc, spmat_exstack.upc, spmat_exstack2.upc, and spmat_conveyor.upc
 * for the source for the kernels.
 * 
 * Usage:
 * permute_matrix [-h][-e prob][-M mask][-n num][-s seed][-Z num]
 * - -h Print usage banner only
 * - -b count is the number of packages in an exstack(2) buffer
 * - -e=p Set the Erdos-Renyi probability to p.
 * - -M=m Set the models mask (1,2,4,8,16,32 for AGI, exstack, exstack2,conveyor,alternate)
 * - -n=n Set the number of rows per PE to n (default = 1000).
 * - -s=s Set a seed for the random number generation.
 * - -Z=z Set the avg number of nonzeros per row to z (default = 10, overrides Erdos-Renyi p).
 *
 */

static void usage(void) {
  T0_fprintf(stderr,"\
This is a demo program that runs various implementations of the permute_matrix kernel.\n\
Usage:\n\
permute_matrix [-h][-e prob][-M mask][-n num][-s seed][-Z num]\n\
 -h Print usage banner only\n\
 -b count is the number of packages in an exstack(2) buffer\n\
 -e=prob Set the Erdos-Renyi probability to p.\n\
 -M=mask Set the models mask (1,2,4,8,16,32 for AGI, exstack, exstack2,conveyor,alternate)\n\
 -n=num Set the number of rows per PE to n (default = 1000).\n\
 -s=seed Set a seed for the random number generation.\n\
 -Z=num  Set the avg number of nonzeros per row to z (default = 10, overrides Erdos-Renyi p).\n\
\n");
  lgp_global_exit(0);
}


int main(int argc, char * argv[]) {
  lgp_init(argc, argv);

  int64_t i;
  int64_t models_mask = 0xF;
  int printhelp = 0;
  double erdos_renyi_prob = 0.0;
  int64_t buf_cnt = 1024;
  int64_t nz_per_row = -1;
  int64_t l_numrows = 10000;
  int64_t numrows;
  int64_t seed = 101892+MYTHREAD;
  int64_t cores_per_node = 1;

  int opt; 
  while( (opt = getopt(argc, argv, "e:hc:n:M:s:Z:")) != -1 ) {
    switch(opt) {
    case 'h': printhelp = 1; break;
    case 'b': sscanf(optarg,"%ld", &buf_cnt);  break;
    case 'c': sscanf(optarg,"%ld" ,&cores_per_node); break;
    case 'e': sscanf(optarg,"%lf", &erdos_renyi_prob);  break;
    case 'n': sscanf(optarg,"%ld", &l_numrows);   break;
    case 'M': sscanf(optarg,"%ld", &models_mask);  break;
    case 's': sscanf(optarg,"%ld", &seed); break;
    case 'Z': sscanf(optarg,"%ld", &nz_per_row);  break;
    default:  break;
    }
  }
  if(printhelp) usage();

  numrows = l_numrows * THREADS;

  /* set erdos_renyi_prob and nz_per_row to be consistent */
  if(nz_per_row == -1 && erdos_renyi_prob == 0.0){
    nz_per_row = 10;
  }else if(nz_per_row == -1){
    nz_per_row = erdos_renyi_prob*numrows;
  }
  erdos_renyi_prob = (2.0*(nz_per_row - 1))/numrows;
  if(erdos_renyi_prob > 1.0)
    erdos_renyi_prob = 1.0;  
  
  T0_fprintf(stderr,"Running permute_matrix on %d PEs\n", THREADS);
  T0_fprintf(stderr,"buf_cnt (stack size)         (-b)= %ld\n", buf_cnt);
  T0_fprintf(stderr,"Erdos-Renyi edge probability (-e)= %lf\n", erdos_renyi_prob);
  T0_fprintf(stderr,"rows per PE (-n)             = %ld\n", l_numrows);
  T0_fprintf(stderr,"Avg # of nonzeros per row    (-Z)= %ld\n", nz_per_row);
  T0_fprintf(stderr,"models_mask (-M)                 = %ld or of 1,2,4,8,16,32 for gets,classic,exstack2,conveyor,ex2cyclic,ex2goto\n", models_mask);
  T0_fprintf(stderr,"seed (-s)                        = %ld\n", seed);


  double t1;
  minavgmaxD_t stat[1];
  int64_t error = 0;
  
  int64_t * rp = rand_permp(numrows, seed);
  int64_t * cp = rand_permp(numrows, seed + 1);  
  
  sparsemat_t * inmat = gen_erdos_renyi_graph_dist(numrows, erdos_renyi_prob, 0, 3, seed + 2);
  if(inmat == NULL){
    T0_printf("ERROR: inmat is null!\n");
    return(-1);
  }

  int64_t use_model;
  sparsemat_t * outmat;
    t1 = wall_seconds();
    outmat = copied_permute_matrix_agi(inmat, rp, cp);
    T0_fprintf(stderr,"permute_matrix_AGI:           \n");
   
    t1 = wall_seconds() - t1;
    lgp_min_avg_max_d( stat, t1, THREADS );
    T0_fprintf(stderr," %8.3lf seconds\n", stat->avg);    
    clear_matrix(outmat);
    
  clear_matrix(inmat);
  lgp_all_free(rp);
  lgp_all_free(cp);
  lgp_barrier();
  lgp_finalize();
  return(error);
}


