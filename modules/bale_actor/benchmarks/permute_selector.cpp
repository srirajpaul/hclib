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

/*! \file permute_matrix.upc
 * \brief Demo program that runs the variants of permute_matrix kernel. This program generates
 * a random square matrix according to the Erdos-Renyi model, two random permutations and then
 * permutes the rows and columns of the matrix according to the two random permutations.
 */

typedef struct pkg_rowcol_t{
  int64_t row;
  int64_t col;
}pkg_rowcol_t;
typedef struct pkg_rowcnt_t{
  int64_t row;
  int64_t cnt;
}pkg_rowcnt_t;
typedef struct pkg_inonz_t{
  int64_t i;
  int64_t nonz;
}pkg_inonz_t;


class PermutationSelectorOne: public hclib::Selector<1, pkg_rowcnt_t> {
    int64_t * tmprowcnts;
  public:
    int64_t lnnz;

    PermutationSelectorOne(int64_t * tmprowcnts)
      :tmprowcnts(tmprowcnts), lnnz(0) {
      mb[0].process = [this](pkg_rowcnt_t pkgrc_p, int sender_rank) {
          this->tmprowcnts[pkgrc_p.row] = pkgrc_p.cnt;
          this->lnnz += pkgrc_p.cnt;
      };
    }
};

class PermutationSelectorTwo: public hclib::Selector<1, pkg_rowcol_t> {
    int64_t * wrkoff;
    sparsemat_t * Ap;
  public:

    PermutationSelectorTwo(int64_t * wrkoff, sparsemat_t * Ap)
      :wrkoff(wrkoff), Ap(Ap) {
      mb[0].process = [this](pkg_rowcol_t pkgnz_p, int sender_rank) {
          this->Ap->lnonzero[ this->Ap->loffset[pkgnz_p.row] + this->wrkoff[pkgnz_p.row] ] = pkgnz_p.col;
          this->wrkoff[pkgnz_p.row]++;
      };
    }
};

class PermutationSelectorThree: public hclib::Selector<2, pkg_inonz_t> {
    int64_t * lcperminv;
    sparsemat_t * Ap;
  public:

    PermutationSelectorThree(int64_t * lcperminv, sparsemat_t * Ap)
      :lcperminv(lcperminv), Ap(Ap) {

      mb[0].process = [this](pkg_inonz_t pkg_p, int sender_rank) {
          pkg_inonz_t pkg_r;
          pkg_r.i = pkg_p.i;
          pkg_r.nonz = this->lcperminv[pkg_p.nonz];
          send(1, pkg_r, sender_rank);
      };

      mb[1].process = [this](pkg_inonz_t pkg_p, int sender_rank) {
          this->Ap->lnonzero[pkg_p.i] = pkg_p.nonz;
      };
    }
};

/*! \brief apply row and column permutations to a sparse matrix using conveyors
 * \param A pointer to the original matrix
 * \param rperminv pointer to the global array holding the inverse of the row permutation
 * \param cperminv pointer to the global array holding the inverse of the column permutation
 * rperminv[i] = j means that row i of A goes to row j in matrix Ap
 * cperminv[i] = j means that col i of A goes to col j in matrix Ap
 * \return a pointer to the matrix that has been produced or NULL if the model can't be used
 */
sparsemat_t * permute_matrix_selector(sparsemat_t * A, int64_t * rperminv, int64_t * cperminv) {
  sparsemat_t * Ap;

  int64_t * lrperminv = lgp_local_part(int64_t, rperminv);
  int64_t * lcperminv = lgp_local_part(int64_t, cperminv);

  //T0_printf("Permuting matrix with conveyors\n");

  /****************************************************************/
  // distribute row counts to the permuted matrix and count the number of nonzeros per thread
  // in the permuted matrix. tmprowcnts holds the post-rperminv rowcounts
  /****************************************************************/
  int64_t * tmprowcnts = (int64_t*)calloc(A->lnumrows + 1, sizeof(int64_t));

  PermutationSelectorOne *ps1_ptr = new PermutationSelectorOne(tmprowcnts);
  hclib::finish([=]() {
    ps1_ptr->start();
    for(int row=0; row<A->lnumrows; row++) {
        int64_t pe = lrperminv[row] % THREADS;
        pkg_rowcnt_t pkg_rc;
        pkg_rc.row = lrperminv[row] / THREADS;
        pkg_rc.cnt = A->loffset[row+1] - A->loffset[row];
        ps1_ptr->send(0, pkg_rc, pe);
    }
    ps1_ptr->done(0);
  });

  lgp_barrier();
  assert(A->nnz = lgp_reduce_add_l(ps1_ptr->lnnz));

  // allocate pmat to the max of the new number of nonzeros per thread
  Ap = init_matrix(A->numrows, A->numcols, ps1_ptr->lnnz);
  delete ps1_ptr;
  if(Ap == NULL) return(NULL);
  lgp_barrier();

  // convert row counts to offsets
  Ap->loffset[0] = 0;
  for(int i = 1; i < Ap->lnumrows+1; i++)
    Ap->loffset[i] = Ap->loffset[i-1] + tmprowcnts[i-1];

  /****************************************************************/
  // re-distribute nonzeros
  // working offset: wrkoff[row] gives the first empty spot on row row
  /****************************************************************/
  int64_t * wrkoff = (int64_t*)calloc(A->lnumrows, sizeof(int64_t));

  PermutationSelectorTwo *ps2_ptr = new PermutationSelectorTwo(wrkoff, Ap);
  hclib::finish([=]() {
    ps2_ptr->start();
    for(int row=0, i=0;i < A->lnnz; i++){
        while( i == A->loffset[row+1] ) // skip empty rows
            row++;
        pkg_rowcol_t pkg_nz;
        pkg_nz.row = lrperminv[row] / THREADS;
        pkg_nz.col = A->lnonzero[i];
        int64_t pe = lrperminv[row] % THREADS;
        ps2_ptr->send(0, pkg_nz, pe);
    }
    ps2_ptr->done(0);
  });

  lgp_barrier();
  delete ps2_ptr;

  /* sanity check */
  int64_t error = 0L;
  for(int i = 0; i < Ap->lnumrows; i++)
    if(wrkoff[i] != tmprowcnts[i]){printf("w[%ld] = %ld trc[%ld] = %ld\n", i, wrkoff[i], i, tmprowcnts[i]);error++;}
  if(error){printf("ERROR! permute_matrix_conveyor: error = %ld\n", error);}

  free(wrkoff);
  free(tmprowcnts);

  /****************************************************************/
  /* do column permutation ... this is essentially an indexgather */
  /****************************************************************/
  PermutationSelectorThree *ps3_ptr = new PermutationSelectorThree(lcperminv, Ap);
  hclib::finish([=]() {
    ps3_ptr->start();
    for(int i=0;i < Ap->lnnz; i++){
        pkg_inonz_t pkg_r;
        pkg_r.i = i;
        pkg_r.nonz = Ap->lnonzero[i] / THREADS;
        int64_t pe = Ap->lnonzero[i] % THREADS;
        ps3_ptr->send(0, pkg_r, pe);
    }
    ps3_ptr->done(0);
  });

  lgp_barrier();
  delete ps3_ptr;

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

  const char *deps[] = { "system", "bale_actor" };
  hclib::launch(deps, 2, [=] {

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
      assert(false);
    }

    int64_t use_model;
    sparsemat_t * outmat;
      t1 = wall_seconds();
      outmat = permute_matrix_selector(inmat, rp, cp);
      T0_fprintf(stderr,"permute_matrix_selector:           \n");
     
      t1 = wall_seconds() - t1;
      lgp_min_avg_max_d( stat, t1, THREADS );
      T0_fprintf(stderr," %8.3lf seconds\n", stat->avg);    
      clear_matrix(outmat);
      
    clear_matrix(inmat);
    lgp_all_free(rp);
    lgp_all_free(cp);
    lgp_barrier();
  });
  return 0;
}


