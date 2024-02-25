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
/*! \file triangle.upc
 * \brief Demo application that counts triangles in a graph.
 */

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

/*!
  \page triangles_page Triangles

This uses matrix algebra approach to counting triangles in a graph.

The adjacency matrix, <b>A</b>, for the graph is a {0,1} matrix
where the rows and cols correspond to the vertices
and \f$a_{ij} \f$ = <b>A[i][j]</b> is 1 exactly when there is a edge between 
vertices <i>v_i</i> and <i>v_j</i>.

The triangle with vertices <i>{v_i, v_j, v_k}</i> has associated 
edges <i>{v_i, v_j}</i>, <i>{v_j, v_k}</i> and <i>{v_k, v_i}</i> 
which correspond to non-zero entries 
\f$a_{ij}\f$,
\f$a_{jk}\f$, and
\f$a_{ki}\f$
in the adjacency matrix.  
Hence the sum
\f$ \sum_{i,j,k} a_{ij}a_{jk}a_{ki} \f$ counts the triangles in the graph.
However, it counts each triangle 6 times according to the 6 symmetries of a triangle
and the 6 symmetric ways to choose the three nonzeros in <b>A</b>.
To count each triangle once, we compute the sum
\f[ \sum_{i=1}^{n}\sum_{j=1}^{i-1}\sum_{k=1}^{j-1} a_{ij}a_{jk}a_{ik} = 
    \sum_{i=1}^{n}\sum_{j=1}^{i-1} a_{ij} \sum_{k=1}^{j-1} a_{jk}a_{ik}. \f]

This picks out a unique labelling from the 6 possible and it means that 
all the information we need about edges is contained in the lower triangular 
part of symmetric adjacency matrix.  We call this matrix <b>L</b>.

The mathematical expression: 
for each nonzero \f$ a_{ij} \f$ compute the dot product 
of row \f$ i\f$ and row \f$ j \f$ becomes
\verbatim
  For each non-zero L[i][j] 
     compute the size of the intersection of the nonzeros in row i and row j
\endverbatim

Usage:
- -a = 0,1: 0 to compute (L & L * U), 1 to compute (L & U * L).
- -e = p: specify the Erdos-Renyi probability p
- -K = str: Generate a Kronecker product graph with specified parameters. See below
- -M mask is the or of 1,2,4,8,16 for the models: agi,exstack,exstack2,conveyor,alternate
- -N = n: Specify the number of rows_per_thread in the matrix (if using the Erdos-Renyi generator).
- -r "file" : Specify a filename containing a matrix in MatrixMarket format to read as input.
- -b = count: Specify the number of packages in an exstack(2) stack

Explanation of -K option. Using a special string format you must specify a mode,
and a sequence of numbers. For example:
"0 3 4 5 9"
The first integer is always the mode and the valid modes are 0, 1, or 2
Mode 0 graphs have no triangles, mode 1 graphs have lots of triangles and mode 2 graphs
have few triangles.
After the first number, the next numbers are the parameters for the two kronecker product graphs. We split
the sequence of numbers in half to get two sequences.
In our example above we would produce the product of K(3,4) and K(5,9).

See "Design, Generation, and Validation of Extreme Scale Power-Law Graphs" by Kepner et. al.
 */

static void usage(void) {
  T0_fprintf(stderr,"\
Usage:\n\
triangle [-h][-a 0,1][-e prob][-K str][-f filename]\n\
- -a = 0,1: 0 to compute (L & L * U), 1 to compute (L & U * L)\n\
- -e = p: specify the Erdos-Renyi probability p\n\
- -h print this message\n\
- -K = str: Generate a Kronecker product graph with specified parameters. See below\n\
- -M mask is the or of 1,2,4,8,16 for the models: agi,exstack,exstack2,conveyor,alternate\n\
- -N = n: Specify the number of rows_per_thread in the matrix (if using the Erdos-Renyi generator)\n\
- -f filename : Specify a filename containing a matrix in MatrixMarket format to read as input\n\
- -b = count: Specify the number of packages in an exstack(2) stack\n\
\n\
Explanation of -K option. Using a special string format you must specify a mode,\n\
and a sequence of numbers. For example:\n\
\"0 3 4 5 9\"\n\
The first integer is always the mode and the valid modes are 0, 1, or 2\n\
Mode 0 graphs have no triangles, mode 1 graphs have lots of triangles and mode 2 graphs\n\
have few triangles.\n\
After the first number, the next numbers are the parameters for the two kronecker product graphs. We split\n\
the sequence of numbers in half to get two sequences.\n\
In our example above we would produce the product of K(3,4) and K(5,9).\n\
\n");
  upc_global_exit(0);
}


/*! \brief Generate a distributed graph that is the product of a
 collection of star graphs. This is done * in two stages. In the first
 stage, the list of stars (parameterized by an integer m, K_{1,m}) is
 split in half and each half-list forms a local adjacency matrix for
 the Kronecker product of the stars (matrices B and C). Then the two
 local matrices are combined to form a distributed adjacency matrix
 for the Kronecker product of B and C.
 *
 * \param B_spec An array of sizes in one half of the list.
 * \param B_num The number of stars in one half of the list.
 * \param C_spec An array of sizes in the other half of the list.
 * \param C_num The number of stars in the other half of the list.
 * \param mode Mode 0 graphs have no triangles, mode 1 graphs have lots of triangles and mode 2 graphs
have few triangles.
* \return A distributed matrix which represents the adjacency matrix for the Kronecker product of all the stars (B and C lists).
 */
#if 0
sparsemat_t * generate_kronecker_graph(int64_t * B_spec, int64_t B_num, int64_t * C_spec, int64_t C_num, int mode){

  T0_fprintf(stderr,"Generating Mode %d Kronecker Product graph (A = B X C) with parameters:  ", mode);
  for(int i = 0; i < B_num; i++) T0_fprintf(stderr,"%ld ", B_spec[i]);
  T0_fprintf(stderr,"X ");
  for(int i = 0; i < C_num; i++) T0_fprintf(stderr,"%ld ", C_spec[i]);   
  T0_fprintf(stderr,"\n");

  sparsemat_t * B = gen_local_mat_from_stars(B_num, B_spec, mode);
  sparsemat_t * C = gen_local_mat_from_stars(C_num, C_spec, mode);   
  if(!B || !C){
    T0_fprintf(stderr,"ERROR: triangles: error generating input!\n"); lgp_global_exit(1);
  }
  
  T0_fprintf(stderr,"B has %ld rows/cols and %ld nnz\n", B->numrows, B->lnnz);
  T0_fprintf(stderr,"C has %ld rows/cols and %ld nnz\n", C->numrows, C->lnnz);
  
  sparsemat_t * A = kron_prod_dist(B, C, 1);
  
  return(A);
}
#endif 

double triangle_upc(int64_t *count, int64_t *sr, sparsemat_t * L, sparsemat_t * U, int64_t alg) {
  int64_t cnt=0;
  int64_t numpulled=0;
  int64_t l_i, ii, k, kk, w, L_i, L_j;
   
  double t1 = wall_seconds();

  if(!L){ T0_printf("ERROR: triangle_agi: NULL L!\n"); return(-1); }
  
  if(alg == 0){    
    //foreach nonzero (i, j) in L
    for(l_i = 0; l_i < L->lnumrows; l_i++){ 
      for(k = L->loffset[l_i] + 1; k < L->loffset[l_i + 1]; k++) {
        L_i = l_i*THREADS + MYTHREAD;
        L_j = L->lnonzero[k];
        assert( L_j < L_i );

        // NOW: get L[L_j,:] and count intersections with L[L_i,:]
        int64_t start = L->loffset[l_i];
        int64_t kbegin = L->offset[L_j]; 
        int64_t kend   = L->offset[L_j + THREADS];
        numpulled+=2;
        for( kk = kbegin; kk < kend; kk++){
          w = L->nonzero[ L_j%THREADS + kk*THREADS];
          numpulled++;
          assert( w < L_j );
          for(ii = start; ii < L->loffset[l_i + 1]; ii++) {
            if( w ==  L->lnonzero[ii] ){ 
              cnt++;
              start = ii + 1;
              break;
            }
            if( w < L->lnonzero[ii] ){ // the rest are all bigger because L is tidy
              start = ii;
              break;
            }
          }
        }
      }
    }
  }else{
    if(!U){ T0_printf("ERROR: triangle_agi: NULL U!\n"); return(-1); }
    
    //foreach nonzero (i, j) in L
    for(l_i = 0; l_i < L->lnumrows; l_i++){ 
      for(k = L->loffset[l_i]; k < L->loffset[l_i + 1]; k++) {
        L_i = l_i*THREADS + MYTHREAD;
        L_j = L->lnonzero[k];
        assert( L_j < L_i );
        
        // NOW: get U[L_j,:] and count intersections with U[L_i,:]
        int64_t start = U->loffset[l_i + 1] - 1;
        int64_t smallest_col_in_row_i = U->lnonzero[U->loffset[l_i]];
        numpulled+=2;
        int64_t kbegin = U->offset[ L_j + THREADS] - 1;
        int64_t kend   = U->offset[L_j];
        for( kk = kbegin; kk >= kend; kk--){
          w = U->nonzero[ L_j%THREADS + kk*THREADS];
          numpulled++;
          if (w < smallest_col_in_row_i) break; // there can't be any more intersections in these rows
          assert( w >= L_j );
          for(ii = start; ii >= U->loffset[l_i]; ii--) {
            if( w ==  U->lnonzero[ii] ){ 
              cnt++;
              start = ii-1;
              break;
            }
            if( w > U->lnonzero[ii] ){ // the rest are all smaller because U is tidy
              start = ii;
              break;
            }
          }
        }
      }
    }
    
  }

  upc_barrier;
  t1 = wall_seconds() - t1;
  double stat = myupc_reduce_add_d( t1)/ THREADS;

  *sr = numpulled; 
  *count = cnt;
  return stat;
}

int main(int argc, char * argv[]) {

  int64_t models_mask = 0xF;  // default is running all models
  int64_t l_numrows = 10000;         // number of a rows per thread
  int64_t nz_per_row = 35;           // target number of nonzeros per row (only for Erdos-Renyi)
  int64_t read_graph = 0L;           // read graph from a file
  char filename[64];
  
  double t1;
  int64_t i, j;
  int64_t alg = 0;
  int64_t gen_kron_graph = 0L;
  int kron_graph_mode = 0;
  char * kron_graph_string;
  double erdos_renyi_prob = 0.0;
  
  int printhelp = 0;
  int opt; 
  while( (opt = getopt(argc, argv, "hM:n:f:a:e:K:")) != -1 ) {
    switch(opt) {
    case 'h': printhelp = 1; break;
    case 'M': sscanf(optarg,"%ld", &models_mask);  break;
    case 'n': sscanf(optarg,"%ld", &l_numrows); break;
    case 'f': read_graph = 1; sscanf(optarg,"%s", filename); break;

    case 'a': sscanf(optarg,"%ld", &alg); break;
    case 'e': sscanf(optarg,"%lg", &erdos_renyi_prob); break;
    case 'K': gen_kron_graph = 1; kron_graph_string = optarg; break;
    default:  break;
    }
  }
  if( printhelp ) usage(); 

  upc_atomic_domain = upc_all_atomicdomain_alloc(UPC_INT64, UPC_ADD|UPC_INC|UPC_CSWAP, 0);
  int64_t numrows = l_numrows * THREADS;
  if(erdos_renyi_prob == 0.0){ // use nz_per_row to get erdos_renyi_prob
    erdos_renyi_prob = (2.0*(nz_per_row - 1))/numrows;
    if(erdos_renyi_prob > 1.0)
      erdos_renyi_prob = 1.0;
  } else {                     // use erdos_renyi_prob to get nz_per_row
    nz_per_row = erdos_renyi_prob * numrows;
  }
  
  T0_fprintf(stderr,"Running triangle on %d PEs\n", THREADS);
  if(!read_graph && !gen_kron_graph){
    T0_fprintf(stderr,"Number of rows per PE   (-N)   %ld\n", l_numrows);
    T0_fprintf(stderr,"Erdos Renyi prob (-e)   %g\n", erdos_renyi_prob);
  }
  T0_fprintf(stderr,"Model mask (M) = %ld (should be 1,2,4,8,16 for agi, exstack, exstack2, conveyors, alternates\n", models_mask);  
  T0_fprintf(stderr,"algorithm (a) = %ld (0 for L & L*U, 1 for L & U*L)\n", alg);
  
  if( printhelp )
    upc_global_exit(0);

  // alg = 0 only needs L
  // alg = 1 needs both U and L

  double correct_answer = -1;
  
  sparsemat_t *A, *L, *U;
  if(read_graph){
    T0_fprintf(stderr, "Reading graph from file not supported\n");
    assert(0);
#if 0
    A = read_matrix_mm_to_dist(filename);
    if(!A)
      lgp_global_exit(1);
    T0_fprintf(stderr,"Reading file %s...\n", filename);
    T0_fprintf(stderr, "A has %ld rows/cols and %ld nonzeros.\n", A->numrows, A->nnz);

    // we should check that A is symmetric!
    
    if(!is_lower_triangular(A, 0)){ //if A is not lower triangular... make it so.      
      T0_fprintf(stderr, "Assuming symmetric matrix... using lower-triangular portion...\n");
      tril(A, -1);
      L = A;
    }else      
      L = A;
    
    sort_nonzeros(L);
#endif
  }else if (gen_kron_graph){
    T0_fprintf(stderr, "kronecker graph not supported\n");
    assert(0);
#if 0
    // string should be <mode> # # ... #
    // we will break the string of numbers (#s) into two groups and create
    // two local kronecker graphs out of them.
    int num;
    char * ptr = kron_graph_string;
    int64_t * kron_specs = (int64_t *)calloc(32, sizeof(int64_t *));
    
    // read the mode
    int ret = sscanf(ptr, "%d ", &kron_graph_mode);
    if(ret == 0)ret = sscanf(ptr, "\"%d ", &kron_graph_mode);
    if(ret == 0){T0_fprintf(stderr,"ERROR reading kron graph string!\n");return(1);}
    //T0_fprintf(stderr,"kron string: %s return = %d\n", ptr, ret);
    //T0_fprintf(stderr,"kron mode: %d\n", kron_graph_mode);
    ptr+=2;
    int mat, num_ints = 0;
    while(sscanf(ptr, "%d %n", &num, &mat) == 1){
      //T0_fprintf(stderr,"%s %d\n", ptr, mat);
      kron_specs[num_ints++] = num;
      ptr+=mat;
    }

    if(num_ints <= 1){
      T0_fprintf(stderr,"ERROR: invalid kronecker product string (%s): must contain at least three integers\n", kron_graph_string); return(-1);}


    /* calculate the number of triangles */
    if(kron_graph_mode == 0){
      correct_answer = 0.0;
    }else if(kron_graph_mode == 1){
      correct_answer = 1;
      for(i = 0; i < num_ints; i++)
        correct_answer *= (3*kron_specs[i] + 1);
      correct_answer *= 1.0/6.0;
      double x = 1;
      for(i = 0; i < num_ints; i++)
        x *= (kron_specs[i] + 1);
      correct_answer = correct_answer - 0.5*x + 1.0/3.0;
    }else if(kron_graph_mode == 2){
      correct_answer = (1.0/6.0)*pow(4,num_ints) - pow(2.0,(num_ints - 1)) + 1.0/3.0;
    }

    correct_answer = round(correct_answer);
    T0_fprintf(stderr, "Pre-calculated answer = %ld\n", (int64_t)correct_answer);
    
    int64_t half = num_ints/2;
    
    L = generate_kronecker_graph(kron_specs, half, &kron_specs[half], num_ints - half, kron_graph_mode);
#endif
  }else{
    L = gen_erdos_renyi_graph_dist(numrows, erdos_renyi_prob, 0, 1, 12345);
  }

  upc_barrier;
  
  if(alg == 1)
    U = transpose_matrix(L);

  upc_barrier;

  T0_fprintf(stderr,"L has %ld rows/cols and %ld nonzeros.\n", L->numrows, L->nnz);
  
  if(!is_lower_triangular(L, 0)){
    T0_fprintf(stderr,"ERROR: L is not lower triangular!\n");
    upc_global_exit(1);
  }
    
  T0_fprintf(stderr,"Run triangle counting ...\n");
  int64_t tri_cnt;           // partial count of triangles on this thread
  int64_t total_tri_cnt;     // the total number of triangles on all threads
  int64_t sh_refs;         // number of shared reference or pushes
  int64_t total_sh_refs;

  
  shared int64_t * cc = upc_all_alloc( L->numrows, sizeof(int64_t));
  int64_t * l_cc = (int64_t*) (cc+MYTHREAD);
  for(i = 0; i < L->lnumrows; i++)
    l_cc[i] = 0;
  upc_barrier;
  
  /* calculate col sums */
  for(i = 0; i < L->lnnz; i++){
    int64_t ret;
    upc_atomic_relaxed(upc_atomic_domain, &ret, UPC_INC, &cc[L->lnonzero[i]], NULL, NULL);
  }
  
  upc_barrier;
  
  int64_t rtimesc_calc = 0;
  for(i = 0; i < L->lnumrows; i++){
    int64_t deg = L->loffset[i + 1] - L->loffset[i];        
    rtimesc_calc += deg*l_cc[i];
  }

  /* calculate sum (r_i choose 2) */
  int64_t rchoose2_calc = 0;
  for(i = 0; i < L->lnumrows; i++){
    int64_t deg = L->loffset[i + 1] - L->loffset[i];
    rchoose2_calc += deg*(deg-1)/2;
  }
  
  /* calculate sum (c_i choose 2) */
  int64_t cchoose2_calc = 0;
  for(i = 0; i < L->lnumrows; i++){
    int64_t deg = l_cc[i];
    cchoose2_calc += deg*(deg-1)/2;
  }
  int64_t pulls_calc = 0;
  int64_t pushes_calc = 0;
  if(alg == 0){
    pulls_calc = myupc_reduce_add_l(rtimesc_calc);
    pushes_calc = myupc_reduce_add_l(rchoose2_calc);
  }else{
    pushes_calc = myupc_reduce_add_l(rtimesc_calc);
    pulls_calc = myupc_reduce_add_l(cchoose2_calc);
  }

  upc_all_free(cc);
  
  T0_fprintf(stderr,"Calculated: Pulls = %ld\n            Pushes = %ld\n\n",pulls_calc, pushes_calc);
  
  int64_t use_model;
  double laptime = 0.0;

  tri_cnt = 0;
  total_tri_cnt = 0;
  sh_refs = 0;
  total_sh_refs = 0;
  T0_fprintf(stderr,"      upc: \n");
  laptime = triangle_upc(&tri_cnt, &sh_refs, L, U, alg);
    
  upc_barrier;
  total_tri_cnt = myupc_reduce_add_l(tri_cnt);
  total_sh_refs = myupc_reduce_add_l(sh_refs);
  T0_fprintf(stderr,"  %8.3lf seconds: %16ld triangles\n", laptime, total_tri_cnt);
  //T0_fprintf(stderr,"%16ld shared refs\n", total_sh_refs);
  if((correct_answer >= 0) && (total_tri_cnt != (int64_t)correct_answer)){
    T0_fprintf(stderr, "ERROR: Wrong answer!\n");
  }
  
  if(correct_answer == -1) {
    correct_answer = total_tri_cnt;
  }
  
  upc_barrier;
  return(0);
}

