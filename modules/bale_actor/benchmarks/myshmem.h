
#ifndef MYSHMEM_H
#define MYSHMEM_H

#define THREADS shmem_n_pes()
#define MYTHREAD shmem_my_pe()

#include <time.h>
#include <sys/time.h>

typedef enum graph_model {FLAT, GEOMETRIC, KRONECKER} graph_model;
typedef enum edge_type {DIRECTED, UNDIRECTED, DIRECTED_WEIGHTED, UNDIRECTED_WEIGHTED} edge_type;
typedef enum self_loops {NOLOOPS, LOOPS} self_loops;

double wall_seconds() {
  struct timeval tp;
  int retVal = gettimeofday(&tp,NULL);
  if (retVal == -1) { perror("gettimeofday:"); fflush(stderr); }
  return ( (double) tp.tv_sec + (double) tp.tv_usec * 1.e-6 );
}

//#define USE_KNUTH   /*!< Default define to set whether we use the Knuth random number generator or rand48 */
#ifdef USE_KNUTH
#define LGP_RAND_MAX 2251799813685248  /*!< max random number depends on which rng we use */
#include "knuth_rng_double_2019.h"
#else
#define LGP_RAND_MAX 281474976710656
#endif

/*! 
 * \brief seed for the random number generator
 * \param seed the seed
 * Note: if all thread call this with the same seed they actually get different seeds.
 */
void rand_seed(int64_t seed){
#ifdef USE_KNUTH
  ranf_start(seed + 1 + MYTHREAD);
#else
  srand48(seed + 1 + MYTHREAD);
#endif
}

/*! 
 * \brief return a random integer mod N.
 * \param N the modulus
 */
int64_t rand_int64(int64_t N){
  assert(N < LGP_RAND_MAX);
#ifdef USE_KNUTH
  return((int64_t)(ranf_arr_next()*N));
#else
  return((int64_t)(drand48()*N));
#endif
}

/*! 
 * \brief return a random double in the interval (0,1]
 */
double rand_double(){
#ifdef USE_KNUTH
  return(ranf_arr_next());
#else
  return(drand48());
#endif
}

#define T0_fprintf if(MYTHREAD==0) fprintf

#define T0_printf if(MYTHREAD==0) printf

#define myshmem_global_malloc(nblock_on_all_threads,blocksize) shmem_malloc( ( ( (size_t)(nblock_on_all_threads)+(shmem_n_pes())-1)/(shmem_n_pes()) ) *(blocksize) ) \

void myshmem_int64_p(int64_t *addr, size_t index, int64_t val) {
    int64_t pe = index % shmem_n_pes();
    int64_t local_index = index / shmem_n_pes();
    shmem_int64_p(addr+local_index, val, pe);
}

int64_t myshmem_int64_g(int64_t *addr, size_t index) {
    int64_t pe = index % shmem_n_pes();
    int64_t local_index = index / shmem_n_pes();
    return shmem_int64_g(addr+local_index, pe);
}

double myshmem_double_g(double *addr, size_t index) {
    int64_t pe = index % shmem_n_pes();
    int64_t local_index = index / shmem_n_pes();
    return shmem_double_g(addr+local_index, pe);
}

int64_t myshmem_fetch_and_inc(int64_t * ptr, int64_t index) {
    int64_t ret;
    int64_t lindex = index/shmem_n_pes();
    int64_t pe = index % shmem_n_pes();
    ret = shmem_int64_atomic_fetch_inc(&ptr[lindex], pe);
    return ret;
}

int64_t myshmem_fetch_and_add(int64_t * ptr, int64_t index, int64_t value) {
    int64_t ret;
    int64_t lindex = index/shmem_n_pes();
    int64_t pe = index % shmem_n_pes();
    ret = shmem_int64_atomic_fetch_add(&ptr[lindex], value, pe);
    return ret;
}

static void *setup_shmem_reduce_workdata(long **psync, size_t xsize) {
  void *work;
  int i;

  work=shmem_malloc(_SHMEM_REDUCE_MIN_WRKDATA_SIZE*xsize);
  *psync=(long *)shmem_malloc(_SHMEM_REDUCE_SYNC_SIZE*sizeof(long));
  for(i=0;i<_SHMEM_REDUCE_SYNC_SIZE;++i) {
    (*psync)[i]=_SHMEM_SYNC_VALUE;
  }
  shmem_barrier_all();
  return work;
}

double myshmem_reduce_add(double myval){
  static double *buff=NULL, *work;
  static long *sync;
  if (buff==NULL) {
    buff=(double*)shmem_malloc(2*sizeof(double));
    work=(double*)setup_shmem_reduce_workdata(&sync,sizeof(double));
  }
  buff[0]=myval;

  shmem_double_sum_to_all(&buff[1],buff,1,0,0,shmem_n_pes(),work,sync);

  shmem_barrier_all();
  return buff[1];
}

int64_t myshmem_reduce_add(int64_t myval){
  static int64_t *buff=NULL, *work;
  static long *sync;
  if (buff==NULL) {
    buff=(int64_t*)shmem_malloc(2*sizeof(int64_t));
    work=(int64_t*)setup_shmem_reduce_workdata(&sync,sizeof(int64_t));
  }
  buff[0]=myval;

  shmem_long_sum_to_all(&buff[1],buff,1,0,0,shmem_n_pes(),work,sync);

  shmem_barrier_all();
  return buff[1];
}

int64_t myshmem_reduce_max(int64_t myval){
  static int64_t *buff=NULL, *work;
  static long *sync;
  if (buff==NULL) {
    buff=(int64_t*)shmem_malloc(2*sizeof(int64_t));
    work=(int64_t*)setup_shmem_reduce_workdata(&sync,sizeof(int64_t));
  }
  buff[0]=myval;

  shmem_long_max_to_all(&buff[1],buff,1,0,0,shmem_n_pes(),work,sync);

  shmem_barrier_all();
  return buff[1];
}

/*! \struct sparsemat_t spmat.h
 * \brief A structure to hold a sparse matrix.
 *
 * We use a distributed version of the standard Compressed Sparse Row
 * (CSR) format.  Since the applications in bale (so far) only need to
 * know whether a matrix entry is zero or nonzero, we don't keep track
 * of the values themselves for the nonzeros in the matrix only the
 * position of the nonzeros.  This reduces the amount of memory needed
 * to store the matrices.  This also saves a lot local sorting and
 * combining that would be required to actually handle values (other
 * than one).
 
 * We store the nonzeros with affinity by rows.  That is, if row i has
 * affinity to PE t, then all the nonzeros in row i have affinity to
 * PE t. The affinity rule for a row is as follows. Row i is
 * assigned to PE (i % NPES).
 *
 * We call a matrix "tidy" if the column indices of nonzeros in every
 * row are sorted in ascending order.
 *
 *  \ingroup spmatgrp
 */
typedef struct sparsemat_t {
  int64_t local;                //!< 0/1 flag specifies whether this is a local or distributed matrix
  int64_t numrows;              //!< the total number of rows in the matrix
  int64_t lnumrows;             //!< the number of rows on this PE 
                                // note lnumrows = (numrows / NPES) + {0 or 1} 
                                //    depending on numrows % NPES
  int64_t numcols;              //!< the nonzeros have values between 0 and numcols
  int64_t nnz;                  //!< total number of nonzeros in the matrix
  int64_t lnnz;                 //!< the number of nonzeros on this PE
  int64_t * offset;      //!< the row offsets into the array of nonzeros
  int64_t * loffset;            //!< the row offsets for the row with affinity to this PE
  int64_t * nonzero;     //!< the global array of nonzeros
  int64_t * lnonzero;           //!< the nonzeros with affinity to this PE
  double * value;        //!< the global array of values of nonzeros. Optional.
  double * lvalue;              //!< local array of values (values for rows on this PE)
}sparsemat_t;

/*! \brief initializes the struct that holds a sparse matrix
 *    given the total number of rows and columns and the local number of non-zeros
 * \param numrows total number of rows
 * \param numcols total number of columns
 * \param nnz_this_thread number of nonzero on this thread
 * \return An initialized sparsemat_t or NULL on error.
 * \ingroup spmatgrp
 */
sparsemat_t * init_matrix(int64_t numrows, int64_t numcols, int64_t nnz_this_thread, int weighted) {
  sparsemat_t * mat = (sparsemat_t *)calloc(1, sizeof(sparsemat_t));
  mat->local = 0;
  mat->numrows  = numrows;
  mat->lnumrows = (numrows + THREADS - MYTHREAD - 1)/THREADS;
  mat->numcols  = numcols;  
  mat->offset   = (int64_t*)myshmem_global_malloc(mat->numrows + THREADS, sizeof(int64_t));
  if(mat->offset == NULL){
    T0_printf("ERROR: init_matrix: could not allocate %ld bytes for offset array\n", mat->numrows*8);
    return(NULL);
  }
  mat->loffset  =  mat->offset;
  int64_t max = myshmem_reduce_max(nnz_this_thread);
  int64_t total = myshmem_reduce_add(nnz_this_thread);
  mat->nonzero = (int64_t*)myshmem_global_malloc(max*THREADS, sizeof(int64_t));
  if(mat->nonzero == NULL){
    T0_printf("ERROR: init_matrix: could not allocate %ld bytes for nonzero array (max = %ld)\n", max*THREADS*8, max);
    return(NULL);
  }
  mat->lnonzero = mat->nonzero;

  if(weighted){
    mat->value = (double*)myshmem_global_malloc(max*THREADS, sizeof(double));
    if(mat->value == NULL){
      T0_printf("ERROR: init_matrix: could not allocate %ld bytes for value array (max = %ld)\n",
                max*THREADS*8, max);
      return(NULL);
    }
    mat->lvalue = mat->value;
  }else{
    mat->value = NULL;
    mat->lvalue = NULL;
  }
 
  mat->nnz = total;
  mat->lnnz = nnz_this_thread;

  return(mat);
}

/*! \brief frees the space allocated for a sparse matrix
 * \param mat pointer to the sparse matrix
 * \ingroup spmatgrp
 */
void clear_matrix(sparsemat_t * mat) {
  if(mat->local){
    free(mat->lnonzero);
    free(mat->loffset);
  }else{
    shmem_free(mat->nonzero);
    shmem_free(mat->offset);
  }
}

/*! \brief produce the transpose of a sparse matrix using UPC
 * \param A  pointer to the original matrix
 * \return a pointer to the matrix that has been produced or NULL if the model can't be used
 * \ingroup spmatgrp
 */
sparsemat_t * transpose_matrix_agi(sparsemat_t *A) {
  int64_t counted_nnz_At;
  int64_t lnnz, i, j, col, row, fromth, idx;
  int64_t pos;
  sparsemat_t * At;
  
  //T0_printf("UPC version of matrix transpose...");
  
  // find the number of nnz.s per thread

  int64_t * shtmp = (int64_t *)myshmem_global_malloc( A->numcols + THREADS, sizeof(int64_t));
  if( shtmp == NULL ) return(NULL);
  int64_t * l_shtmp =  shtmp;
  int64_t lnc = (A->numcols + THREADS - MYTHREAD - 1)/THREADS;
  for(i=0; i < lnc; i++)
    l_shtmp[i] = 0;
  shmem_barrier_all();

  for( i=0; i< A->lnnz; i++) {                   // histogram the column counts of A
    assert( A->lnonzero[i] < A->numcols );
    assert( A->lnonzero[i] >= 0 ); 
    int64_t index = A->lnonzero[i];
    int64_t lindex = index/shmem_n_pes();
    int64_t pe = index % shmem_n_pes();
    pos = shmem_int64_atomic_fetch_inc(shtmp+lindex, pe);
  }
  shmem_barrier_all();


  lnnz = 0;
  for( i = 0; i < lnc; i++) {
    lnnz += l_shtmp[i];
  }
  
  At = init_matrix(A->numcols, A->numrows, lnnz, 0);
  if(!At){printf("ERROR: transpose_matrix_upc: init_matrix failed!\n");return(NULL);}

  int64_t sum = myshmem_reduce_add(lnnz);      // check the histogram counted everything
  assert( A->nnz == sum ); 

  // compute the local offsets
  At->loffset[0] = 0;
  for(i = 1; i < At->lnumrows+1; i++)
    At->loffset[i] = At->loffset[i-1] + l_shtmp[i-1];

  // get the global indices of the start of each row of At
  for(i = 0; i < At->lnumrows; i++)
    l_shtmp[i] = MYTHREAD + THREADS * (At->loffset[i]);
    
  shmem_barrier_all();

  //redistribute the nonzeros 
  for(row=0; row<A->lnumrows; row++) {
    for(j=A->loffset[row]; j<A->loffset[row+1]; j++){
      int64_t index = A->lnonzero[j];
      int64_t lindex = index/shmem_n_pes();
      int64_t pe = index % shmem_n_pes();
      pos = shmem_int64_atomic_fetch_add(shtmp+lindex,  (int64_t) THREADS, pe);
      myshmem_int64_p(At->nonzero, pos, row*THREADS + MYTHREAD);
    }
  }
  shmem_barrier_all();

  //shmem_barrier_all();
  //if(!MYTHREAD)printf("done\n");
  shmem_free(shtmp);

  return(At);
}

/*! \brief comparison function to support qsort
 */
int nz_comp(const void *a, const void *b) {
  return( *(uint64_t *)a - *(uint64_t *)b );
}

/*! \brief sort the non-zeros in each row of a sparse matrix (make it tidy)
 * \param mat pointer to the sparse matrix
 * \ingroup spmatgrp
 */
int sort_nonzeros( sparsemat_t *mat) {
  int i;
  for(i = 0; i < mat->lnumrows; i++){
     qsort( &(mat->lnonzero[mat->loffset[i]]), mat->loffset[i+1] - mat->loffset[i], sizeof(int64_t), nz_comp );
  }
  shmem_barrier_all();
  return(0);
}

/*! \brief produce the transpose of a sparse matrix using UPC
 * \param omat  pointer to the original matrix
 * \return a pointer to the matrix that has be computed or NULL on failure
 *
 * This is wrapper for implementations written in the different models.
 * It is interest enough to be its own apps, one should experiment with it
 * within the apps framework. 
 *
 * \ingroup spmatgrp
 */
sparsemat_t * transpose_matrix(sparsemat_t *omat) {
  sparsemat_t * A;
  A = transpose_matrix_agi(omat);
  //A = transpose_matrix_exstack(omat, 1024);
  //A = transpose_matrix_exstack2(omat, 1024);
  //A = transpose_matrix_conveyor(omat);
  if(!A){return(NULL);}
  
  sort_nonzeros(A);
  return(A);
}

/*! \brief Generates the lower half of the adjacency matrix (non-local) for an Erdos-Renyi random
 * graph. This subroutine uses ALG1 from the paper "Efficient Generation of Large Random Networks" 
 * by Batageli and Brandes appearing in Physical Review 2005. 
 * Instead of flipping a coin for each potential edge this algorithm generates a sequence of 
 * "gaps" between 1s in the upper or lower triangular portion of the adjacency matrix using 
 * a geometric random variable.
 *
 * We parallelized this algorithm by noting that the geometric random variable is memory-less. This means, we 
 * can start the first row on each PE independently. In fact, we could start each row from scratch if we
 * we wanted to! This makes this routine embarrassingly parallel.
 *
 * \param n The total number of vertices in the graph.
 * \param p The probability that each non-loop edge is present.
 * \param edgetype See enum edge_type enum. DIRECTED, UNDIRECTED, DIRECTED_WEIGHTED, UNDIRECTED_WEIGHTED
 * \param loops See self_loops enum. Does all or no vertices have self loops.
 * \param seed A random seed. This should be a single across all PEs (it will be modified by each PE individually).
 * \return A distributed sparsemat_t that holds the graph.  It is either a lower triangular  matrix
           or square matrix, weighted or not with or without a diagonal.
 */
  sparsemat_t * erdos_renyi_random_graph(int64_t n, double p, edge_type edgetype, self_loops loops, uint64_t seed){
  
  int64_t row, col, i;
  int64_t ln = (n + THREADS - MYTHREAD - 1)/THREADS;
  int64_t lnnz = 0, lnnz_orig=0;
  double D = log(1 - p);
  int64_t r;
  int64_t end = n;
  int64_t ndiag = ln;
  
  rand_seed(seed);

  /* count lnnz so we can allocate A correctly */   
  row = MYTHREAD;
  col = 1 + floor(log(1 - rand_double()) / D);
  while(row < n){
    if(edgetype == UNDIRECTED || edgetype == UNDIRECTED_WEIGHTED)
      end = row;
    while(col < end){
      // if we just hit a diagonal entry (we don't have to generate this one later)
      if(col == row) ndiag--; 
      lnnz_orig++;
      col += 1 + floor(log(1 - rand_double()) / D);
    }

    row += THREADS;
    col -= end;
  }
  if(loops == LOOPS) lnnz_orig += ndiag;

  shmem_barrier_all();

  int weighted = (edgetype == UNDIRECTED_WEIGHTED || edgetype == DIRECTED_WEIGHTED);
  sparsemat_t * A = init_matrix(n, n, lnnz_orig, weighted);
  if(!A){T0_printf("ERROR: erdos_renyi_random_graph: init_matrix failed!\n"); return(NULL);}

  /* reset the seed so we get the same sequence of coin flips */
  rand_seed(seed);

  /* now go through the same sequence of random events and fill in nonzeros */
  A->loffset[0] = 0;
  row = MYTHREAD;
  col = 1 + floor(log(1 - rand_double()) / D);
  while(row < n){
    int need_diag = (loops == LOOPS);
    if(edgetype == UNDIRECTED || edgetype == UNDIRECTED_WEIGHTED)
      end = row;
    while(col < end){
      if(col == row) need_diag = 0;
      A->lnonzero[lnnz++] = col;
      col += 1 + floor(log(1 - rand_double()) / D);
    }
    if(need_diag) {
      A->lnonzero[lnnz++] = row;
    }
    row+=THREADS;
    col -= end;
    A->loffset[row/THREADS] = lnnz;
  }
  
  if(lnnz != lnnz_orig){
    printf("ERROR: lnnz (%ld) != lnnz_orig (%ld)\n", lnnz, lnnz_orig);
    return(NULL);
  }

  // fill in the weights
  if(weighted){
    for(i = 0; i < lnnz; i++){
      A->lvalue[i] = rand_double();
    }
  }

  if(loops == LOOPS && (edgetype == DIRECTED || edgetype == DIRECTED_WEIGHTED))
    sort_nonzeros(A); // to get the diagonal entry sorted correctly
  
  return(A);
}

/*!
  \brief Compute partial sums across threads.  
  In the formulas below, \a m represents <tt>MYTHREAD</tt>.
  \note This function must be called on all threads.
  \param x input value \f$x_m\f$
  \return \f$\sum_{i<=m} x_i\f$ 
* \ingroup libgetputgrp
*/
int64_t myshmem_partial_add(int64_t x) {

  int64_t * tmp = (int64_t*)myshmem_global_malloc(THREADS, sizeof(int64_t));
  int64_t out = 0;
  
  myshmem_int64_p(tmp, MYTHREAD, x);

  shmem_barrier_all();

  for (int i = 0; i <= MYTHREAD; i++) {
    out += myshmem_int64_g(tmp, i);
  }

  shmem_barrier_all();

  return out;
}

/*! 
  \brief Compute prior partial sums (not including this value) across threads.  
  In the formulas below, \a m represents <tt>MYTHREAD</tt>.
  \note This function must be called on all threads.
  \param x input value \f$x_m\f$
  \return \f$\sum_{i<m} x_i\f$ 
  * \ingroup libgetputgrp
  */
int64_t myshmem_prior_add(int64_t x) {
  return myshmem_partial_add(x) - x;
}

/*! \brief checks that a global array is in fact a permutation
 * \param perm SHARED pointer to the global array
 * \param N the length of the global array
 * \return 1 if it is permutation
 * \ingroup spmatgrp
 */
int is_perm(int64_t * perm, int64_t N) {
  int64_t i;
  int64_t * lperm = perm;
  int64_t * flag = (int64_t*)myshmem_global_malloc(N, sizeof(int64_t));
  if( flag == NULL ) return(0);
  int64_t * lflag = flag;
  int64_t l_N = (N + THREADS - MYTHREAD - 1)/THREADS;

  for(i = 0; i < l_N; i++) 
    lflag[i] = 0;
  shmem_barrier_all();
  for(i = 0; i < l_N; i++) 
   myshmem_int64_p(flag, lperm[i], 1);
  shmem_barrier_all();
  int64_t err = 0L;
  for(i = 0; i < l_N; i++) 
    if(lflag[i] == 0) 
      err++;
  shmem_free(flag);
  err = myshmem_reduce_add(err);
  return(err == 0L);
}

/******************************************************************************/
/*! \brief create a global int64_t array with a uniform random permutation
 * \param N the length of the global array
 * \param seed seed for the random number generator
 * \return the permutation
 * 
 * This is a collective call.
 * this implements the random dart algorithm to generate the permutation.
 * Each thread throws its elements of the perm array randomly at large target array.
 * Each element claims a unique entry in the large array using compare_and_swap.
 * This gives a random permutation with spaces in it, then you squeeze out the spaces.
 * \ingroup spmatgrp
 */
int64_t * rand_permp_agi(int64_t N, int seed) {  
  int64_t * ltarget, *lperm;
  int64_t r, i, j;
  int64_t pos, numdarts, numtargets, lnumtargets;

  if( seed != 0 ) srand48( seed );

  //T0_printf("Entering rand_permp_atomic...");fflush(0);

  int64_t * perm = (int64_t*)myshmem_global_malloc(N, sizeof(int64_t));
  if( perm == NULL ) return(NULL);
  lperm = perm;

  int64_t l_N = (N + THREADS - MYTHREAD - 1)/THREADS;
  int64_t M = 2*N;
  int64_t l_M = (M + THREADS - MYTHREAD - 1)/THREADS;

  int64_t * target = (int64_t*)myshmem_global_malloc(M, sizeof(int64_t));
  if( target == NULL ) return(NULL);
  ltarget = target;
  
  for(i=0; i<l_M; i++)
    ltarget[i] = -1L;
  shmem_barrier_all();

  i=0;
  while(i < l_N){                // throw the darts until you get l_N hits
    r = lrand48() % M;
    long lindex = r/shmem_n_pes();
    long pe = r % shmem_n_pes();
    if( shmem_long_atomic_compare_swap(target+lindex, -1L, (i*THREADS + MYTHREAD), pe) == (-1L) ){
      i++;
    }
  }
  shmem_barrier_all();

  numdarts = 0;
  for(i = 0; i < l_M; i++)    // count how many darts I got
    numdarts += (ltarget[i] != -1L );

  pos = myshmem_prior_add(numdarts);    // my first index in the perm array is the number 
                                      // of elements produce by the smaller threads
  for(i = 0; i < l_M; i++){
    if(ltarget[i] != -1L ) {
       myshmem_int64_p(perm, pos, ltarget[i]);
       pos++;
    }
  }

  shmem_free(target);
  shmem_barrier_all();
  //T0_printf("done!\n");
  return(perm);
}

/*! \brief Produce a global array the holds a uniform random permutation.
 * \param N the global length of the permutaion
 * \param seed the seed for the random number generator
 * \return a pointer to the matrix that has been produced or NULL if the model can't be used
 *
 * This is wrapper for implementations written in the different models.
 * It is interest enough to be its own apps, one should experiment with it
 * within the apps framework. 
 *
 * This is a collective call.
 *
 * this implements the random dart algorithm to generate the permutation.
 * Each thread throws its elements of the perm array randomly at large target array.
 * Each element claims a unique entry in the large array using compare_and_swap.
 * This gives a random permutation with spaces in it, then you squeeze out the spaces.
 * \ingroup spmatgrp
 */
int64_t * rand_permp(int64_t N, int seed) {
  int64_t * p;
  p = rand_permp_agi(N, seed);
  //p = rand_permp_exstack(N, seed, 1024);
  //p = rand_permp_exstack2(N, seed, 1024);
  
  if(!is_perm(p, N)){
    T0_printf("ERROR: rand_permp: not a permutation!\n");fflush(0);
    return(NULL);
  }
  return(p);
}

/*! \brief apply row and column permutations to a sparse matrix using straight UPC
 * \param A pointer to the original matrix
 * \param rperminv pointer to the global array holding the inverse of the row permutation
 * \param cperminv pointer to the global array holding the inverse of the column permutation
 * rperminv[i] = j means that row i of A goes to row j in matrix Ap
 * cperminv[i] = j means that col i of A goes to col j in matrix Ap
 * \return a pointer to the matrix that has been produced or NULL if the model can't be used
 * \ingroup spmatgrp
 */
sparsemat_t * permute_matrix_agi(sparsemat_t *A, int64_t *rperminv, int64_t *cperminv) {
  //T0_printf("Permuting matrix with single puts\n");
  int64_t i, j, col, row, pos;
  int64_t * lrperminv = rperminv;
  int64_t * rperm = (int64_t*)myshmem_global_malloc(A->numrows, sizeof(int64_t));
  if( rperm == NULL ) return(NULL);
  int64_t *lrperm = rperm;

  //compute rperm from rperminv 
  for(i=0; i < A->lnumrows; i++){
    myshmem_int64_p(rperm, lrperminv[i], i*THREADS + MYTHREAD);
  }

  shmem_barrier_all();
  
  int64_t cnt = 0, off, nxtoff;
  for(i = 0; i < A->lnumrows; i++){
    row = lrperm[i];
    off    = myshmem_int64_g(A->offset, row);
    nxtoff = myshmem_int64_g(A->offset, row + THREADS);
    cnt += nxtoff - off;
  }
  shmem_barrier_all();

  sparsemat_t * Ap = init_matrix(A->numrows, A->numcols, cnt, 0);
  
  // fill in permuted rows
  Ap->loffset[0] = pos = 0;
  for(i = 0; i < Ap->lnumrows; i++){
    row = lrperm[i];
    off    = myshmem_int64_g(A->offset, row);
    nxtoff = myshmem_int64_g(A->offset, row + THREADS);
    for(j = off; j < nxtoff; j++){
      Ap->lnonzero[pos++] = myshmem_int64_g(A->nonzero, j*THREADS + row%THREADS);
    }
    Ap->loffset[i+1] = pos;
  }
  
  assert(pos == cnt);

  shmem_barrier_all();
  
  // finally permute column indices
  for(i = 0; i < Ap->lnumrows; i++){
    for(j = Ap->loffset[i]; j < Ap->loffset[i+1]; j++){
      Ap->lnonzero[j] = myshmem_int64_g(cperminv, Ap->lnonzero[j]);      
    }
  }
  shmem_barrier_all();

  shmem_free(rperm);

  T0_printf("done\n");
  return(Ap);
}

/*! \brief apply row and column permutations to a sparse matrix using straight UPC
 * \param omat pointer to the original matrix
 * \param rperminv pointer to the global array holding the inverse of the row permutation
 * \param cperminv pointer to the global array holding the inverse of the column permutation
 * rperminv[i] = j means that row i of A goes to row j in matrix Ap
 * cperminv[i] = j means that col i of A goes to col j in matrix Ap
 * \param mode which buffering model to use
 * \return a pointer to the matrix that has be computed or NULL on failure
 *
 * This is wrapper for implementations written in the different models.
 * It is interest enough to be its own apps, one should experiment with it
 * within the apps framework.
 *
 * \ingroup spmatgrp
 */
sparsemat_t * permute_matrix(sparsemat_t *omat, int64_t *rperminv, int64_t *cperminv) {
  return( permute_matrix_agi(omat, rperminv, cperminv) );
  // return( permute_matrix_exstack(omat, rperminv, cperminv, 1024) );
  //return( permute_matrix_exstack2(omat, rperminv, cperminv, 1024) );
  //return( permute_matrix_conveyor(omat, rperminv, cperminv) );
}

/*! \brief checks that the sparse matrix is lower triangluar
 * \param A pointer to the sparse matrix
 * \param unit_diagonal set to 1 to make sure all pivots are nonzero
 * \return 0 on success, non-0 on error.
 * kind of a speciality routine to check that toposort might of worked
 * \ingroup spmatgrp
 */
int is_lower_triangular(sparsemat_t *A, int64_t unit_diagonal) {
  int64_t i,j, row, * ltri_rperm, * ltri_cperm;
  int64_t err = 0, err2 = 0;

  shmem_barrier_all();

  /* we are hoping for an lower triangular matrix here */
  for(i=0; i < A->lnumrows; i++){
    int64_t global_row = i*THREADS + MYTHREAD;
    int pivot = 0;
    for(j=A->loffset[i]; j < A->loffset[i+1];j++){
      if( A->lnonzero[j] > global_row ) {
        err++;
      }else if( A->lnonzero[j] == global_row){
        pivot = 1;
      }
    }
    if(!pivot)
      err2++;
  }

  err = myshmem_reduce_add(err);
  err2 = (unit_diagonal ? myshmem_reduce_add(err2) : 0);
  if( err || err2 ){
    //if(!MYTHREAD)printf("\nThere are %ld nz above diag. and %ld missing pivots in lower.\n", err, err2);
    fflush(0);
  }

  shmem_barrier_all();

  return(!(err || err2));
}

/*! \brief checks that the sparse matrix is upper triangluar
 * \param A pointer to the sparse matrix
 * \param unit_diagonal set to 1 to make sure all pivots are nonzero
 * \return 0 on success, non-0 on error.
 *
 * \ingroup spmatgrp
 */
int is_upper_triangular(sparsemat_t *A, int64_t unit_diagonal) {
  int64_t i,j, row, * ltri_rperm, * ltri_cperm;
  int64_t err = 0, err2 = 0;
  shmem_barrier_all();

  /* we are hoping for an upper triangular matrix here */
  for(i=0; i < A->lnumrows; i++){
    int64_t global_row = i*THREADS + MYTHREAD;
    int pivot = 0;
    for(j=A->loffset[i]; j < A->loffset[i+1];j++){
      if( A->lnonzero[j] < global_row ) {
        err++;
      }else if( A->lnonzero[j] == global_row){
        pivot = 1;
      }
    }
    if(!pivot)
      err2++;
  }

  err = myshmem_reduce_add(err);
  err2 = (unit_diagonal ? myshmem_reduce_add(err2) : 0);
  if( err || err2 ){
    //if(!MYTHREAD)printf("\nThere are %ld nz below diag. and %ld missing pivots in upper.\n", err, err2);
    fflush(0);
  }

  shmem_barrier_all();

  return(!(err || err2));
}

/*! \brief compare the structs that hold two sparse matrices
 * \param lmat pointer to the left sparse matrix
 * \param rmat pointer to the right sparse matrix
 * \return 0 on success
 * \ingroup spmatgrp
 */
int compare_matrix(sparsemat_t *lmat, sparsemat_t *rmat) {
  int i,j;

  if( lmat->numrows != rmat->numrows ){
    if(!MYTHREAD)printf("(lmat->numrows = %ld)  != (rmat->numrows = %ld)", lmat->numrows, rmat->numrows );
    return(1);
  }
  if( lmat->lnumrows != rmat->lnumrows ){
    fprintf(stderr,"THREAD %03d: (lmat->lnumrows = %ld)  != (rmat->lnumrows = %ld)", 
            MYTHREAD, lmat->lnumrows, rmat->lnumrows );
    return(1);
  }
  if( lmat->numcols != rmat->numcols ){
    if(!MYTHREAD)printf("(lmat->numcols = %ld)  != (rmat->numcols = %ld)", lmat->numcols, rmat->numcols );
    return(1);
  }
  if( lmat->nnz != rmat->nnz ){
    if(!MYTHREAD)printf("(lmat->nnz = %ld)  != (rmat->nnz = %ld)", lmat->nnz, rmat->nnz );
    return(1);
  }
  if( lmat->lnnz != rmat->lnnz ){
    fprintf(stderr,"THREAD %03d: (lmat->lnnz = %ld)  != (rmat->lnnz = %ld)", 
            MYTHREAD, lmat->lnnz, rmat->lnnz );
    return(1);
  }

  if( lmat->loffset[0] != 0 || rmat->loffset[0] != 0 
    || (lmat->loffset[0] != rmat->loffset[0] ) ){
    if(!MYTHREAD)printf("THREAD %03d: (lmat->loffset[0] = %ld)  != (rmat->loffset[0] = %ld)", 
       MYTHREAD, lmat->loffset[0], rmat->loffset[0] );
    return(1);
  }

  
  for(i = 0; i < lmat->lnumrows; i++){
    if( lmat->loffset[i+1] != rmat->loffset[i+1] ){
       if(!MYTHREAD)printf("THREAD %03d: (lmat->loffset[%d] = %ld)  != (rmat->loffset[%d] = %ld)", 
          MYTHREAD, i+1, lmat->loffset[i+1], i+1, rmat->loffset[i+1] );
       return(1);
    }
  }
  
  for(j=0; j< lmat->lnnz; j++) {
    if( lmat->lnonzero[j] != rmat->lnonzero[j] ){
      if(!MYTHREAD)printf("THREAD %03d: (lmat->lnonzero[%d] = %ld)  != (rmat->lnonzero[%d] = %ld)", 
                MYTHREAD, j, lmat->lnonzero[j], j, rmat->lnonzero[j] );
      return(1);
    }
  }

  return(0);
}

#endif //MYSHMEM_H

