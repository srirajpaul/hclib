/******************************************************************
//
// 
 *****************************************************************/ 
/*! \file tricent_agi.cpp
 * \brief Demo application that calculates triangle centrality for vertices in a graph. 
 *      Triangle centrality involves finding important vertices in graphs according to the 
 *      concentration of triangles in the neighborhood of a vertex.
 *      
 *      NOTE: this program is created for full symmetric L+U matrices
 */

#include <math.h>
#include <shmem.h>
extern "C" {
#include <spmat.h>
}

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
  lgp_global_exit(0);
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
/*
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
*/

double tricount_agi(int64_t *count, sparsemat_t * L, int64_t* vertex_tri_count, int64_t* lvertex_tri_count) {
  int64_t cnt=0;
  int64_t l_i, ii, k, kk, w, L_i, L_j;
   
  double t1 = wall_seconds();

  //foreach nonzero (i, j) in L
  for(l_i = 0; l_i < L->lnumrows; l_i++){ 
    for(k = L->loffset[l_i] + 1; k < L->loffset[l_i + 1]; k++) {
      L_i = l_i*THREADS + MYTHREAD;

      L_j = L->lnonzero[k]; if (L_j >= L->numrows) L_j = L_j - L->numrows; // conversion due to negating triangle neighbor nnz's
      if (L_i < L_j) continue;  // add check due to full matrix

      // NOW: get L[L_j,:] and count intersections with L[L_i,:]
      int64_t start = L->loffset[l_i];
      int64_t kbegin = lgp_get_int64(L->offset, L_j); 
      int64_t kend   = lgp_get_int64(L->offset, L_j + THREADS);
      for(kk = kbegin; kk < kend; kk++){
        w = lgp_get_int64(L->nonzero, L_j%THREADS + kk*THREADS); if (w >= L->numrows) w = w - L->numrows; // conversion due to negating triangle neighbor nnz's
        if (L_j < w) continue;  // add check due to full matrix

        for(ii = start; ii < L->loffset[l_i + 1]; ii++) {
          int64_t nonzero_ii = L->lnonzero[ii]; if (nonzero_ii >= L->numrows) nonzero_ii = nonzero_ii - L->numrows; // conversion due to negating triangle neighbor nnz's
          if( w == nonzero_ii ){ 
            cnt++;

            lvertex_tri_count[l_i]++; //vertex u is local
            lgp_fetch_and_inc(vertex_tri_count, L_j + L_j%THREADS); //vertex v is possibly remote 
            lgp_fetch_and_inc(vertex_tri_count, w + w%THREADS); // vertex w is possibly remote

            // update tri_N(u)
            // loop through u's nonzeros, first_neigh is w, second_neigh is v
            for (int64_t bu = L->loffset[l_i]; bu < L->loffset[l_i + 1]; bu++) {
              if (w == L->lnonzero[bu] || L_j == L->lnonzero[bu]) {
                L->lnonzero[bu] = L->lnonzero[bu] + L->numrows;
              }
            }

            // update tri_N(v)
            // loop through v's nonzeros, first_neigh is u, second_neigh is w
            for(int64_t bv = kbegin; bv < kend; bv++){
              int64_t nnz_v = lgp_get_int64(L->nonzero, L_j%THREADS + bv*THREADS);
              if ((l_i * THREADS + MYTHREAD) == nnz_v || w == nnz_v) {
                lgp_fetch_and_add(L->nonzero, L_j%THREADS + bv*THREADS, L->numrows);
              }
            }

            // update tri_N(w)
            // loop through w's nonzeros, first_neigh is u, second_neigh is v
            int64_t wbegin = lgp_get_int64(L->offset, w); 
            int64_t wend   = lgp_get_int64(L->offset, w + THREADS);
            for(int64_t bw = wbegin; bw < wend; bw++){
              int64_t nnz_w = lgp_get_int64(L->nonzero, w%THREADS + bw*THREADS);
              if ((l_i * THREADS + MYTHREAD) == nnz_w || L_j == nnz_w) {
                lgp_fetch_and_add(L->nonzero, w%THREADS + bw*THREADS, L->numrows);
              }
            }

            start = ii + 1;
            break;
          }
          if( w < nonzero_ii ){ // the rest are all bigger because L is tidy
            start = ii;
            break;
          }
        }
      }
    }
  }

  lgp_barrier();
  t1 = wall_seconds() - t1;
  *count = cnt;
  return(t1);
}

double tricent_agi(int64_t total_tri_cnt, double* tricent_counts, sparsemat_t* L, int64_t* vertex_tri_count, int64_t* lvertex_tri_count) {
  // Start timing
  double t2 = wall_seconds();
  lgp_barrier();

  for (int64_t v = 0; v < L->lnumrows; v++) {
    for (int64_t k = L->loffset[v]; k < L->loffset[v+1]; k++) {
      int64_t u = L->lnonzero[k];
      if (u >= L->numrows) {
        // nnz is in tri_N(v) set, update (1/3)*T(u) / T(G)
        u = u - L->numrows;
        int64_t tri_count_neigh = lgp_get_int64(vertex_tri_count, u);
        tricent_counts[v] += (((double)1/3) * ((double)tri_count_neigh)) / (double)total_tri_cnt;
      } else {
        // nnz is in N(v)\tri_N(v), update (1)*T(u) / T(G)
        int64_t tri_count_non_neigh = lgp_get_int64(vertex_tri_count, u);
        tricent_counts[v] += ((double)tri_count_non_neigh) / (double)total_tri_cnt;
      }
    }
    tricent_counts[v] += (((double)1/3) * ((double)lvertex_tri_count[v])) / (double)total_tri_cnt;
  }

  lgp_barrier();
  t2 = wall_seconds() - t2;
  return t2;
}

int main(int argc, char * argv[]) {

  lgp_init(argc, argv);

  int64_t models_mask = 0;//ALL_Models;  // default is running all models
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

  int64_t numrows = l_numrows * THREADS;
  if(erdos_renyi_prob == 0.0){ // use nz_per_row to get erdos_renyi_prob
    erdos_renyi_prob = (2.0*(nz_per_row - 1))/numrows;
    if(erdos_renyi_prob > 1.0)
      erdos_renyi_prob = 1.0;
  } else {                     // use erdos_renyi_prob to get nz_per_row
    nz_per_row = erdos_renyi_prob * numrows;
  }
  
  T0_fprintf(stderr,"Running triangle centrality on %d PEs\n", THREADS);
  if(!read_graph && !gen_kron_graph){
    T0_fprintf(stderr,"Number of rows per PE   (-N)   %ld\n", l_numrows);
    T0_fprintf(stderr,"Erdos Renyi prob (-e)   %g\n", erdos_renyi_prob);
  }
  T0_fprintf(stderr,"Model mask (M) = %ld (should be 1,2,4,8,16 for agi, exstack, exstack2, conveyors, alternates\n", models_mask);  
  T0_fprintf(stderr,"algorithm (a) = %ld (0 for L & L*U, 1 for L & U*L)\n", alg);
  
  if( printhelp )
    lgp_global_exit(0);

  // alg = 0 only needs L
  // alg = 1 needs both U and L

  double correct_answer = -1;
  
    sparsemat_t *A, *L, *U, *A2;
/*    
    if (read_graph) {
        A = read_matrix_mm_to_dist(filename);
        if (!A) lgp_global_exit(1);

        T0_fprintf(stderr,"Reading file %s...\n", filename);
        T0_fprintf(stderr, "A has %ld rows/cols and %ld nonzeros.\n", A->numrows, A->nnz);

        // we should check that A is symmetric!
    
        if (!is_lower_triangular(A, 0)) { //if A is not lower triangular... make it so.      
            T0_fprintf(stderr, "Assuming symmetric matrix... using lower-triangular portion...\n");
            tril(A, -1);
            L = A;
        } else {
            L = A;
        }
    
        sort_nonzeros(L);

    }  else if (gen_kron_graph) {
        // string should be <mode> # # ... #
        // we will break the string of numbers (#s) into two groups and create
        // two local kronecker graphs out of them.
        int num;
        char* ptr = kron_graph_string;
        int64_t* kron_specs = (int64_t*)calloc(32, sizeof(int64_t*));
    
        // read the mode
        int ret = sscanf(ptr, "%d ", &kron_graph_mode);
        if (ret == 0) ret = sscanf(ptr, "\"%d ", &kron_graph_mode);
        if (ret == 0) { T0_fprintf(stderr, "ERROR reading kron graph string!\n"); return(1); }
        T0_fprintf(stderr,"kron string: %s return = %d\n", ptr, ret);
        T0_fprintf(stderr,"kron mode: %d\n", kron_graph_mode);
        ptr += 2;
        int mat, num_ints = 0;
        while (sscanf(ptr, "%d %n", &num, &mat) == 1) {
            T0_fprintf(stderr,"%s %d\n", ptr, mat);
            kron_specs[num_ints++] = num;
            ptr+=mat;
        }

        if (num_ints <= 1) {
            T0_fprintf(stderr,"ERROR: invalid kronecker product string (%s): must contain at least three integers\n", kron_graph_string); return(-1);
        }

        // calculate the number of triangles 
        if (kron_graph_mode == 0) {
            correct_answer = 0.0;
        } else if (kron_graph_mode == 1) {
            correct_answer = 1;
            for (i = 0; i < num_ints; i++)
                correct_answer *= (3 * kron_specs[i] + 1);
            correct_answer *= 1.0 / 6.0;
            double x = 1;
            for (i = 0; i < num_ints; i++)
                x *= (kron_specs[i] + 1);
            correct_answer = correct_answer - 0.5 * x + 1.0 / 3.0;
        } else if (kron_graph_mode == 2) {
            correct_answer = (1.0 / 6.0) * pow(4, num_ints) - pow(2.0, (num_ints - 1)) + 1.0 / 3.0;
        }

        correct_answer = round(correct_answer);
        T0_fprintf(stderr, "Pre-calculated answer = %ld\n", (int64_t)correct_answer);

        int64_t half = num_ints / 2;

        L = generate_kronecker_graph(kron_specs, half, &kron_specs[half], num_ints - half, kron_graph_mode);
    } else {
*/
        int n = numrows;
        L = erdos_renyi_random_graph(numrows, erdos_renyi_prob, UNDIRECTED, NOLOOPS, 12345);
        U = transpose_matrix(L);
        if(!U){T0_printf("ERROR: gen_er_graph_dist: U is NULL!\n"); lgp_global_exit(1);}
        int64_t lnnz = L->lnnz + U->lnnz;
        A2 = init_matrix(n, n, lnnz, 0);
        if(!A2){T0_printf("ERROR: gen_er_graph_dist: A2 is NULL!\n"); lgp_global_exit(1);}

        A2->loffset[0] = 0;
        lnnz = 0;
        for(i = 0; i < L->lnumrows; i++){
            int64_t global_row = i*THREADS + MYTHREAD;
            for(j = L->loffset[i]; j < L->loffset[i+1]; j++)
                A2->lnonzero[lnnz++] = L->lnonzero[j];
            for(j = U->loffset[i]; j < U->loffset[i+1]; j++){
                if(U->lnonzero[j] != global_row) // don't copy the diagonal element twice!
                    A2->lnonzero[lnnz++] = U->lnonzero[j];
            }
            A2->loffset[i+1] = lnnz;
        }
        lgp_barrier();
        clear_matrix(U);
        clear_matrix(L);
        free(U); free(L);
        // A2 is a full symmetric adjacency matrix
    //}

  lgp_barrier();

  T0_fprintf(stderr, "A2 has %ld rows/cols and %ld nonzeros.\n", A2->numrows, A2->nnz);

  T0_fprintf(stderr, "Run triangle centrality ...\n");
  int64_t tri_cnt = 0;           // partial count of triangles on this thread
  int64_t total_tri_cnt = 0;     // the total number of triangles on all threads
  double laptime_tricount = 0.0;
  double laptime_tricent = 0.0;
  int64_t ltc_size = A2->lnumrows;
  int64_t tc_size = ltc_size * THREADS;
  int64_t* vertex_tri_count = (int64_t*)lgp_all_alloc(tc_size, sizeof(int64_t));
  int64_t* lvertex_tri_count = lgp_local_part(int64_t, vertex_tri_count);
  for (int64_t x = 0; x < ltc_size; x++) lvertex_tri_count[x] = 0;
  double* tricent_counts = (double*)calloc(ltc_size, sizeof(double));
  for (int64_t x = 0; x < ltc_size; x++) tricent_counts[x] = 0.0;

  // Running agi model for triangle counting
  T0_fprintf(stderr, "\nrunning tricount (agi) ... \n");
  laptime_tricount = tricount_agi(&tri_cnt, A2, vertex_tri_count, lvertex_tri_count);
  lgp_barrier();
  total_tri_cnt = lgp_reduce_add_l(tri_cnt);
  lgp_barrier();
  //T0_fprintf(stderr, "  %8.5lf seconds: %16ld triangles\n", laptime_tricount, total_tri_cnt);
  
  // Running agi model for triangle centrality
  T0_fprintf(stderr, "\nrunning tricent (agi) ... \n");
  laptime_tricent = tricent_agi(total_tri_cnt, tricent_counts, A2, vertex_tri_count, lvertex_tri_count);
  lgp_barrier();
  T0_fprintf(stderr, "  %8.5lf seconds\n", laptime_tricount+laptime_tricent);
  //for (int64_t x = 0; x < ltc_size; x++) printf("tricent: %f\n",tricent_counts[x]);
  lgp_barrier();

  lgp_finalize();

  return(0);
}

