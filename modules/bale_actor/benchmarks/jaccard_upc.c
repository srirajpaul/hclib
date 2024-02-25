/******************************************************************
//
// 
 *****************************************************************/ 
/*! \file jaccard_upc.c
 * \brief Demo application that calculates jaccard index for edges in a graph. 
 *      The Jaccard index gives a value that represents the similarity 
 *      of the neighborhoods of two vertices connected by an edge.
 *
 *      NOTE: this program is created for full symmetric L+U matrices
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

double get_edge_degrees(sparsemat_t* L, int64_t* edge_degrees) {
    int64_t pos = 0;
    for (int64_t y = 0; y < L->lnumrows; y++) {
        for (int64_t k = L->loffset[y]; k < L->loffset[y+1]; k++) {
            if ((y*THREADS + MYTHREAD) < L->lnonzero[k]) continue;
            int64_t deg_u_upper = L->offset[L->lnonzero[k] + THREADS];
            int64_t deg_u_lower = L->offset[L->lnonzero[k]];
            edge_degrees[pos] = deg_u_upper - deg_u_lower;
            pos++;
        }
    }
    upc_barrier;

    return 0;
}

double jaccard_agi(sparsemat_t* A2, double* jaccard_index, int64_t* edge_degrees) {
    // Start timing
    double t1 = wall_seconds();
    upc_barrier;

    int64_t pos = 0;
    for (int64_t v = 0; v < A2->lnumrows; v++) {             //vertex v (local)
        int64_t deg_v = A2->loffset[v+1] - A2->loffset[v];
        for (int64_t k = A2->loffset[v]; k < A2->loffset[v+1]; k++) {
            int64_t u = A2->lnonzero[k];                     //vertex u (possibly remote)
            if ((v*THREADS + MYTHREAD) < u) continue;

            // calculate intersection
            for (int64_t kk = A2->loffset[v]; kk < A2->loffset[v+1]; kk++) {
                int64_t x = A2->lnonzero[kk];
                if (x == u) continue;
                int64_t row_begin = A2->offset[u];
                int64_t row_end = A2->offset[u + THREADS];
                for (int64_t uk = row_begin; uk < row_end; uk++) {
                    int64_t nonzero_x = A2->nonzero[u%THREADS + uk*THREADS];
                    if (x == nonzero_x) {
                        jaccard_index[pos]++;
                    }
                }
            }

            pos++;
        }

    }
    upc_barrier;

    pos = 0;
    for (int64_t v = 0; v < A2->lnumrows; v++) {
        int64_t deg_v = A2->loffset[v+1] - A2->loffset[v];
        for (int64_t k = A2->loffset[v]; k < A2->loffset[v+1]; k++) {
            if ((v*THREADS + MYTHREAD) < A2->lnonzero[k]) continue;

            // calculate jaccard index using: intersection / (deg(u) + deg(v) - intersection)
            jaccard_index[pos] = jaccard_index[pos] / ((double)edge_degrees[pos] + (double)deg_v - jaccard_index[pos]);
            
            //increase position in jaccard index local array
            pos++;
        }
    }

    upc_barrier;
    t1 = wall_seconds() - t1;
    return t1;
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
    T0_fprintf(stderr,"ERROR: triangles: error generating input!\n"); upc_global_exit(1);
  }
  
  T0_fprintf(stderr,"B has %ld rows/cols and %ld nnz\n", B->numrows, B->lnnz);
  T0_fprintf(stderr,"C has %ld rows/cols and %ld nnz\n", C->numrows, C->lnnz);
  
  sparsemat_t * A = kron_prod_dist(B, C, 1);
  
  return(A);
}
*/

int main(int argc, char* argv[]) {

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
    while ((opt = getopt(argc, argv, "hM:n:f:a:e:K:")) != -1) {
        switch (opt) {
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

    upc_atomic_domain = upc_all_atomicdomain_alloc(UPC_INT64, UPC_ADD|UPC_INC|UPC_CSWAP, 0);
    int64_t numrows = l_numrows * THREADS;
    if (erdos_renyi_prob == 0.0) { // use nz_per_row to get erdos_renyi_prob
        erdos_renyi_prob = (2.0 * (nz_per_row - 1)) / numrows;
        if (erdos_renyi_prob > 1.0) erdos_renyi_prob = 1.0;
    } else {                     // use erdos_renyi_prob to get nz_per_row
        nz_per_row = erdos_renyi_prob * numrows;
    }
  
    T0_fprintf(stderr,"Running jaccard on %d threads\n", THREADS);
    if (!read_graph && !gen_kron_graph) {
        T0_fprintf(stderr,"Number of rows per thread   (-N)   %ld\n", l_numrows);
        T0_fprintf(stderr,"Erdos Renyi prob (-e)   %g\n", erdos_renyi_prob);
    }
    T0_fprintf(stderr,"Model mask (M) = %ld (should be 1,2,4,8,16 for agi, exstack, exstack2, conveyors, alternates\n", models_mask);  
    T0_fprintf(stderr,"algorithm (a) = %ld (0 for L & L*U, 1 for L & U*L)\n", alg);
  
    if (printhelp)
        upc_global_exit(0);

    // alg = 0 only needs L
    // alg = 1 needs both U and L

    double correct_answer = -1;
  
    sparsemat_t *A, *L, *U, *A2;
/*    
    if (read_graph) {
        A = read_matrix_mm_to_dist(filename);
        if (!A) upc_global_exit(1);

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
        if(!U){T0_printf("ERROR: gen_er_graph_dist: U is NULL!\n"); upc_global_exit(1);}
        int64_t lnnz = L->lnnz + U->lnnz;
        A2 = init_matrix(n, n, lnnz, 0);
        if(!A2){T0_printf("ERROR: gen_er_graph_dist: A2 is NULL!\n"); upc_global_exit(1);}

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
        upc_barrier;
        clear_matrix(U);
        clear_matrix(L);
        free(U); free(L);
        // A2 is a full symmetric adjacency matrix
    //}

    upc_barrier;

    T0_fprintf(stderr, "A2 has %ld rows/cols and %ld nonzeros.\n", A2->numrows, A2->nnz);

    T0_fprintf(stderr, "\nRunning Jaccard Index (upc): \n");

    int64_t ljaccard_size = myupc_reduce_max_l(A2->lnnz);
    double* jaccard_index = calloc(ljaccard_size, sizeof(double));
    for (int64_t x = 0; x < ljaccard_size; x++) jaccard_index[x] = 0;
    // get the degrees of each vertex
    int64_t* edge_degrees = calloc(ljaccard_size, sizeof(int64_t));
    for (int64_t x = 0; x < ljaccard_size; x++) edge_degrees[x] = 0;
    upc_barrier;
    get_edge_degrees(A2, edge_degrees);
    
    upc_barrier;
    //for (int64_t i = 0; i < ljaccard_size; i++) printf("\npe: %d, deg: %d", MYTHREAD, edge_degrees[i]);
    double laptime_jaccard = 0.0;

    // Running agi model for jaccard index
    laptime_jaccard = jaccard_agi(A2, jaccard_index, edge_degrees);
    upc_barrier;

    //for (int64_t i = 0; i < ljaccard_size; i++) printf("\npe: %d, jaccard index: %f", MYTHREAD, jaccard_index[i]);
    //printf("\n");

    T0_fprintf(stderr, " %8.5lf seconds\n", laptime_jaccard);

    upc_barrier;

    return 0;

}

