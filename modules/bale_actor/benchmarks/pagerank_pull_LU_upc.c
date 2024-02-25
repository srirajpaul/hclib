/******************************************************************
//
// 
 *****************************************************************/ 
/*! \file pagerank_upc.c
 * \brief Demo application that calculates page rank for vertices in a graph. 
 *      PageRank works by counting the number and quality of links to a page to 
 *      determine a rough estimate of how important the website is. The underlying 
 *      assumption is that more important websites are likely to receive more links 
 *      from other websites.
 *      
 *      NOTE: this program is created for full symmetric L+U matrices
 *      NOTE: this is the pull version of pagerank
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

double pagerank_upc(sparsemat_t* L, shared  double* page_rank, double* lpage_rank, double* page_rank_curr_iter) {
    // Start timing
    double t1 = wall_seconds();
    upc_barrier;


    double damping = 0.85;
    int64_t num_iterations = 1;
    double error_delta = 0.00001;

    while (1) {
        double error_current = 0.0;
        double error_all_pes = 0.0;

        // reset page rank array for current iteration
        for (int64_t x = 0; x < L->lnumrows; x++) page_rank_curr_iter[x] = 0.0;

        // calculate current iteration page rank (basic calculation)
        for (int64_t y = 0; y < L->lnumrows; y++) {
            for (int64_t k = L->loffset[y]; k < L->loffset[y+1]; k++) {
                int64_t pe_src = L->lnonzero[k] % THREADS;
                int64_t src_idx = L->lnonzero[k] / THREADS;

                double pr_y = page_rank[L->lnonzero[k]];
                int64_t deg_y_upper = L->offset[L->lnonzero[k] + THREADS];
                int64_t deg_y_lower = L->offset[L->lnonzero[k]];
                page_rank_curr_iter[y] += pr_y / (double)(deg_y_upper - deg_y_lower);

            }
        }
        upc_barrier;

        // calculate page rank (full calculation) and current error
        for (int64_t i = 0; i < L->lnumrows; i++) {
            page_rank_curr_iter[i] = ((1-damping) / L->numrows) + (damping * page_rank_curr_iter[i]);
            error_current = error_current + fabs(page_rank_curr_iter[i] - lpage_rank[i]);
        }
        upc_barrier;

        error_all_pes = myupc_reduce_add_d(error_current);
        upc_barrier;
        //printf("error: %f\n", error_all_pes);

        // check for termination (error)
        if (error_all_pes < error_delta) break;

        // update new page rank values
        for (int64_t i = 0; i < L->lnumrows; i++) lpage_rank[i] = page_rank_curr_iter[i];
        upc_barrier;

        // increment num of iterations
        num_iterations++;
    }
    upc_barrier;

    //T0_fprintf(stderr, "Number of iterations: %d\n", num_iterations);

    t1 = wall_seconds() - t1;
    return t1;
}

/*
sparsemat_t* generate_kronecker_graph(
    int64_t* B_spec,
    int64_t B_num,
    int64_t* C_spec,
    int64_t C_num,
    int mode)
{

    T0_fprintf(stderr, "Generating Mode %d Kronecker Product graph (A = B X C) with parameters:  ", mode);
        for(int i = 0; i < B_num; i++) T0_fprintf(stderr, "%ld ", B_spec[i]);
    T0_fprintf(stderr, "X ");
        for(int i = 0; i < C_num; i++) T0_fprintf(stderr, "%ld ", C_spec[i]);   
    T0_fprintf(stderr, "\n");

    sparsemat_t* B = gen_local_mat_from_stars(B_num, B_spec, mode);
    sparsemat_t* C = gen_local_mat_from_stars(C_num, C_spec, mode);   
    if(!B || !C) {
        T0_fprintf(stderr,"ERROR: triangles: error generating input!\n"); upc_global_exit(1);
    }

    T0_fprintf(stderr, "B has %ld rows/cols and %ld nnz\n", B->numrows, B->lnnz);
    T0_fprintf(stderr, "C has %ld rows/cols and %ld nnz\n", C->numrows, C->lnnz);

    sparsemat_t* A = kron_prod_dist(B, C, 1);
  
    return A;
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
    // if (printhelp) usage(); // Skipping print help
    int64_t numrows = l_numrows * THREADS;
    if (erdos_renyi_prob == 0.0) { // use nz_per_row to get erdos_renyi_prob
        erdos_renyi_prob = (2.0 * (nz_per_row - 1)) / numrows;
        if (erdos_renyi_prob > 1.0) erdos_renyi_prob = 1.0;
    } else {                     // use erdos_renyi_prob to get nz_per_row
        nz_per_row = erdos_renyi_prob * numrows;
    }

    T0_fprintf(stderr,"Running page rank on %d threads\n", THREADS);
    if (!read_graph && !gen_kron_graph) {
        T0_fprintf(stderr,"Number of rows per thread   (-N)   %ld\n", l_numrows);
        T0_fprintf(stderr,"Erdos Renyi prob (-e)   %g\n", erdos_renyi_prob);
    }

    T0_fprintf(stderr,"Model mask (M) = %ld (should be 1,2,4,8,16 for agi, exstack, exstack2, conveyors, alternates\n", models_mask);  
    T0_fprintf(stderr,"algorithm (a) = %ld (0 for L & L*U, 1 for L & U*L)\n", alg);

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

    T0_fprintf(stderr,"A2 has %ld rows/cols and %ld nonzeros.\n", A2->numrows, A2->nnz);

    T0_fprintf(stderr, "\nRunning PageRank (upc): \n");

    // create page rank array for final values
    int64_t lpr_size = A2->lnumrows;
    int64_t pr_size = lpr_size * THREADS;
    // create and initialize new page rank array
    shared double* page_rank = upc_all_alloc(pr_size, sizeof(double));
    double* lpage_rank = (double*)(page_rank + MYTHREAD);
    for (int64_t x = 0; x < lpr_size; x++) lpage_rank[x] = 1.0/pr_size; // Initalize page ranks to 1/N
    // create and initialize new page rank array for current iteration values
    double* page_rank_curr_iter = (double*)calloc(lpr_size, sizeof(double));//pr_size, sizeof(double));
    for (int64_t x = 0; x < lpr_size; x++) page_rank_curr_iter[x] = 0.0;
    // create and initialize degrees of each vertex
    //int64_t* degrees_total = upc_all_alloc(pr_size, sizeof(int64_t));
    //for (int64_t x = 0; x < lpr_size; x++) degrees_total[x] = A2->loffset[x+1] - A2->loffset[x];

    upc_barrier;

    double laptime_pagerank = 0.0;

    // Running agi model for pagerank
    laptime_pagerank = pagerank_upc(A2, page_rank, lpage_rank, page_rank_curr_iter);
    upc_barrier;
    //for (int64_t x = 0; x < lpr_size; x++) printf("page rank: %f\n", page_rank[x]);
    T0_fprintf(stderr, "  %8.5lf seconds\n", laptime_pagerank);

    /*      CHECKING FOR ANSWER CORRECTNESS      */
    // sum of all page ranks should equal 1
    double sum_pageranks = 0.0;
    double total_sum_pageranks = 0.0;
    for (int64_t i = 0; i < lpr_size; i++) sum_pageranks += page_rank[i];
    total_sum_pageranks = myupc_reduce_add_d(sum_pageranks);
    //T0_fprintf(stderr, "sum: %f\n", sum_pageranks);
    T0_fprintf(stderr, "Page rank answer is %s.\n", fabs(total_sum_pageranks - 1.0) < 0.001 ? "correct" : "not correct");

    upc_barrier;

    return 0;
}

