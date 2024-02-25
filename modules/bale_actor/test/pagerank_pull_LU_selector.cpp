/******************************************************************
//
// 
 *****************************************************************/ 
/*! \file pagerank_pull_LU_selector.cpp
 * \brief Demo application that calculates page rank for vertices in a graph. 
 *      PageRank works by counting the number and quality of links to a page to 
 *      determine a rough estimate of how important the website is. The underlying 
 *      assumption is that more important websites are likely to receive more links 
 *      from other websites.
 *      
 *      NOTE: this program is created for full symmetric L+U matrices
 *      NOTE: this is the pull version of pagerank
 */

#include <math.h>
#include <shmem.h>
extern "C" {
#include "spmat.h"
}
#include "selector.h"

#define THREADS shmem_n_pes()
#define MYTHREAD shmem_my_pe()

enum MailBoxType {REQUEST, RESPONSE};

typedef struct PageRankPkt {
    double pr_y;
    double deg_y;
    int64_t dest_idx;
    int64_t src_idx;
    int64_t pe_dest;
} PageRankPkt;

class PageRankSelector: public hclib::Selector<2, PageRankPkt> {
public:
    PageRankSelector(sparsemat_t* mat, double* pg_rnk, double* pg_rnk_cur_itr, int64_t* deg_ttl) : mat_(mat), pg_rnk_(pg_rnk), pg_rnk_cur_itr_(pg_rnk_cur_itr), deg_ttl_(deg_ttl) {
        mb[REQUEST].process = [this] (PageRankPkt pkt, int sender_rank) { 
            this->req_process(pkt, sender_rank);
        };
        mb[RESPONSE].process = [this] (PageRankPkt pkt, int sender_rank) { 
            this->resp_process(pkt, sender_rank);
        };
    }

private:
    //shared variables
    double* pg_rnk_;
    double* pg_rnk_cur_itr_;
    int64_t* deg_ttl_;
    sparsemat_t* mat_;

    void req_process(PageRankPkt pkg, int sender_rank) {
        PageRankPkt pkg2;
        pkg2.pr_y = pg_rnk_[pkg.src_idx];
        pkg2.deg_y = (double)(mat_->loffset[pkg.src_idx+1] - mat_->loffset[pkg.src_idx]);//(double)deg_ttl_[pkg.src_idx];
        pkg2.dest_idx = pkg.dest_idx;
        send(RESPONSE, pkg2, pkg.pe_dest);
    }

    void resp_process(PageRankPkt pkg, int sender_rank) {
        pg_rnk_cur_itr_[pkg.dest_idx] += pkg.pr_y / (double)pkg.deg_y;
        //printf("page rank: %f, deg: %f\n", pkg.pr_y, pkg.deg_y);
    }
};

double calculate_curr_iter_pagerank(sparsemat_t* L, double* page_rank, double* page_rank_curr_iter, int64_t* degrees_total) {
    PageRankSelector* pagerankSelector = new PageRankSelector(L, page_rank, page_rank_curr_iter, degrees_total);

    hclib::finish([=]() {
        pagerankSelector->start();
        PageRankPkt prPKG;

        // calculate current iteration page rank (basic calculation)
        for (int64_t y = 0; y < L->lnumrows; y++) {
            for (int64_t k = L->loffset[y]; k < L->loffset[y+1]; k++) {
                prPKG.src_idx = L->lnonzero[k] / THREADS;
                prPKG.pe_dest = MYTHREAD;
                prPKG.dest_idx = y;
                pagerankSelector->send(REQUEST, prPKG, L->lnonzero[k] % THREADS);
            }
        }

        pagerankSelector->done(REQUEST);
    });

    return 0;
}

double pagerank_selector(sparsemat_t* L, double* page_rank, double* page_rank_curr_iter, int64_t* degrees_total) {
    // Start timing
    double t1 = wall_seconds();
    lgp_barrier();

    hclib::finish([=]() {
        PageRankPkt prPKG;
        int64_t pe;

        double damping = 0.85;
        int64_t num_iterations = 1;
        double error_delta = 0.00001;

        while (true) {
            double error_current = 0.0;

            // reset page rank array for current iteration
            for (int64_t x = 0; x < L->lnumrows; x++) page_rank_curr_iter[x] = 0.0;

            // calculate current iteration page rank (basic calculation)
            calculate_curr_iter_pagerank(L, page_rank, page_rank_curr_iter, degrees_total);

            // calculate page rank (full calculation) and current error
            for (int64_t i = 0; i < L->lnumrows; i++) {
                page_rank_curr_iter[i] = ((1-damping) / L->numrows) + (damping * page_rank_curr_iter[i]);
                error_current = error_current + fabs(page_rank_curr_iter[i] - page_rank[i]);
            }
            lgp_barrier();

            // reduce error val from all pe's
            double error_all_pes = lgp_reduce_add_d(error_current);
            //printf("error: %f\n", error_all_pes);
            lgp_barrier();

            // check for termination (error)
            if (error_all_pes < error_delta) break;

            // update new page rank values
            for (int64_t i = 0; i < L->lnumrows; i++) page_rank[i] = page_rank_curr_iter[i];
            lgp_barrier();
    
            // increment num of iterations
            num_iterations++;
        }

        //T0_fprintf(stderr, "Number of iterations: %d\n", num_iterations);

    });

    lgp_barrier();

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
        T0_fprintf(stderr,"ERROR: triangles: error generating input!\n"); lgp_global_exit(1);
    }

    T0_fprintf(stderr, "B has %ld rows/cols and %ld nnz\n", B->numrows, B->lnnz);
    T0_fprintf(stderr, "C has %ld rows/cols and %ld nnz\n", C->numrows, C->lnnz);

    sparsemat_t* A = kron_prod_dist(B, C, 1);
  
    return A;
}
*/


int main(int argc, char* argv[]) {

    const char *deps[] = { "system", "bale_actor" };
    hclib::launch(deps, 2, [=] {

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

        sparsemat_t *A, *L, *U;

        int n = numrows;
        L = erdos_renyi_random_graph(numrows, erdos_renyi_prob, UNDIRECTED, NOLOOPS, 12345);
        U = transpose_matrix(L);
        if(!U){T0_printf("ERROR: gen_er_graph_dist: U is NULL!\n"); lgp_global_exit(1);}
        int64_t lnnz = L->lnnz + U->lnnz;
        sparsemat_t * A2 = init_matrix(n, n, lnnz, 0);
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


        lgp_barrier();

        T0_fprintf(stderr, "A2 has %ld rows/cols and %ld nonzeros.\n", A2->numrows, A2->nnz);


        /* Using input of full matrix graph created by gen_erdos_renyi_graph_dist */
        //sparsemat_t* A2 = gen_erdos_renyi_graph_dist(1000, 0.5, 0, 0, 12345);
        //clear_matrix(L); free(L);
        //T0_fprintf(stderr, "A2 has %ld rows/cols and %ld nonzeros.\n", A2->numrows, A2->nnz);
        //lgp_barrier();

        T0_fprintf(stderr, "\nRunning PageRank (selector): \n");

        // create page rank array for final values
        int64_t lpr_size = A2->lnumrows;
        int64_t pr_size = lpr_size * THREADS;
        // create and initialize new page rank array
        double* page_rank = (double*)lgp_all_alloc(pr_size, sizeof(double));
        for (int64_t x = 0; x < lpr_size; x++) page_rank[x] = 1.0/pr_size; // Initalize page ranks to 1/N
        // create and initialize new page rank array for current iteration values
        double* page_rank_curr_iter = (double*)lgp_all_alloc(pr_size, sizeof(double));
        for (int64_t x = 0; x < lpr_size; x++) page_rank_curr_iter[x] = 0.0;
        // create and initialize degrees of each vertex
        int64_t* degrees_total = (int64_t*)lgp_all_alloc(pr_size, sizeof(int64_t));
        for (int64_t x = 0; x < lpr_size; x++) degrees_total[x] = A2->loffset[x+1] - A2->loffset[x];

        lgp_barrier();

        double laptime_pagerank = 0.0;

        // Running selector model for pagerank
        laptime_pagerank = pagerank_selector(A2, page_rank, page_rank_curr_iter, degrees_total);
        lgp_barrier();

        T0_fprintf(stderr, "  %8.5lf seconds\n", laptime_pagerank);

        /*      CHECKING FOR ANSWER CORRECTNESS      */
        // sum of all page ranks should equal 1
        double sum_pageranks = 0.0;
        double total_sum_pageranks = 0.0;
        for (int64_t i = 0; i < lpr_size; i++) sum_pageranks += page_rank[i];
        total_sum_pageranks = lgp_reduce_add_d(sum_pageranks);
        T0_fprintf(stderr, "Page rank answer is %s.\n", abs(total_sum_pageranks - 1.0) < 0.001 ? "correct" : "not correct");

        lgp_barrier();
    
    });

    return 0;
}

