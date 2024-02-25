/******************************************************************
//
// 
 *****************************************************************/ 
/*! \file jaccard_selector.cpp
 * \brief Demo application that calculates jaccard index for edges in a graph. 
 *      The Jaccard index gives a value that represents the similarity 
 *      of the neighborhoods of two vertices connected by an edge.
 *
 *      NOTE: this program is created for full symmetric L+U matrices
 */

#include <math.h>
#include <shmem.h>
extern "C" {
#include "spmat.h"
}
#include "selector.h"

enum MailBoxType {RESPONSE, REQUEST};

typedef struct JaccardPkt {
    int64_t x;
    int64_t pos;
    int64_t index_u;
} JaccardPkt;

class JaccardSelector: public hclib::Selector<2, JaccardPkt> {
public:
    JaccardSelector(sparsemat_t* mat, double* jac_idx) : mat_(mat), jac_idx_(jac_idx) {
        mb[REQUEST].process = [this] (JaccardPkt pkt, int sender_rank) { 
            this->req_process(pkt, sender_rank);
        };
        mb[RESPONSE].process = [this] (JaccardPkt pkt, int sender_rank) { 
            this->resp_process(pkt, sender_rank);
        };
    }

private:
    //shared variables
    double* jac_idx_;
    sparsemat_t* mat_;

    void req_process(JaccardPkt pkg, int sender_rank) {
        JaccardPkt pkg2;
        for (int64_t uk = mat_->loffset[pkg.index_u]; uk < mat_->loffset[pkg.index_u+1]; uk++) {
            if (pkg.x == mat_->lnonzero[uk]) {
                pkg2.pos = pkg.pos;
                send(RESPONSE, pkg2, sender_rank);
            }
        }
    }

    void resp_process(JaccardPkt pkg, int sender_rank) {
        jac_idx_[pkg.pos]++;
    }

};

typedef struct DegreePkt {
    int64_t i;
    int64_t src_idx;
} DegreePkt;

class DegreeSelector: public hclib::Selector<2, DegreePkt> {
public:
    DegreeSelector(sparsemat_t* mat, int64_t* deg_ttl) : mat_(mat), deg_ttl_(deg_ttl) {
        mb[REQUEST].process = [this] (DegreePkt pkt, int sender_rank) { 
            this->req_process(pkt, sender_rank);
        };
        mb[RESPONSE].process = [this] (DegreePkt pkt, int sender_rank) { 
            this->resp_process(pkt, sender_rank);
        };
    }

private:
    //shared variables
    int64_t* deg_ttl_;
    sparsemat_t* mat_;

    void req_process(DegreePkt pkg, int sender_rank) {
        DegreePkt pkg2;
        pkg2.src_idx = pkg.src_idx;
        pkg2.i = mat_->loffset[pkg.i + 1] - mat_->loffset[pkg.i];
        send(RESPONSE, pkg2, sender_rank);
    }

    void resp_process(DegreePkt pkg, int sender_rank) {
        deg_ttl_[pkg.src_idx] = pkg.i;
    }

};

double get_edge_degrees(sparsemat_t* L, int64_t* edge_degrees) {
    DegreeSelector* degreeSelector = new DegreeSelector(L, edge_degrees);

    hclib::finish([=]() {
        degreeSelector->start();
        DegreePkt degPKG;
    
        int64_t pos = 0;
        for (int64_t y = 0; y < L->lnumrows; y++) {
            for (int64_t k = L->loffset[y]; k < L->loffset[y+1]; k++) {
                if ((y*THREADS + MYTHREAD) < L->lnonzero[k]) continue;
                degPKG.i = L->lnonzero[k] / THREADS;
                degPKG.src_idx = pos;
                degreeSelector->send(REQUEST, degPKG, L->lnonzero[k] % THREADS);
                pos++;
            }
        }

        degreeSelector->done(REQUEST);
    });

    return 0;
}

double jaccard_selector(sparsemat_t* A2, double* jaccard_index, int64_t* edge_degrees) {
    // Start timing
    double t1 = wall_seconds();
    lgp_barrier();

    JaccardSelector* jacSelector = new JaccardSelector(A2, jaccard_index);

    hclib::finish([=]() {
        jacSelector->start();

        int64_t pos = 0;
        JaccardPkt pkg;

        for (int64_t v = 0; v < A2->lnumrows; v++) {             //vertex v (local)
            for (int64_t k = A2->loffset[v]; k < A2->loffset[v+1]; k++) {
                int64_t u = A2->lnonzero[k];                     //vertex u (possibly remote)
                if ((v * THREADS + MYTHREAD) < u) continue;

                // calculate intersection
                for (int64_t kk = A2->loffset[v]; kk < A2->loffset[v+1]; kk++) {
                    int64_t x = A2->lnonzero[kk];
                    if (x == u) continue;
                    pkg.index_u = u / THREADS;
                    pkg.x = x;
                    pkg.pos = pos;
                    jacSelector->send(REQUEST, pkg, u % THREADS);
                }
                
                //increase position in jaccard index local array
                pos++;
            }
        }
        jacSelector->done(REQUEST);
        
    });

    lgp_barrier();


    int64_t pos = 0;
    for (int64_t v = 0; v < A2->lnumrows; v++) {
        int64_t deg_v = A2->loffset[v+1] - A2->loffset[v];
        for (int64_t k = A2->loffset[v]; k < A2->loffset[v+1]; k++) {
            if ((v*THREADS + MYTHREAD) < A2->lnonzero[k]) continue;

            // calculate jaccard index using: intersection / (deg(u) + deg(v) - intersection)
            //printf("edge deg: %d\n", edge_degrees[pos]);
            jaccard_index[pos] = jaccard_index[pos] / ((double)edge_degrees[pos] + (double)deg_v - jaccard_index[pos]);
            
            //increase position in jaccard index local array
            pos++;
        }
    }

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

        T0_fprintf(stderr,"Running jaccard on %d threads\n", THREADS);
        if (!read_graph && !gen_kron_graph) {
            T0_fprintf(stderr,"Number of rows per thread   (-N)   %ld\n", l_numrows);
            T0_fprintf(stderr,"Erdos Renyi prob (-e)   %g\n", erdos_renyi_prob);
        }

        T0_fprintf(stderr,"Model mask (M) = %ld (should be 1,2,4,8,16 for agi, exstack, exstack2, conveyors, alternates\n", models_mask);  
        T0_fprintf(stderr,"algorithm (a) = %ld (0 for L & L*U, 1 for L & U*L)\n", alg);

        double correct_answer = -1;

        sparsemat_t *A, *L, *U;
/*
        if (read_graph) {
            A = read_matrix_mm_to_dist(filename);
            if (!A) assert(false);
            
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

        } else if (gen_kron_graph) {
            // string should be <mode> # # ... #
            // we will break the string of numbers (#s) into two groups and create
            // two local kronecker graphs out of them.
            int num;
            char* ptr = kron_graph_string;
            int64_t* kron_specs = (int64_t*)calloc(32, sizeof(int64_t *));

            // read the mode
            int ret = sscanf(ptr, "%d ", &kron_graph_mode);
            if (ret == 0) ret = sscanf(ptr, "\"%d ", &kron_graph_mode);
            if (ret == 0) { T0_fprintf(stderr, "ERROR reading kron graph string!\n"); assert(false); }
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
                T0_fprintf(stderr, "ERROR: invalid kronecker product string (%s): must contain at least three integers\n", kron_graph_string); 
                assert(false);
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
                for (i = 0; i < num_ints; i++) {
                    x *= (kron_specs[i] + 1);
                }

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
        //}

        lgp_barrier();

        T0_fprintf(stderr, "A2 has %ld rows/cols and %ld nonzeros.\n", A2->numrows, A2->nnz);

        /* Using input of full matrix graph created by gen_erdos_renyi_graph_dist */
        //sparsemat_t* A2 = gen_erdos_renyi_graph_dist(500, 0.5, 0, 0, 12345);
        //sparsemat_t* L = copy_matrix(A2);
        //tril(L, -1);
        //lgp_barrier();

        T0_fprintf(stderr, "\nRunning Jaccard Index (selector): \n");

        int64_t ljaccard_size = lgp_reduce_max_l(A2->lnnz);
        double* jaccard_index = (double*)lgp_all_alloc(ljaccard_size*THREADS, sizeof(double));
        for (int64_t x = 0; x < ljaccard_size; x++) jaccard_index[x] = 0;
        // get the degrees of each vertex
        int64_t* edge_degrees = (int64_t*)lgp_all_alloc(ljaccard_size*THREADS, sizeof(int64_t));
        for (int64_t x = 0; x < ljaccard_size; x++) edge_degrees[x] = 0;
        lgp_barrier();
        get_edge_degrees(A2, edge_degrees);
        
        lgp_barrier();
        //for (int64_t i = 0; i < ljaccard_size; i++) printf("\npe: %d, deg: %d", MYTHREAD, edge_degrees[i]);
        double laptime_jaccard = 0.0;

        // Running selector model for jaccard index
        laptime_jaccard = jaccard_selector(A2, jaccard_index, edge_degrees);
        lgp_barrier();

        //for (int64_t i = 0; i < ljaccard_size; i++) printf("\npe: %d, jaccard index: %f", MYTHREAD, jaccard_index[i]);
        //printf("\n");

        T0_fprintf(stderr, " %8.5lf seconds\n", laptime_jaccard);

        lgp_barrier();

        lgp_finalize();

    });

    return 0;
}
