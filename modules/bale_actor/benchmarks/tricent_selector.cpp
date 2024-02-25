/******************************************************************
//
// 
 *****************************************************************/ 
/*! \file tricent_selector.cpp
 * \brief Demo application that calculates triangle centrality for vertices in a graph. 
 *      Triangle centrality involves finding important vertices in graphs according to the 
 *      concentration of triangles in the neighborhood of a vertex.
 *      
 *      NOTE: this program is created for full symmetric L+U matrices
 */

#include <math.h>
#include <shmem.h>
extern "C" {
#include "spmat.h"
}
#include "selector.h"

enum MailBoxType {REQUEST, RESPONSE};

typedef struct TrianglePkt {
    int64_t u;
    int64_t v;
    int64_t w;
    int64_t pos;
    int64_t first_neigh;
    int64_t second_neigh;
} TrianglePkt;

typedef struct TricentPkt {
    int64_t index_u;
    int64_t index_src;
    int64_t pe_src;
    double tri_count;
} TricentPkt;

class TriangleSelector: public hclib::Selector<2, TrianglePkt> {
public:
    TriangleSelector(int64_t* cnt, sparsemat_t* mat, int64_t* cnt_uvw) : cnt_(cnt), mat_(mat), cnt_uvw_(cnt_uvw) {
        mb[REQUEST].process = [this] (TrianglePkt pkt, int sender_rank) { 
            this->req_process(pkt, sender_rank);
        };
        mb[RESPONSE].process = [this] (TrianglePkt pkt, int sender_rank) { 
            this->resp_process(pkt, sender_rank);
        };
    }

private:
    //shared variables
    int64_t* cnt_;
    sparsemat_t* mat_;
    int64_t* cnt_uvw_;

    void req_process(TrianglePkt pkg, int sender_rank) {
        int64_t tempCount = 0;
        
        for (int64_t k = mat_->loffset[pkg.v/THREADS]; k < mat_->loffset[pkg.v/THREADS + 1]; k++) {
            int64_t nnz_x = mat_->lnonzero[k]; if (nnz_x >= mat_->numrows) nnz_x = nnz_x - mat_->numrows; // conversion due to negating triangle neighbor nnz's

            if (pkg.w == nnz_x) {
                tempCount++;
                TrianglePkt triPKG;
                // update T(u) and tri_N(u)
                triPKG.pos = pkg.u;
                triPKG.first_neigh = pkg.w;
                triPKG.second_neigh = pkg.v;
                send(RESPONSE, triPKG, sender_rank);

                // update T(v) and tri_N(v)
                triPKG.pos = pkg.v/THREADS;
                triPKG.first_neigh = pkg.u * THREADS + sender_rank;
                triPKG.second_neigh = pkg.w;
                send(RESPONSE, triPKG, MYTHREAD);

                // update T(w) and tri_N(w)
                triPKG.pos = pkg.w / THREADS;
                triPKG.first_neigh = pkg.u * THREADS + sender_rank;
                triPKG.second_neigh = pkg.v;
                send(RESPONSE, triPKG, pkg.w%THREADS);
                
                break;
            }
            if (pkg.w < nnz_x) { // requires that nonzeros are increasing
                break;
            }
        }

        // update the share variable
        *cnt_ += tempCount;
    }

    void resp_process(TrianglePkt pkg, int sender_rank) {
        cnt_uvw_[pkg.pos]++; // update T(vertex)

        // update tri_N(vertex) for w,u , w,v  ,  u,w , u,v  ,  v,w , v,u
        for (int64_t k = mat_->loffset[pkg.pos]; k < mat_->loffset[pkg.pos + 1]; k++) {
            if (pkg.first_neigh == mat_->lnonzero[k] || pkg.second_neigh == mat_->lnonzero[k]) {
                mat_->lnonzero[k] = mat_->lnonzero[k] + mat_->numrows;
            }
        }
    }
};

class TricentSelector: public hclib::Selector<2, TricentPkt> {
public:
    TricentSelector(double* tricnts, int64_t totaltricnt, int64_t* vrtx_tri_cnt) : tricnts_(tricnts), totaltricnt_(totaltricnt), vrtx_tri_cnt_(vrtx_tri_cnt) {
        mb[REQUEST].process = [this] (TricentPkt pkt, int sender_rank) { 
            this->req_process(pkt, sender_rank);
        };
        mb[RESPONSE].process = [this] (TricentPkt pkt, int sender_rank) { 
            this->resp_process(pkt, sender_rank);
        };
    }  

private:
    //shared variables
    double* tricnts_;
    int64_t totaltricnt_;
    int64_t* vrtx_tri_cnt_;
    
    void req_process(TricentPkt pkg, int sender_rank) {
        TricentPkt pkg2;
        pkg2.index_src = pkg.index_src;
        pkg2.tri_count = (pkg.tri_count * (double)vrtx_tri_cnt_[pkg.index_u]) / (double)totaltricnt_;
        send(RESPONSE, pkg2, pkg.pe_src);
    }

    void resp_process(TricentPkt pkg, int sender_rank) {
        tricnts_[pkg.index_src] += pkg.tri_count;
    }

};

double tricount_selector(int64_t* count, sparsemat_t* L, int64_t* vertex_tri_count) {
    // Start timing
    double t1 = wall_seconds();

    TriangleSelector* tricountSelector = new TriangleSelector(count, L, vertex_tri_count);

    hclib::finish([=]() {
        tricountSelector->start();
        int64_t k,kk, pe;
        int64_t l_i, L_i, L_j;
        TrianglePkt pkg;

        //foreach nonzero (i, j) in L
        for (l_i = 0; l_i < L->lnumrows; l_i++) {
            for (k = L->loffset[l_i]; k < L->loffset[l_i + 1]; k++) {
                L_i = l_i * THREADS + MYTHREAD;

                L_j = L->lnonzero[k]; if (L_j >= L->numrows) L_j = L_j - L->numrows; // conversion due to negating triangle neighbor nnz's

                if (L_i < L_j) continue; // check due to full matrix

                pkg.u = l_i;
                pkg.v = L_j;
                pe = L_j % THREADS;

                for (kk = L->loffset[l_i]; kk < L->loffset[l_i + 1]; kk++) {
                    pkg.w = L->lnonzero[kk]; if (pkg.w >= L->numrows) pkg.w = pkg.w - L->numrows; // conversion due to negating triangle neighbor nnz's

                    if (pkg.w > L_j) {
                        break;
                    }

                    tricountSelector->send(REQUEST, pkg, pe);
                }
            }
        }
        // Indicate that we are done with sending messages to the REQUEST mailbox
        tricountSelector->done(REQUEST);
    });

    lgp_barrier();

    t1 = wall_seconds() - t1;

    return t1;
}

double tricent_selector(int64_t total_tri_cnt, double* tricent_counts, sparsemat_t* L, int64_t* vertex_tri_count) {
    // Start timing
    double t2 = wall_seconds();
    
    TricentSelector* tricentSelector = new TricentSelector(tricent_counts, total_tri_cnt, vertex_tri_count);

    hclib::finish([=]() {
        tricentSelector->start();
        int64_t v, k, u;
        TricentPkt pkg;
        
        for (v = 0; v < L->lnumrows; v++) {
            pkg.index_src = v;
            pkg.pe_src = MYTHREAD;
            for (k = L->loffset[v]; k < L->loffset[v+1]; k++) {
                u = L->lnonzero[k];
                
                if (u >= L->numrows) {
                    // nnz is in tri_N(v) set, update (1/3)*T(u) / T(G)
                    u = u - L->numrows;
                    pkg.index_u = u / THREADS;
                    pkg.tri_count = (double)1/3;
                    tricentSelector->send(REQUEST, pkg, u%THREADS);
                } else {
                    // nnz is in N(v)\tri_N(v), update (1)*T(u) / T(G)
                    pkg.index_u = u / THREADS;
                    pkg.tri_count = (double)1;
                    tricentSelector->send(REQUEST, pkg, u%THREADS);
                }
            }
            // add T(v) from tri_N^+(v)
            tricent_counts[v] += (((double)1/3) * (double)vertex_tri_count[v]) / (double)total_tri_cnt;
        }
        
        tricentSelector->done(REQUEST);
    });
    lgp_barrier();
    t2 = wall_seconds() - t2;
    return t2;
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

        T0_fprintf(stderr,"Running tricent on %d threads\n", THREADS);
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

        T0_fprintf(stderr, "Run triangle centrality ...\n");
        int64_t tri_cnt = 0;           // partial count of triangles on this thread
        int64_t total_tri_cnt = 0;     // the total number of triangles on all threads
        double laptime_tricount = 0.0;
        double laptime_tricent = 0.0;
        int64_t ltc_size = A2->lnumrows;
        int64_t tc_size = ltc_size * THREADS;
        int64_t* vertex_tri_count = (int64_t*)lgp_all_alloc(tc_size, sizeof(int64_t));
        for (int64_t x = 0; x < ltc_size; x++) vertex_tri_count[x] = 0;
        double* tricent_counts = (double*)lgp_all_alloc(tc_size, sizeof(double));
        for (int64_t x = 0; x < ltc_size; x++) tricent_counts[x] = 0.0;

        // Running selector model for triangle counting
        T0_fprintf(stderr, "\nrunning tricount (selector) ... \n");
        laptime_tricount = tricount_selector(&tri_cnt, A2, vertex_tri_count);
        lgp_barrier();
        total_tri_cnt = lgp_reduce_add_l(tri_cnt);
        lgp_barrier();
        //T0_fprintf(stderr, "  %8.5lf seconds: %16ld triangles\n", laptime_tricount, total_tri_cnt);
        
        // Running selector model for triangle centrality
        T0_fprintf(stderr, "\nrunning tricent (selector) ... \n");
        laptime_tricent = tricent_selector(total_tri_cnt, tricent_counts, A2, vertex_tri_count);
        lgp_barrier();
        T0_fprintf(stderr, "  %8.5lf seconds\n", laptime_tricount+laptime_tricent);
        //for (int64_t x = 0; x < ltc_size; x++) printf("tricent: %f\n",tricent_counts[x]);
        lgp_barrier();
        
    });

    return 0;
}

