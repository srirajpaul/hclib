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

#include <math.h>
#include <shmem.h>
extern "C" {
#include "spmat.h"
}
#include <std_options.h>
#undef USE_LAMBDA
#include "selector.h"

#define THREADS shmem_n_pes()
#define MYTHREAD shmem_my_pe()

int64_t *count_global = nullptr;
int64_t *get_cnt() { return count_global; }

sparsemat_t *L_global = nullptr;
sparsemat_t *get_mat() { return L_global; }

double triangle_selector(int64_t *count, int64_t *sr, sparsemat_t *L,
                         sparsemat_t *U, int64_t alg) {
  int64_t numpushed = 0;

  if (!L) {
    T0_printf("ERROR: triangle_selector: NULL L!\n");
  }

  // Start timing
  double t1 = wall_seconds();

  if (alg == 0) {
    // TriangleSelector* triSelector = new TriangleSelector(count, L);
    // hclib::Actor<>* triSelector = new hclib::Actor<>();

    struct packet0 {
      int64_t vj;
      int64_t w;
    };
    class SynthesizedActor0 : public hclib::Selector<1, packet0> {
    public:
      void process0(packet0 pkt, int sender_rank) {
        int64_t tempCount = 0;

        for (int64_t k = get_mat()->loffset[pkt.vj];
             k < get_mat()->loffset[pkt.vj + 1]; k++) {
          if (pkt.w == get_mat()->lnonzero[k]) {
            tempCount++;
            break;
          }
          if (pkt.w <
              get_mat()->lnonzero[k]) { // requires that nonzeros are increasing
            break;
          }
        }

        // update the share variable
        *(get_cnt()) += tempCount;
      }
      SynthesizedActor0() {
        mb[0].process = [this](packet0 pkt, int sender_rank) {
          this->process0(pkt, sender_rank);
        };
      }
    };
    SynthesizedActor0 *triSelector = new SynthesizedActor0();
    // hclib::Selector<1> *triSelector = new hclib::Selector<1>();
    count_global = count;
    L_global = L;

    hclib::finish([=, &numpushed]() {
      triSelector->start();
      int64_t k, kk, pe;
      int64_t l_i, L_i, L_j;

      // foreach nonzero (i, j) in L
      for (l_i = 0; l_i < L->lnumrows; l_i++) {
        for (k = L->loffset[l_i]; k < L->loffset[l_i + 1]; k++) {
          L_i = l_i * THREADS + MYTHREAD;
          L_j = L->lnonzero[k];

          pe = L_j % THREADS;
          int64_t vj = L_j / THREADS;
          for (kk = L->loffset[l_i]; kk < L->loffset[l_i + 1]; kk++) {
            int64_t w = L->lnonzero[kk];

            if (w > L_j) {
              break;
            }

            numpushed++;
            packet0 pkt0;
            pkt0.vj = vj;
            pkt0.w = w;
            triSelector->send(0, pkt0, pe);
          }
        }
      }
      // Indicate that we are done with sending messages to the REQUEST mailbox
      triSelector->done(0);
    });
  } else {
    assert(false);
#if 0
        if (!U) {
            T0_printf("ERROR: triangle_selector: NULL U!\n");
            assert(false);
        }

        TriangleSelector* triSelector = new TriangleSelector(count, U);

        hclib::finish([=, &numpushed] () {
            triSelector->start();
            int64_t k,kk, pe;
            int64_t l_i, L_i, L_j;
            TrianglePkt pkg;

            //foreach nonzero (i, j) in L
            for (l_i = 0; l_i < L->lnumrows; l_i++) {
                for (k = L->loffset[l_i]; k < L->loffset[l_i + 1]; k++) {
                    L_i = l_i * THREADS + MYTHREAD;
                    L_j = L->lnonzero[k];

                    pe = L_j % THREADS;
                    pkg.vj = L_j / THREADS;

                    for (kk = U->loffset[l_i]; kk < U->loffset[l_i + 1]; kk++) {
                        pkg.w = U->lnonzero[kk];
                        numpushed++;
                        triSelector->send(REQUEST, pkg, pe);
                    }
                }
            }
            // Indicate that we are done with sending messages to the REQUEST mailbox
            triSelector->done(REQUEST);
        });
#endif
  }

  lgp_barrier();
  *sr = numpushed;
  // Updating count is not necessary since that is taken cared by the mailbox
  // logic *count = cnt;
  minavgmaxD_t stat[1];
  t1 = wall_seconds() - t1;
  lgp_min_avg_max_d(stat, t1, THREADS);

  return t1;
}

sparsemat_t *generate_kronecker_graph(int64_t *B_spec, int64_t B_num,
                                      int64_t *C_spec, int64_t C_num,
                                      int mode) {

  T0_fprintf(stderr,
             "Generating Mode %d Kronecker Product graph (A = B X C) with "
             "parameters:  ",
             mode);
  for (int i = 0; i < B_num; i++)
    T0_fprintf(stderr, "%ld ", B_spec[i]);
  T0_fprintf(stderr, "X ");
  for (int i = 0; i < C_num; i++)
    T0_fprintf(stderr, "%ld ", C_spec[i]);
  T0_fprintf(stderr, "\n");

  sparsemat_t* B = kronecker_product_of_stars(B_num, B_spec, mode);
  sparsemat_t* C = kronecker_product_of_stars(C_num, C_spec, mode);
  if (!B || !C) {
    T0_fprintf(stderr, "ERROR: triangles: error generating input!\n");
    lgp_global_exit(1);
  }

  T0_fprintf(stderr, "B has %ld rows/cols and %ld nnz\n", B->numrows, B->lnnz);
  T0_fprintf(stderr, "C has %ld rows/cols and %ld nnz\n", C->numrows, C->lnnz);

  sparsemat_t* A = kronecker_product_graph_dist(B, C);

  return A;
}

int main(int argc, char *argv[]) {
  const char *deps[] = {"system", "bale_actor"};
  hclib::launch(deps, 2, [=] {
    int64_t buf_cnt = 1024;
    int64_t models_mask = ALL_Models; // default is running all models
    int64_t l_numrows = 10000;        // number of a rows per thread
    int64_t nz_per_row =
        35; // target number of nonzeros per row (only for Erdos-Renyi)
    int64_t read_graph = 0L; // read graph from a file
    char filename[64];
    int64_t cores_per_node = 0;

    double t1;
    int64_t i, j;
    int64_t alg = 0;
    int64_t gen_kron_graph = 0L;
    int kron_graph_mode = 0;
    char *kron_graph_string;
    double erdos_renyi_prob = 0.0;

    int printhelp = 0;
    int opt;
    while ((opt = getopt(argc, argv, "hb:c:M:n:f:a:e:K:")) != -1) {
      switch (opt) {
      case 'h':
        printhelp = 1;
        break;
      case 'b':
        sscanf(optarg, "%ld", &buf_cnt);
        break;
      case 'c':
        sscanf(optarg, "%ld", &cores_per_node);
        break;
      case 'M':
        sscanf(optarg, "%ld", &models_mask);
        break;
      case 'n':
        sscanf(optarg, "%ld", &l_numrows);
        break;
      case 'f':
        read_graph = 1;
        sscanf(optarg, "%s", filename);
        break;

      case 'a':
        sscanf(optarg, "%ld", &alg);
        break;
      case 'e':
        sscanf(optarg, "%lg", &erdos_renyi_prob);
        break;
      case 'K':
        gen_kron_graph = 1;
        kron_graph_string = optarg;
        break;
      default:
        break;
      }
    }

    // if (printhelp) usage(); // Skipping print help
    int64_t numrows = l_numrows * THREADS;
    if (erdos_renyi_prob == 0.0) { // use nz_per_row to get erdos_renyi_prob
      erdos_renyi_prob = (2.0 * (nz_per_row - 1)) / numrows;
      if (erdos_renyi_prob > 1.0)
        erdos_renyi_prob = 1.0;
    } else { // use erdos_renyi_prob to get nz_per_row
      nz_per_row = erdos_renyi_prob * numrows;
    }

    T0_fprintf(stderr, "Running triangle on %d threads\n", THREADS);
    if (!read_graph && !gen_kron_graph) {
      T0_fprintf(stderr, "Number of rows per thread   (-N)   %ld\n", l_numrows);
      T0_fprintf(stderr, "Erdos Renyi prob (-e)   %g\n", erdos_renyi_prob);
    }

    T0_fprintf(stderr,
               "Model mask (M) = %ld (should be 1,2,4,8,16 for agi, exstack, "
               "exstack2, conveyors, alternates\n",
               models_mask);
    T0_fprintf(stderr, "algorithm (a) = %ld (0 for L & L*U, 1 for L & U*L)\n",
               alg);

    double correct_answer = -1;

    sparsemat_t *A, *L, *U;

    if (read_graph) {
      A = read_matrix_mm_to_dist(filename);
      if (!A)
        assert(false);

      T0_fprintf(stderr, "Reading file %s...\n", filename);
      T0_fprintf(stderr, "A has %ld rows/cols and %ld nonzeros.\n", A->numrows,
                 A->nnz);

      // we should check that A is symmetric!

      if (!is_lower_triangular(
              A, 0)) { // if A is not lower triangular... make it so.
        T0_fprintf(
            stderr,
            "Assuming symmetric matrix... using lower-triangular portion...\n");
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
      char *ptr = kron_graph_string;
      int64_t *kron_specs = (int64_t *)calloc(32, sizeof(int64_t *));

      // read the mode
      int ret = sscanf(ptr, "%d ", &kron_graph_mode);
      if (ret == 0)
        ret = sscanf(ptr, "\"%d ", &kron_graph_mode);
      if (ret == 0) {
        T0_fprintf(stderr, "ERROR reading kron graph string!\n");
        assert(false);
      }
      T0_fprintf(stderr, "kron string: %s return = %d\n", ptr, ret);
      T0_fprintf(stderr, "kron mode: %d\n", kron_graph_mode);
      ptr += 2;
      int mat, num_ints = 0;
      while (sscanf(ptr, "%d %n", &num, &mat) == 1) {
        T0_fprintf(stderr, "%s %d\n", ptr, mat);
        kron_specs[num_ints++] = num;
        ptr += mat;
      }

      if (num_ints <= 1) {
        T0_fprintf(stderr,
                   "ERROR: invalid kronecker product string (%s): must contain "
                   "at least three integers\n",
                   kron_graph_string);
        assert(false);
      }

      /* calculate the number of triangles */
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
        correct_answer = (1.0 / 6.0) * pow(4, num_ints) -
                         pow(2.0, (num_ints - 1)) + 1.0 / 3.0;
      }

      correct_answer = round(correct_answer);
      T0_fprintf(stderr, "Pre-calculated answer = %ld\n",
                 (int64_t)correct_answer);

      int64_t half = num_ints / 2;

      L = generate_kronecker_graph(kron_specs, half, &kron_specs[half],
                                   num_ints - half, kron_graph_mode);
    } else {
      L = erdos_renyi_random_graph(numrows, erdos_renyi_prob, UNDIRECTED, NOLOOPS, 12345);
    }

    lgp_barrier();
    if (alg == 1)
      U = transpose_matrix(L);

    lgp_barrier();

    T0_fprintf(stderr, "L has %ld rows/cols and %ld nonzeros.\n", L->numrows,
               L->nnz);

    if (!is_lower_triangular(L, 0)) {
      T0_fprintf(stderr, "ERROR: L is not lower triangular!\n");
      assert(false);
    }

    T0_fprintf(stderr, "Run triangle counting ...\n");
    int64_t tri_cnt;       // partial count of triangles on this thread
    int64_t total_tri_cnt; // the total number of triangles on all threads
    int64_t sh_refs;       // number of shared reference or pushes
    int64_t total_sh_refs;

    int64_t *cc = (int64_t *)lgp_all_alloc(L->numrows, sizeof(int64_t));
    int64_t *l_cc = lgp_local_part(int64_t, cc);
    for (i = 0; i < L->lnumrows; i++)
      l_cc[i] = 0;

    lgp_barrier();

    /* calculate col sums */
    for (i = 0; i < L->lnnz; i++) {
      long lindex = L->lnonzero[i] / THREADS;
      long pe = L->lnonzero[i] % THREADS;
      lgp_fetch_and_inc(&cc[lindex], pe);
    }

    lgp_barrier();

    int64_t rtimesc_calc = 0;
    for (i = 0; i < L->lnumrows; i++) {
      int64_t deg = L->loffset[i + 1] - L->loffset[i];
      rtimesc_calc += deg * l_cc[i];
    }

    /* calculate sum (r_i choose 2) */
    int64_t rchoose2_calc = 0;
    for (i = 0; i < L->lnumrows; i++) {
      int64_t deg = L->loffset[i + 1] - L->loffset[i];
      rchoose2_calc += deg * (deg - 1) / 2;
    }

    /* calculate sum (c_i choose 2) */
    int64_t cchoose2_calc = 0;
    for (i = 0; i < L->lnumrows; i++) {
      int64_t deg = l_cc[i];
      cchoose2_calc += deg * (deg - 1) / 2;
    }

    int64_t pulls_calc = 0;
    int64_t pushes_calc = 0;
    if (alg == 0) {
      pulls_calc = lgp_reduce_add_l(rtimesc_calc);
      pushes_calc = lgp_reduce_add_l(rchoose2_calc);
    } else {
      pushes_calc = lgp_reduce_add_l(rtimesc_calc);
      pulls_calc = lgp_reduce_add_l(cchoose2_calc);
    }

    lgp_all_free(cc);

    T0_fprintf(stderr, "Calculated: Pulls = %ld\n            Pushes = %ld\n\n",
               pulls_calc, pushes_calc);

    int64_t use_model;
    double laptime = 0.0;

    tri_cnt = 0;
    total_tri_cnt = 0;
    sh_refs = 0;
    total_sh_refs = 0;

    // only running selector model
    T0_fprintf(stderr, "Running Selector: \n");
    laptime = triangle_selector(&tri_cnt, &sh_refs, L, U, alg);
    lgp_barrier();

    total_tri_cnt = lgp_reduce_add_l(tri_cnt);
    total_sh_refs = lgp_reduce_add_l(sh_refs);
    T0_fprintf(stderr, "  %8.3lf seconds: %16ld triangles\n", laptime,
               total_tri_cnt);
    // T0_fprintf(stderr, "%16ld shared refs\n", total_sh_refs);
    if ((correct_answer >= 0) && (total_tri_cnt != (int64_t)correct_answer)) {
      T0_fprintf(stderr, "ERROR: Wrong answer!\n");
    }

    if (correct_answer == -1) {
      correct_answer = total_tri_cnt;
    }

    lgp_barrier();
  });

  return 0;
}
