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

#include "shmem.h"
extern "C" {
#include <spmat.h>
}
#include "selector.h"

#define THREADS shmem_n_pes()
#define MYTHREAD shmem_my_pe()

/*! \file transpose_matrix.upc
 * \brief Demo program that runs the variants of transpose_matrix kernel.
 */

/*! 
 * \page transpose_matrix_page Transpose Matrix
 *
 * Demo program that runs the variants of transpose matrix kernel. This program
 * generates a random square asymmetrix (via the Erdos-Renyi model) and then transposes
 * it in parallel.
 *
 * See files spmat_agi.upc, spmat_exstack.upc, spmat_exstack2.upc, and spmat_conveyor.upc
 * for the source for the kernels.
 * 
 * Usage:
 * transpose_matrix [-h][-b count][-M mask][-n num][-T tabsize][-c num]
 * - -h Print usage banner only
 * - -b count is the number of packages in an exstack(2) buffer
 * - -e=p Set the Erdos-Renyi probability to p.
 * - -M=m Set the models mask (1,2,4,8,16,32 for AGI, exstack, exstack2,conveyor,alternate)
 * - -n=n Set the number of rows per PE to n (default = 1000).
 * - -s=s Set a seed for the random number generation.
 * - -Z=z Set the avg number of nonzeros per row to z (default = 10, overrides Erdos-Renyi p).
 */

typedef struct pkg_rowcol_t {
    int64_t row;    
    int64_t col;
} pkg_rowcol_t;

enum MailBoxType {REQUEST};

class TransposeMatrixSelectorPhase1 : public hclib::Selector<1, int64_t> {
public:
    TransposeMatrixSelectorPhase1(int64_t* lcounts, int64_t& lnnz) : lcounts_(lcounts), lnnz_(lnnz) {
        mb[REQUEST].process = [this] (int64_t col, int sender_rank) { 
            this->precessMsg(col, sender_rank);
        };
    }

private:
    int64_t* lcounts_;
    int64_t& lnnz_;

    void precessMsg(int64_t col, int sender_rank) {
        lcounts_[col]++;
        lnnz_++;
    }
};

class TransposeMatrixSelectorPhase2 : public hclib::Selector<1, pkg_rowcol_t> {
public:
    TransposeMatrixSelectorPhase2(uint64_t& numtimespop, sparsemat_t* At, int64_t* wrkoff) : numtimespop_(numtimespop), At_(At), wrkoff_(wrkoff) {
        mb[REQUEST].process = [this] (pkg_rowcol_t pkg, int sender_rank) { 
            this->precessMsg(pkg, sender_rank);
        };
    }

private:
    uint64_t& numtimespop_;
    sparsemat_t* At_;
    int64_t* wrkoff_;

    void precessMsg(pkg_rowcol_t pkg, int sender_rank) {
        numtimespop_++;
        At_->lnonzero[At_->loffset[pkg.col] + wrkoff_[pkg.col]] = pkg.row;
        wrkoff_[pkg.col]++;
    }
};

sparsemat_t* transpose_matrix_selector(sparsemat_t* A) {
    int64_t lnnz = 0;

    /* get the colcnts */
    int64_t lnumcols = (A->numrows + THREADS - MYTHREAD - 1) / THREADS;  
    int64_t* lcounts = (int64_t*)calloc(lnumcols, sizeof(int64_t));
    lgp_barrier();

    TransposeMatrixSelectorPhase1* phase1Selector = new TransposeMatrixSelectorPhase1(lcounts, lnnz);

    hclib::finish([A, phase1Selector]() {
        phase1Selector->start();
        for (int64_t i = 0; i < A->lnnz; i++) {
            int64_t col = A->lnonzero[i] / THREADS;
            int64_t pe = A->lnonzero[i] % THREADS;

            phase1Selector->send(REQUEST ,col, pe);
        }
        phase1Selector->done(REQUEST);
    });
    delete phase1Selector;

    int64_t sum = lgp_reduce_add_l(lnnz);
    assert(A->nnz == sum); 
  
    sparsemat_t* At = init_matrix(A->numcols, A->numrows, lnnz);
    if (!At) {
        printf("ERROR: transpose_matrix: init_matrix failed!\n");
        
        return nullptr;
    }

    /* convert colcounts to offsets */
    At->loffset[0] = 0;  
    for (int64_t i = 1; i < At->lnumrows+1; i++)
        At->loffset[i] = At->loffset[i-1] + lcounts[i-1];

    lgp_barrier();
  
    /* redistribute the nonzeros */
    int64_t* wrkoff = (int64_t*)calloc(At->lnumrows, sizeof(int64_t));
    if (!wrkoff) {
        printf("ERROR: transpose_matrix: wrkoff alloc fail!\n");
        
        return nullptr;
    }

    uint64_t numtimespop = 0;
    TransposeMatrixSelectorPhase2* phase2Selector = new TransposeMatrixSelectorPhase2(numtimespop, At, wrkoff);

    hclib::finish([A, phase2Selector]() {
        phase2Selector->start();
        int64_t row = 0;
        pkg_rowcol_t pkg_nz;

        for (int64_t i = 0; i < A->lnnz; i++) {
            while (i == A->loffset[row + 1]) row++;

            pkg_nz.row = row * THREADS + MYTHREAD;
            pkg_nz.col = A->lnonzero[i] / THREADS;
            int64_t pe = A->lnonzero[i] % THREADS;

            phase2Selector->send(REQUEST, pkg_nz, pe);
        }
        phase2Selector->done(REQUEST);
    });
  
    lgp_barrier();
    delete phase2Selector;

    numtimespop = lgp_reduce_add_l(numtimespop);
    if (numtimespop != A->nnz) {
        printf("ERROR: numtimespop %ld \n", numtimespop);
        printf("%d wrkoff %ld\n", MYTHREAD, wrkoff[0]);

        return nullptr;
    }

    for (int64_t i = 0; i < At->lnumrows; i++) {
        if (wrkoff[i] != lcounts[i]) {
            printf("ERROR: %d wrkoff[%ld] = %ld !=  %ld = lcounts[%ld]\n", MYTHREAD, i, wrkoff[i], lcounts[i], i);
        
            return nullptr;
        }
    }

    free(wrkoff);
    free(lcounts);

    return At;
}

int main(int argc, char** argv) {
    const char *deps[] = { "system", "bale_actor" };
    hclib::launch(deps, 2, [=] {
        int64_t i;
        int64_t check = 1;
        int64_t models_mask = 0xF;
        int printhelp = 0;
        double erdos_renyi_prob = 0.0;
        int64_t nz_per_row = -1;
        int64_t l_numrows = 10000;
        int64_t numrows;
        int64_t seed = 101892+MYTHREAD;
        sparsemat_t * inmat, * outmat;

        int opt; 
        while( (opt = getopt(argc, argv, "hCe:n:M:s:Z:")) != -1 ) {
            switch(opt) {
                case 'h': printhelp = 1; break;
                case 'C': check = 0; break;
                case 'e': sscanf(optarg,"%lf", &erdos_renyi_prob); break;
                case 'n': sscanf(optarg,"%ld", &l_numrows);   break;
                case 'M': sscanf(optarg,"%ld", &models_mask);  break;
                case 's': sscanf(optarg,"%ld", &seed); break;
                case 'Z': sscanf(optarg,"%ld", &nz_per_row);  break;
                default:  break;
            }
        }
    
        numrows = l_numrows * THREADS;

        /* set erdos_renyi_prob and nz_per_row to be consistent */
        if (nz_per_row == -1 && erdos_renyi_prob == 0) {
            nz_per_row = 10;
        } else if (nz_per_row == -1) {
            nz_per_row = erdos_renyi_prob*numrows;
        }
        
        erdos_renyi_prob = (2.0 * (nz_per_row - 1)) / numrows;
        if (erdos_renyi_prob > 1.0)
            erdos_renyi_prob = 1.0;

        T0_fprintf(stderr,"Running transpose_matrix on %d PEs\n", THREADS);
        T0_fprintf(stderr,"Erdos-Renyi edge probability(-e) = %lf\n", erdos_renyi_prob);
        T0_fprintf(stderr,"rows per PE (-n)                 = %ld\n", l_numrows);
        T0_fprintf(stderr,"models_mask (-M)                 = %ld or of 1,2,4,8,16 for gets,classic,exstack2,conveyor,alternate\n", models_mask);
        T0_fprintf(stderr,"seed (-s)                        = %ld\n", seed);
        T0_fprintf(stderr,"Avg # of nonzeros per row   (-Z) = %ld\n", nz_per_row);


        double t1;
        minavgmaxD_t stat[1];
        
        inmat = gen_erdos_renyi_graph_dist(numrows, erdos_renyi_prob, 0, 3, seed + 2);
        if (!inmat) {
            T0_printf("ERROR: inmat is null!\n");
            assert(false);
        }
    
        t1 = wall_seconds();
        outmat = transpose_matrix_selector(inmat);
        T0_fprintf(stderr, "Selector:     \n");
        t1 = wall_seconds() - t1;

        lgp_min_avg_max_d(stat, t1, THREADS);
        T0_fprintf(stderr, " %8.3lf seconds\n", stat->avg);

#if USE_ERROR_CHECK
        /* correctness check */
        if (check) {      
            sparsemat_t * outmatT = transpose_matrix(outmat);
            if (compare_matrix(outmatT, inmat)) {
                T0_fprintf(stderr,"ERROR: transpose of transpose does not match!\n");
            }
            clear_matrix(outmatT);
        }
#endif // USE_ERROR_CHECK

        clear_matrix(outmat);
        clear_matrix(inmat);
        lgp_barrier();
    });

    return 0;
}

