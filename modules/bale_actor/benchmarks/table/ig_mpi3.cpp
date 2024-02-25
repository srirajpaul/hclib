
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

/*! \file histo_shmem.cpp
 * \brief An implementation of histogram that uses OpenSHMEM global atomics.
 */
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <assert.h>
#include <mpi.h>

int THREADS;
int MYTHREAD;

#include <time.h>
#include <sys/time.h>

#define T0_fprintf if(MYTHREAD==0) fprintf

double wall_seconds() {
  struct timeval tp;
  int retVal = gettimeofday(&tp,NULL);
  if (retVal == -1) { perror("gettimeofday:"); fflush(stderr); }
  return ( (double) tp.tv_sec + (double) tp.tv_usec * 1.e-6 );
}

/*!
 * \brief This routine implements straight forward,
 *         single word atomic updates to implement histogram.
 * \param *index array of indices into the shared array of counts.
 * \param l_num_ups the local length of the index array
 * \param *counts SHARED pointer to the count array.
 * \return average run time
 *
 */
double ig_mpi3(int64_t *tgt, int64_t *pckindx, int64_t l_num_req,  int64_t *table, MPI_Win win) {
  int64_t i;
  int64_t pe, col;
  double tm;
  double stat;

  MPI_Win_fence(0, win); // start epoch
  tm = wall_seconds();
  for(i = 0; i < l_num_req; i++){
    col = pckindx[i] >> 16;
    pe  = pckindx[i] & 0xffff;
    MPI_Get(&tgt[i], 1, MPI_INT64_T, pe, col, 1, MPI_INT64_T, win);
  }
  MPI_Win_fence(0, win); // end epoch
  tm = wall_seconds() - tm;
  double sum;
  MPI_Reduce(&tm, &sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  stat = sum / THREADS;
  return stat;
}

int64_t ig_check_and_zero(int64_t use_model, int64_t *tgt, int64_t *index, int64_t l_num_req) {
  int64_t errors=0;
  int64_t i;
  MPI_Barrier(MPI_COMM_WORLD);
  for(i=0; i<l_num_req; i++){
    if(tgt[i] != (-1)*(index[i] + 1)){
      errors++;
      if(errors < 5)  // print first five errors, report all the errors
        fprintf(stderr,"ERROR: model %ld: Thread %d: tgt[%ld] = %ld != %ld)\n",
                use_model,  MYTHREAD, i, tgt[i], (-1)*(index[i] + 1));
               //use_model,  MYTHREAD, i, tgt[i],(-1)*(i*THREADS+MYTHREAD + 1) );
    }
    tgt[i] = 0;
  }
  if( errors > 0 )
    fprintf(stderr,"ERROR: %ld: %ld total errors on thread %d\n", use_model, errors, MYTHREAD);
  MPI_Barrier(MPI_COMM_WORLD);
  return(errors);
}

int main(int argc, char * argv[]) {
  MPI_Init(&argc, &argv);

  // Get the number of processes
  MPI_Comm_size(MPI_COMM_WORLD, &THREADS);

  // Get the rank of the process
  MPI_Comm_rank(MPI_COMM_WORLD, &MYTHREAD);
  if (1) {
    // Get the name of the processor
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    int name_len;
    MPI_Get_processor_name(processor_name, &name_len);

    // Print off a hello world message
    fprintf(stderr,"Hello world from processor %s, rank %d out of %d processors\n",
           processor_name, MYTHREAD, THREADS);
  }

  int64_t i;
  int64_t ltab_siz = 100000;
  int64_t l_num_req  = 1000000;      // number of requests per thread
  int64_t num_errors = 0L, total_errors = 0L;
  int64_t printhelp = 0;

  int opt;
  while( (opt = getopt(argc, argv, "hn:T:")) != -1 ) {
    switch(opt) {
    case 'h': printhelp = 1; break;
    case 'n': sscanf(optarg,"%ld" ,&l_num_req);   break;
    case 'T': sscanf(optarg,"%ld" ,&ltab_siz);   break;
    default:  break;
    }
  }

  T0_fprintf(stderr,"Running ig on %d PEs\n", THREADS);
  T0_fprintf(stderr,"Number of Request / PE           (-n)= %ld\n", l_num_req );
  T0_fprintf(stderr,"Table size / PE                  (-T)= %ld\n", ltab_siz);
  fflush(stderr);

  int tsize;
  MPI_Type_size(MPI_INT64_T,&tsize);
  T0_fprintf(stderr, "sizeof(int64_t) = %d, MPI_Type_size(MPI_INT64_T) = %d\n", sizeof(int64_t), tsize);
  assert(sizeof(int64_t) == tsize);

  // Allocate and zero out the counts array
  int64_t tab_siz = ltab_siz*THREADS;
  int64_t *table;
  MPI_Win win;
  int err = MPI_Win_allocate(ltab_siz* sizeof(int64_t), sizeof(int64_t), MPI_INFO_NULL, MPI_COMM_WORLD, &table, &win);
  if (err != MPI_SUCCESS || !table) {
    T0_fprintf(stderr, "Table allocation failed\n");
  }

  int64_t *ltable = table;
  // fill the table with the negative of its shared index
  // so that checking is easy
  for(i=0; i<ltab_siz; i++)
    ltable[i] = (-1)*(i*THREADS + MYTHREAD + 1);
  // As in the histo example, index is used by the * version.
  // pckindx is used my the buffered versions
  int64_t *index   =  (int64_t*)calloc(l_num_req, sizeof(int64_t)); assert(index != NULL);
  int64_t *pckindx =  (int64_t*)calloc(l_num_req, sizeof(int64_t)); assert(pckindx != NULL);
  int64_t indx, lindx, pe;
  srand(MYTHREAD+ 5 );
  for(i = 0; i < l_num_req; i++){
    indx = rand() % tab_siz;
    index[i] = indx;
    lindx = indx / THREADS;      // the distributed version of indx
    pe  = indx % THREADS;
    pckindx[i] = (lindx << 16) | (pe & 0xffff); // same thing stored as (local index, thread) "shmem style"
  }

  int64_t *tgt  =  (int64_t*)calloc(l_num_req, sizeof(int64_t)); assert(tgt != NULL);

  MPI_Barrier(MPI_COMM_WORLD);

  int64_t use_model;
  double laptime = 0.0;
  int64_t num_models = 0L;               // number of models that are executed
  {
    laptime = ig_mpi3(tgt, pckindx, l_num_req, ltable, win);
    num_models++;
    T0_fprintf(stderr,"  %8.3lf seconds\n", laptime);
    num_errors += ig_check_and_zero(use_model, tgt, index, l_num_req);
  }

  total_errors = num_errors;
  if( total_errors ) {
    T0_fprintf(stderr,"YOU FAILED!!!!\n");
  }

  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Win_free(&win);
  free(index);
  free(pckindx);
  free(tgt);

  MPI_Finalize();

  return 0;
}
