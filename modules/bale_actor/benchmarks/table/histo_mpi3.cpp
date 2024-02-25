
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
double histo_mpi3(int64_t *pckindx, int64_t T, int64_t *counts, MPI_Win win) {
  int ret;
  int64_t i;
  double tm;
  int64_t pe, col;
  double stat;
  int64_t one = 1;
  int64_t result;

  MPI_Win_fence(0, win); // start epoch
  tm = wall_seconds();
  i = 0UL;
  for(i = 0; i < T; i++) {
    col = pckindx[i] >> 16;
    pe  = pckindx[i] & 0xffff;
#ifdef USE_FETCH_AND_OP
    MPI_Fetch_and_op(&one, &result, MPI_INT64_T, pe, col, MPI_SUM, win);
#elif USE_ACCUMULATE
    MPI_Accumulate(&one, 1, MPI_INT64_T, pe, col, 1, MPI_INT64_T, MPI_SUM, win);
#elif USE_GET_ACCUMULATE
    MPI_Get_accumulate(&one, 1, MPI_INT64_T, &result, 1, MPI_INT64_T, pe, col, 1, MPI_INT64_T, MPI_SUM, win);
#else
    MPI_Abort(MPI_COMM_WORLD, 99);
#endif
  }
  MPI_Win_fence(0, win); // end epoch
  tm = wall_seconds() - tm;
  double sum;
  MPI_Reduce(&tm, &sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  stat = sum / THREADS;
  return stat;
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

  int64_t l_num_ups  = 1000000;     // per thread number of requests (updates)
  int64_t lnum_counts = 1000;       // per thread size of the table

  int64_t i;

  int printhelp = 0;
  int opt;
  while( (opt = getopt(argc, argv, "hn:T:")) != -1 ) {
    switch(opt) {
    case 'h': printhelp = 1; break;
    case 'n': sscanf(optarg,"%ld" ,&l_num_ups);  break;
    case 'T': sscanf(optarg,"%ld" ,&lnum_counts);  break;
    default:  break;
    }
  }

  T0_fprintf(stderr,"Running histo on %d PEs\n", THREADS);
  T0_fprintf(stderr,"Number updates / PE              (-n)= %ld\n", l_num_ups);
  T0_fprintf(stderr,"Table size / PE                  (-T)= %ld\n", lnum_counts);
  fflush(stderr);

  int tsize;
  MPI_Type_size(MPI_INT64_T,&tsize);
  T0_fprintf(stderr, "sizeof(int64_t) = %d, MPI_Type_size(MPI_INT64_T) = %d\n", sizeof(int64_t), tsize);
  assert(sizeof(int64_t) == tsize);

  // Allocate and zero out the counts array
  int64_t num_counts = lnum_counts*THREADS;
  int64_t *counts;
  MPI_Win win;
  int err = MPI_Win_allocate(lnum_counts* sizeof(int64_t), sizeof(int64_t), MPI_INFO_NULL, MPI_COMM_WORLD, &counts, &win);
  if (err != MPI_SUCCESS || !counts) {
    T0_fprintf(stderr, "Table allocation failed\n");
  }

  int64_t *lcounts = counts;
  for(i = 0; i < lnum_counts; i++)
    lcounts[i] = 0L;

  // index is a local array of indices into the shared counts array.
  // This is used by the * version.
  // To avoid paying the UPC tax of computing index[i]/THREADS and index[i]%THREADS
  // when using the exstack and conveyor models
  // we also store a packed version that holds the pe (= index%THREADS) and lindx (=index/THREADS)
  int64_t *index   = (int64_t *) calloc(l_num_ups, sizeof(int64_t)); assert(index != NULL);
  int64_t *pckindx = (int64_t *) calloc(l_num_ups, sizeof(int64_t)); assert(pckindx != NULL);
  int64_t indx, lindx, pe;

  srand(MYTHREAD + 120348);
  for(i = 0; i < l_num_ups; i++) {
    //indx = i % num_counts;          //might want to do this for debugging
    indx = rand() % num_counts;
    index[i] = indx;
    lindx = indx / THREADS;
    pe  = indx % THREADS;
    pckindx[i]  =  (lindx << 16L) | (pe & 0xffff);
  }

  MPI_Barrier(MPI_COMM_WORLD);

  int64_t use_model;
  double laptime = 0.0;
  int64_t num_models = 0L;               // number of models that are executed
  {
    laptime = histo_mpi3(pckindx, l_num_ups, counts, win);
    num_models++;
    T0_fprintf(stderr,"  %8.3lf seconds\n", laptime);
  }
  MPI_Barrier(MPI_COMM_WORLD);

  // Check the results
  // Assume that the atomic add version will correctly zero out the counts array
  int64_t minus_num_models = - num_models;
  int64_t result;
  MPI_Win_fence(0, win); // start epoch
  for(i = 0; i < l_num_ups; i++) {
    lindx = index[i] / THREADS;
    pe  = index[i] % THREADS;
    MPI_Get_accumulate(&minus_num_models, 1, MPI_INT64_T, &result, 1, MPI_INT64_T, pe, lindx, 1, MPI_INT64_T, MPI_SUM, win);
  }
  MPI_Win_fence(0, win); // end epoch
  MPI_Barrier(MPI_COMM_WORLD);

  int64_t num_errors = 0, totalerrors = 0;
  for(i = 0; i < lnum_counts; i++) {
    if(lcounts[i] != 0L) {
      num_errors++;
      if(num_errors < 5)  // print first five errors, report number of errors below
        fprintf(stderr,"ERROR: Thread %d error at %ld (= %ld)\n", MYTHREAD, i, lcounts[i]);
    }
  }
  MPI_Reduce(&num_errors, &totalerrors, 1, MPI_INT64_T, MPI_SUM, 0, MPI_COMM_WORLD);
  if(totalerrors) {
     T0_fprintf(stderr,"FAILED!!!! total errors = %ld\n", totalerrors);
  }
  MPI_Win_free(&win);
  free(index);
  free(pckindx);

  MPI_Finalize();

  return 0;
}
