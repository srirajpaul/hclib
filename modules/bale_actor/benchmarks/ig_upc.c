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

/*! \file ig_upc.cpp
 * \brief An implementation of indexgather that uses single word gets to shared addresses.
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

/*!
 * \brief This routine implements the single word get version indexgather
 * \param *tgt array of target locations for the gathered values
 * \param *index array of indices into the global array of counts
 * \param l_num_req the length of the index array
 * \param *table shared pointer to the shared table array.
 * \return average run time
 *
 */
double ig_upc(int64_t *tgt, int64_t *index, int64_t l_num_req, shared int64_t *table) {
  int64_t i;
  int64_t pe, col;
  double tm;
  double stat;

  upc_barrier;
  tm = wall_seconds();

  for(i = 0; i < l_num_req; i++){
    #pragma pgas defer_sync
    tgt[i] = table[index[i]];
  }

  tm = wall_seconds() - tm;
  upc_barrier;

  stat = myupc_reduce_add_d(tm)/THREADS;

  return stat;
}

int64_t ig_check_and_zero(int64_t use_model, int64_t *tgt, int64_t *index, int64_t l_num_req) {
  int64_t errors=0;
  int64_t i;
  upc_barrier;
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
  upc_barrier;
  return(errors);
}

int main(int argc, char * argv[]) {
  //char hostname[1024];
  //hostname[1023] = '\0';
  //gethostname(hostname, 1023);
  //fprintf(stderr,"Hostname: %s rank: %d\n", hostname, MYTHREAD);
  
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
  
  // Allocate and populate the shared table array 
  int64_t tab_siz = ltab_siz*THREADS;
  shared int64_t * table  = upc_all_alloc(tab_siz, sizeof(int64_t)); assert(table != NULL);
  int64_t *ltable  = (int64_t*)(table+MYTHREAD);
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

  upc_barrier;

  int64_t use_model;
  double laptime = 0.0;

        laptime = ig_upc(tgt, index, l_num_req, table);

     T0_fprintf(stderr,"  %8.3lf seconds\n", laptime);
   
     num_errors += ig_check_and_zero(use_model, tgt, index, l_num_req);

  total_errors = num_errors;
  if( total_errors ) {
    T0_fprintf(stderr,"YOU FAILED!!!!\n");
  } 

  upc_barrier;
  upc_all_free(table);
  free(index);
  free(pckindx);
  free(tgt);
  return 0;
}

