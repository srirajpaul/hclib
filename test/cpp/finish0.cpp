/* Copyright (c) 2013, Rice University

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

1.  Redistributions of source code must retain the above copyright
     notice, this list of conditions and the following disclaimer.
2.  Redistributions in binary form must reproduce the above
     copyright notice, this list of conditions and the following
     disclaimer in the documentation and/or other materials provided
     with the distribution.
3.  Neither the name of Rice University
     nor the names of its contributors may be used to endorse or
     promote products derived from this software without specific
     prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

 */

/**
 * DESC: top-level finish for a bunch of asyncs
 */
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

#include "hclib_cpp.h"

#define NB_ASYNC 127

int ran[NB_ASYNC];

void init_ran(int *ran, int size) {
    while (size >= 0) {
        ran[size] = -1;
        size--;
    }
}

void assert_done(int start, int end) {
    while(start < end) {
        assert(ran[start] == start);
        start++;
    }
}
int main (int argc, char ** argv) {
    printf("Call Init\n");
    int mid = NB_ASYNC/2;
    hclib::launch([&mid]() {
        int i = 0;
        int indices [NB_ASYNC];

        hclib::finish([=, &i, &indices]() {

            // This is ok to have these on stack because this
            // code is alive until the end of the program.
            init_ran(ran, NB_ASYNC);
            printf("Go over [%d:%d]\n", i, mid);
            while(i < mid) {
                indices[i] = i;
                hclib::async([=]() { int index = indices[i];
                    assert(ran[index] == -1); ran[index] = index; });
                i++;
            }
        });

        printf("Midway\n");
        assert_done(0, mid);
        printf("Go over [%d:%d]\n", i, NB_ASYNC);
        while(i < NB_ASYNC) {
            indices[i] = i;
            //Note: Forcefully pass the address we want to write to as a void **
            hclib::async([=]() { int index = indices[i];
                assert(ran[index] == -1); ran[index] = index; });
            i++;
        }
    });
    printf("Check results: ");
    assert_done(mid, NB_ASYNC);
    printf("OK\n");
    return 0;
}
