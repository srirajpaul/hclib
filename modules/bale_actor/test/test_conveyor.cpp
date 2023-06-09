
#include <shmem.h>
#include <stdio.h>
extern "C" {
#include <convey.h>
}

int main() {
    shmem_init();
    convey_t* conveyor = convey_new_elastic(ELASTIC_BUFFER_SIZE, SIZE_MAX, 0, NULL, 0);
    //convey_new(SIZE_MAX, 0, NULL, 0);
    if(!conveyor){printf("ERROR: histo_conveyor: convey_new failed!\n"); return(-1.0);}
    int ret = convey_begin(conveyor, sizeof(int64_t), 0);
    if(ret < 0){printf("ERROR: histo_conveyor: begin failed!\n"); return(-1.0);}

    printf("Finished inititialization\n");
    int num = 10, i=0, pe = (shmem_my_pe() + 1)%shmem_n_pes();
    while(convey_advance(conveyor, i == num)) {
        for(; i< num; i++){
            int64_t val = shmem_my_pe() * 1000 + i;
            //if( !convey_push(conveyor, &val, pe)) break;
            if(!convey_epush(conveyor, 8, &val, pe)) break;
        }
        int64_t pop;
        while( convey_epull(conveyor, &pop) == convey_OK) {
            printf("In rank %d, val %d\n", shmem_my_pe(), pop);
        }
    }
    convey_free(conveyor);
    shmem_finalize();
    return 0;
}
