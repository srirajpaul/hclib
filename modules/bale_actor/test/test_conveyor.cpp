
#include <shmem.h>
#include <stdio.h>
extern "C" {
#include <convey.h>
}

int main() {
    shmem_init();
    convey_t* conveyor = convey_new(SIZE_MAX, 0, NULL, 0);
    if(!conveyor){printf("ERROR: histo_conveyor: convey_new failed!\n"); return(-1.0);}
    int ret = convey_begin(conveyor, sizeof(int64_t), 0);
    if(ret < 0){printf("ERROR: histo_conveyor: begin failed!\n"); return(-1.0);}

    printf("Finished inititialization\n");
    int num = 10, i=0, pe = (shmem_my_pe() + 1)%shmem_n_pes();
    while(convey_advance(conveyor, i == num)) {
        for(; i< num; i++){
            int64_t val = shmem_my_pe() * 1000 + i;
            if( !convey_push(conveyor, &val, pe)) break;
        }
        int64_t pop;
        while( convey_pull(conveyor, &pop, NULL) == convey_OK) {
            printf("In rank %d, val %d\n", shmem_my_pe(), pop);
        }
    }
    convey_free(conveyor);
    shmem_finalize();
    return 0;
}
