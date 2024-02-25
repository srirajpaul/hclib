//Stub for custom BFS implementations
#include "common.h"
#include "aml.h"
#ifdef __cplusplus
extern "C" {
#endif
#include "csr_reference.h"
#ifdef __cplusplus
}
#endif
#include "bitmap_reference.h"
#include <stdint.h>
#include <inttypes.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <limits.h>
#include <assert.h>
#include <stdint.h>

#include "libgetput.h"
#include "selector.h"

int *q1,*q2;
int qc,q2c; //pointer to first free element

//VISITED bitmap parameters
unsigned long *visited;
int64_t visited_size;

int64_t *pred_glob,*column;
int *rowstarts;
oned_csr_graph g;

typedef struct visitmsg {
	//both vertexes are VERTEX_LOCAL components as we know src and dest PEs to reconstruct VERTEX_GLOBAL
	int vloc;
	int vfrom;
} visitmsg;

class BFSSelector: public hclib::Selector<1, visitmsg> {
    int *q2;
    int *q2c;
    int64_t *pred_glob;
    void process(visitmsg m, int sender_rank) {
        if (!TEST_VISITEDLOC(m.vloc)) {
            //printf("visited: m.vloc:%d\n", m.vloc);
            SET_VISITEDLOC(m.vloc);
            q2[*q2c] = m.vloc;
            *q2c = *q2c + 1;
            pred_glob[m.vloc] = VERTEX_TO_GLOBAL(sender_rank, m.vfrom);
        }
    }

public:
BFSSelector(int* _q2, int* _q2c, int64_t* _pred_glob): q2(_q2), q2c(_q2c), pred_glob(_pred_glob) {
        mb[0].process = [this](visitmsg pkt, int sender_rank) { this->process(pkt, sender_rank); };
    }
};

//user should provide this function which would be called once to do kernel 1: graph convert
void make_graph_data_structure(const tuple_graph* const tg) {
	//graph conversion, can be changed by user by replacing oned_csr.{c,h} with new graph format 
	convert_graph_to_oned_csr(tg, &g);
	column=g.column;
    rowstarts=(int*)g.rowstarts;

	visited_size = (g.nlocalverts + ulong_bits - 1) / ulong_bits;
    q1 = (int*)xmalloc(g.nlocalverts*sizeof(int)); //100% of vertexes
	q2 = (int*)xmalloc(g.nlocalverts*sizeof(int));
	for(int i=0;i<g.nlocalverts;i++) q1[i]=0,q2[i]=0; //touch memory
	visited = (long unsigned int*)xmalloc(visited_size*sizeof(unsigned long));
	//user code to allocate other buffers for bfs
}


void run_bfs_selector(int64_t root, int64_t* pred) {
    int64_t nvisited;
    long sum;
    unsigned int i,j,k,lvl=1;
    int rownext = 0;
    int rowlast = 1;

    pred_glob=pred;

    CLEAN_VISITED();

    qc=0; sum=1; q2c=0;

    nvisited=1;
    if(VERTEX_OWNER(root) == rank) {
        pred[VERTEX_LOCAL(root)]=root;
        SET_VISITED(root);
        q1[0]=VERTEX_LOCAL(root);
        qc=1;
    }

    // While there are vertices in current level
    while (sum) {
        BFSSelector *bfss_ptr = new BFSSelector(q2, &q2c, pred_glob);
        hclib::finish([=, &q2c]() {
                bfss_ptr->start();
                //for all vertices in current level send visit AMs to all neighbours
                for(int i=0;i<qc;i++)
                    for(int j=rowstarts[q1[i]];j<rowstarts[q1[i]+1];j++) {
                        //printf("j = %d, dest = %d, m = {%d, %d}\n", j, VERTEX_OWNER(COLUMN(j)), VERTEX_LOCAL(COLUMN(j)), q1[i]);
                        int dest = VERTEX_OWNER(COLUMN(j));
                        visitmsg m = {VERTEX_LOCAL(COLUMN(j)), q1[i]};
                        bfss_ptr->send(0, m, dest);
                    }
                bfss_ptr->done(0);
            });
        lgp_barrier();
        delete bfss_ptr;

        qc=q2c;int *tmp=q1;q1=q2;q2=tmp;
		sum=qc;
		aml_long_allsum(&sum);
        //printf("sum = %d\n", sum);
		nvisited+=sum;

		q2c=0;
    }
    lgp_barrier();
}

//user should provide this function which would be called several times to do kernel 2: breadth first search
//pred[] should be root for root, -1 for unrechable vertices
//prior to calling run_bfs pred is set to -1 by calling clean_pred
void run_bfs(int64_t root, int64_t* pred) {
    /*
    static bool first = true;
    if (first) {
        shmem_init();
        first = false;
    }
    */
    
	pred_glob=pred;
	//user code to do bfs
    const char *deps[] = { "system", "bale_actor" };
    //hclib::launch(deps, 2, [=] {
            hclib::finish([=] {
                    run_bfs_selector(root, pred);
                });
    //    });
    //shmem_finalize();
}

//we need edge count to calculate teps. Validation will check if this count is correct
//user should change this function if another format (not standart CRS) used
void get_edge_count_for_teps(int64_t* edge_visit_count) {
	long i,j;
	long edge_count=0;
	for(i=0;i<g.nlocalverts;i++)
		if(pred_glob[i]!=-1) {
			for(j=g.rowstarts[i];j<g.rowstarts[i+1];j++)
				if(COLUMN(j)<=VERTEX_TO_GLOBAL(my_pe(),i))
					edge_count++;
		}
	aml_long_allsum(&edge_count);
	*edge_visit_count=edge_count;
}

//user provided function to initialize predecessor array to whatevere value user needs
void clean_pred(int64_t* pred) {
	int i;
	for(i=0;i<g.nlocalverts;i++) pred[i]=-1;
}

//user provided function to be called once graph is no longer needed
void free_graph_data_structure(void) {
	free_oned_csr_graph(&g);
	free(visited);
}

//user should change is function if distribution(and counts) of vertices is changed
size_t get_nlocalverts_for_pred(void) {
	return g.nlocalverts;
}
