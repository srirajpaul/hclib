CFLAGS = -Drestrict=__restrict__ -O3 -DGRAPH_GENERATOR_MPI -DREUSE_CSR_FOR_VALIDATION -I../aml
LDFLAGS = -lpthread
MPICC = mpicc
MPICXX = mpicxx

all: graph500_reference_bfs_sssp graph500_reference_bfs graph500_custom_bfs graph500_custom_bfs_sssp

GENERATOR_SOURCES = ../generator/graph_generator.c ../generator/make_graph.c ../generator/splittable_mrg.c ../generator/utils.c
SOURCES = main.c utils.c validate.c ../aml/aml.c
HEADERS = common.h csr_reference.h bitmap_reference.h

graph500_reference_bfs_sssp: bfs_reference.c $(SOURCES) $(HEADERS) $(GENERATOR_SOURCES) csr_reference.c sssp_reference.c
	$(MPICC) $(CFLAGS) -DSSSP $(LDFLAGS) -o graph500_reference_bfs_sssp bfs_reference.c sssp_reference.c csr_reference.c $(SOURCES) $(GENERATOR_SOURCES) -lm $(LDFLAGS)

graph500_reference_bfs: bfs_reference.c $(SOURCES) $(HEADERS) $(GENERATOR_SOURCES) csr_reference.c
	$(MPICC) $(CFLAGS) $(LDFLAGS) -o graph500_reference_bfs bfs_reference.c csr_reference.c $(SOURCES) $(GENERATOR_SOURCES) -lm $(LDFLAGS)

graph500_custom_bfs: bfs_custom.c $(SOURCES) $(HEADERS) $(GENERATOR_SOURCES) csr_reference.c
	$(MPICC) $(CFLAGS) -shared -fPIC utils.c validate.c ../aml/aml.c $(GENERATOR_SOURCES) csr_reference.c -o libbfs_custom.so
	#$(MPICC) $(CFLAGS) -c main.c
	$(MPICXX) $(CFLAGS) $(LDFLAGS) -DUSE_SHMEM=1 -I${HOME}/local/include -I${HCLIB_ROOT}/../modules/bale_actor/inc/ -I${HCLIB_ROOT}/include -I${BALE_INSTALL}/include -o graph500_custom_bfs main.c bfs_custom.c -L./ -lbfs_custom -lm -lpthread -L${HCLIB_ROOT}/lib -lhclib -L${HCLIB_ROOT}/../modules/bale_actor/lib/ -lhclib_bale_actor -L${BALE_INSTALL}/lib -lconvey -llibgetput -L${HOME}/local/lib -lsma
	#$(MPICXX) $(CFLAGS) $(LDFLAGS) -DUSE_SHMEM=1 -I${HOME}/local/include -I${HCLIB_ROOT}/inc/ -I${HCLIB_ROOT}/include -I${BALE_INSTALL}/include -o graph500_custom_bfs main.o bfs_custom.c -L./ -lbfs_custom -lm -lpthread -L${HCLIB_ROOT}/lib -lhclib -L${HCLIB_ROOT}/lib/ -lhclib_actor -L${BALE_INSTALL}/lib -lconvey -llibgetput -L${HOME}/local/lib -lsma -lpmi_simple

graph500_custom_bfs_sssp: bfs_custom.c sssp_custom.c $(SOURCES) $(HEADERS) $(GENERATOR_SOURCES) csr_reference.c
	$(MPICC) $(CFLAGS) $(LDFLAGS) -o graph500_custom_bfs_sssp bfs_custom.c sssp_custom.c csr_reference.c $(SOURCES) $(GENERATOR_SOURCES) -lm $(LDFLAGS)

clean:
	-rm -f graph500_*
