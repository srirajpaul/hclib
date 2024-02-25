#module swap PrgEnv-intel PrgEnv-gnu; module load cray-shmem bupc openmpi

CXX ?= CC
OSHCC ?=oshcc
#OSHCC ?=CC
UPCC ?=upcc
MPICC ?=mpicc
#MPICC ?=CC
SRUN ?=srun

TARGETS=histo_shmem histo_conveyor histo_upc histo_mpi3\
		ig_shmem ig_nbi_shmem ig_conveyor ig_upc ig_nbi_upc ig_mpi3

TARGETS_CHARM=histo_charm++ ig_charm++

all: $(TARGETS) $(TARGETS_CHARM)

%_shmem: %_shmem.cpp
	$(OSHCC) -g -O3 -o $@ $^

histo_mpi3: histo_mpi3.cpp
	#$(MPICC) -g -O3 -DUSE_FETCH_AND_OP -o histo_mpi3.fa $^
	#$(MPICC) -g -O3 -DUSE_GET_ACCUMULATE -o histo_mpi3.gaccum $^
	$(MPICC) -g -O3 -DUSE_ACCUMULATE -o histo_mpi3 $^
	$(CXX) -g -O3 -DUSE_ACCUMULATE -o histo_cc_mpi3 $^

ig_mpi3: ig_mpi3.cpp
	$(MPICC) -g -O3 -o ig_mpi3 $^
	$(CXX) -g -O3 -o ig_cc_mpi3 $^

%_conveyor: %_conveyor.cpp
	$(OSHCC) -g -O3 -I$(BALE_INSTALL)/include -L$(BALE_INSTALL)/lib -o $@ $^ -lspmat -lconvey -lexstack -llibgetput

%_upc: %_upc.c
	$(UPCC) -cupc2c -O -o $@ $^

.PHONY: $(TARGETS_CHARM)
$(TARGETS_CHARM):
	make -C $@
	#cd $@ && make

test:
	for target in $(TARGETS); do $(SRUN) ./$$target ; done

clean:
	make -C histo_charm++ clean && make -C ig_charm++ clean && rm -f $(TARGETS)

