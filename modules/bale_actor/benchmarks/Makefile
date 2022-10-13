include $(HCLIB_ROOT)/../modules/bale_actor/inc/hclib_bale_actor.pre.mak
include $(HCLIB_ROOT)/include/hclib.mak
include $(HCLIB_ROOT)/../modules/bale_actor/inc/hclib_bale_actor.post.mak

CXX ?= CC
UPCC ?= upcc
SRUN ?= srun
BALE_FLAGS ?= -DUSE_SHMEM=1

TARGETS=histo_agi histo_shmem histo_conveyor histo_selector \
        ig_agi ig_shmem ig_conveyor ig_selector \
        toposort_agi toposort_shmem toposort_conveyor toposort_selector \
        triangle_agi triangle_shmem triangle_conveyor triangle_selector \
        randperm_agi randperm_shmem  randperm_conveyor randperm_selector \
        permute_agi permute_shmem permute_conveyor permute_selector \
        transpose_agi transpose_shmem transpose_conveyor transpose_selector \
        jaccard_agi jaccard_shmem jaccard_conveyor jaccard_selector \
        pagerank_pull_LU_agi pagerank_pull_LU_shmem pagerank_pull_LU_conveyor pagerank_pull_LU_selector \
        tricent_agi tricent_shmem tricent_conveyor tricent_selector

TARGETS_UPC=histo_upc ig_upc toposort_upc triangle_upc randperm_upc permute_upc transpose_upc jaccard_upc pagerank_pull_LU_upc tricent_upc

all: $(TARGETS)

upc: $(TARGETS_UPC)

%_agi: %_agi.cpp
	$(CXX) -g -O3 -std=c++11 $(BALE_FLAGS) -I$(BALE_INSTALL)/include -L$(BALE_INSTALL)/lib -o $@ $^ -lspmat -lconvey -lexstack -llibgetput -lm

%_conveyor: %_conveyor.cpp
	$(CXX) -g -O3 -std=c++11 $(BALE_FLAGS) -I$(BALE_INSTALL)/include -L$(BALE_INSTALL)/lib -o $@ $^ -lspmat -lconvey -lexstack -llibgetput -lm
	#$(CXX) -g -O3 -std=c++11 -I $(HPCTOOLKIT_HOME)/include -I$(BALE_INSTALL)/include -L$(BALE_INSTALL)/lib -o $@ $^ -lspmat -lconvey -lexstack -llibgetput -lm -L $(HPCTOOLKIT_HOME)/lib/hpctoolkit -lhpctoolkit

%_selector: %_selector.cpp
	$(CXX) -g -O3 -std=c++11 $(BALE_FLAGS) $(HCLIB_CFLAGS) $(HCLIB_LDFLAGS) -o $@ $^ $(HCLIB_LDLIBS) -lspmat -lconvey -lexstack -llibgetput -lhclib_bale_actor -lm
	#$(CXX) -g -O3 -std=c++11 -I $(HPCTOOLKIT_HOME)/include $(HCLIB_CFLAGS) $(HCLIB_LDFLAGS) -o $@ $^ $(HCLIB_LDLIBS) -lspmat -lconvey -lexstack -llibgetput -lhclib_bale_actor -lm -L $(HPCTOOLKIT_HOME)/lib/hpctoolkit -lhpctoolkit

%_shmem: %_shmem.cpp
	$(CXX) -g -O3 -o $@ $^ -lm

%_upc: %_upc.c
	$(UPCC) -cupc2c -O -o $@ $^ -lm

test:
	for target in $(TARGETS); do $(SRUN) ./$$target; done

clean:
	rm -f $(TARGETS) $(TARGETS_UPC)

