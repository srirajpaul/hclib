include $(HCLIB_ROOT)/../modules/bale_actor/inc/hclib_bale_actor.pre.mak
include $(HCLIB_ROOT)/include/hclib.mak
include $(HCLIB_ROOT)/../modules/bale_actor/inc/hclib_bale_actor.post.mak

SRUN ?= oshrun

TARGETS=test_conveyor test_selector test_selector3 \
        histo_agi histo_conveyor histo_selector histo_selector_lambda histo_selector_lambda2  histo_agi_selector\
        ig_agi ig_conveyor ig_selector ig_selector_lambda ig_selector_lambda2 ig_agi_selector\
        toposort_agi toposort_conveyor toposort_selector toposort_selector_lambda \
        triangle_agi triangle_conveyor triangle_selector triangle_selector_lambda\
        randperm_agi randperm_conveyor randperm_selector randperm_selector_lambda \
        permute_agi permute_conveyor permute_selector permute_selector_lambda \
        transpose_agi transpose_conveyor transpose_selector transpose_selector_lambda \
        put_selector_shmem put_nbi_shmem \
        get_selector_shmem get_nbi_shmem atomic_inc_selector_shmem atomic_add_selector_shmem

TARGETS_PURE=toposort_pure_selector randperm_pure_selector

all: $(TARGETS)

%: %.cpp
	$(CXX) -g -O3 -std=c++11 -DUSE_SHMEM=1 $(HCLIB_CFLAGS) $(HCLIB_LDFLAGS) -o $@ $^ $(HCLIB_LDLIBS) -lspmat -lconvey -lexstack -llibgetput -lhclib_bale_actor -lm
	#$(CXX) -g -O3 -std=c++11 $(HCLIB_CFLAGS) $(HCLIB_LDFLAGS) -o $@ $^ $(HCLIB_LDLIBS) -llibgetput
	#$(CXX) -g -O3 -std=c++11 -I $(HPCTOOLKIT_HOME)/include $(HCLIB_CFLAGS) $(HCLIB_LDFLAGS) -o $@ $^ $(HCLIB_LDLIBS) -L $(HPCTOOLKIT_HOME)/lib/hpctoolkit/ -lspmat -lconvey -lexstack -llibgetput  -lhclib_bale_actor -lm -lhpctoolkit
	#$(CXX) -g -O3 -std=c++11 -DUSE_LOCK $(HCLIB_CFLAGS) $(HCLIB_LDFLAGS) -o $@_lock $^ $(HCLIB_LDLIBS) -lconvey -lhclib_bale_actor_lock

test:
	for target in $(TARGETS); do $(SRUN) ./$$target; done

clean:
	rm -f $(TARGETS)
