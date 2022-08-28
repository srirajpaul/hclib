include $(HCLIB_ROOT)/include/hclib.mak

#DEFINES=-DHC_ASSERTION_CHECK -DUSE_CRAY_SHMEM_7 -DUSE_LOCK
#DEFINES=-DHC_ASSERTION_CHECK -DUSE_CRAY_SHMEM_7
DEFINES=-DHC_ASSERTION_CHECK
# DEFINES=-DVERBOSE

OBJS=obj/hclib_bale_actor.o

TARGETS=lib/libhclib_bale_actor.so lib/libhclib_bale_actor.a

all: $(TARGETS)

lib/libhclib_bale_actor.a: $(OBJS)
	$(AR) rvs $@ $(OBJS)

lib/libhclib_bale_actor.so: $(OBJS)
	$(CXX) -shared -fPIC $(HCLIB_LDLIBS) $(HCLIB_LDFLAGS) -o ${CURDIR}/$@ $(OBJS)

obj/hclib_bale_actor.o: src/hclib_bale_actor.cpp
	$(CXX) -g -O3 -c -fPIC -std=c++11 -Iinc $^ -o $@ $(HCLIB_CFLAGS) $(DEFINES)

clean:
	rm -f obj/*.o lib/*

install: $(TARGETS)
	cp lib/*.a lib/*.so $(HCLIB_ROOT)/lib
	#cp inc/*.h $(HCLIB_ROOT)/include/
