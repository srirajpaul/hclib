-include $(CHARM_HOME)/common.mk
CHARMC = $(CHARM_HOME)/bin/charmc
BINARY = ig

CHARMCFLAGS = $(OPTS) -O3

all: $(BINARY)

ig: ig.o
	$(CHARMC) $(CHARMCFLAGS) -language charm++ -o $@ $+ -module NDMeshStreamer

.SECONDARY: $(patsubst %.C,%.decl.h,$(wildcard *.C))
.SECONDARY: $(patsubst %.C,%.def.h,$(wildcard *.C))

ig.def.h ig.decl.h: ig.ci.stamp

ig.ci.stamp: ig.ci
	$(CHARMC) $(CHARMCFLAGS) $<
	touch $@

ig.o: ig.C ig.decl.h ig.def.h
	$(CHARMC) $(CHARMCFLAGS) -c ig.C

test: $(BINARY)
	$(call run, +p4 ./histo 14 8 )

clean:
	rm -f *.o *.decl.h *.def.h $(BINARY) charmrun* *.stamp
