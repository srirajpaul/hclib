-include $(CHARM_HOME)/common.mk
CHARMC = $(CHARM_HOME)/bin/charmc
BINARY = histo

CHARMCFLAGS = $(OPTS) -O3

all: $(BINARY)

histo: histo.o
	$(CHARMC) $(CHARMCFLAGS) -language charm++ -o $@ $+ -module NDMeshStreamer

.SECONDARY: $(patsubst %.C,%.decl.h,$(wildcard *.C))
.SECONDARY: $(patsubst %.C,%.def.h,$(wildcard *.C))

histo.def.h histo.decl.h: histo.ci.stamp

histo.ci.stamp: histo.ci
	$(CHARMC) $(CHARMCFLAGS) $<
	touch $@

histo.o: histo.C histo.decl.h histo.def.h
	$(CHARMC) $(CHARMCFLAGS) -c histo.C

test: $(BINARY)
	$(call run, +p4 ./histo 14 8 )

clean:
	rm -f *.o *.decl.h *.def.h $(BINARY) charmrun* *.stamp
