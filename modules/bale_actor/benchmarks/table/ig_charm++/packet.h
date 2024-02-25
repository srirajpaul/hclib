#ifndef _PACKET_H_
#define _PACKET_H_
class packet1 {
 public:
  CmiInt8 val;
  CmiInt8 idx;
  CmiInt8 pe;

  void pup(PUP::er &p) {
    // remember to pup your superclass if there is one
    p|val;
    p|idx;
    p|pe;
  }
};
class packet2 {
 public:
  CmiInt8 val;
  CmiInt8 idx;

  void pup(PUP::er &p) {
    // remember to pup your superclass if there is one
    p|val;
    p|idx;
  }
};
#endif
