
#ifndef SELECTOR_SHMEM_H
#define SELECTOR_SHMEM_H

#include<selector.h>

template<typename T>
struct PutPkt {
    T *loc;
    T val;
};

template<typename T>
struct AtomicIncPkt {
    T *loc;
};

template<typename T>
using AtomicAddPkt = PutPkt<T>;

template<typename T>
struct GetPkt {
    T *dest;
    union {
        T *src;
        T val;
    };
};

enum MailBoxType{REQUEST, RESPONSE};

template<typename T, typename P=PutPkt<T>>
class Put_nbi : public hclib::Actor<P> {

  void process(P pkt, int sender_rank) {
      *(pkt.loc) = pkt.val;
  }

  public:
    Put_nbi() {
        hclib::Actor<P>::mb[0].process = [this](P pkt, int sender_rank) { this->process(pkt, sender_rank);};
        //hclib::Actor<P>::start();
    }

    void operator()(T *dest, T *src, size_t nelems,  int pe) {
        for(int i=0;i<nelems;i++) {
            hclib::Actor<P>::send({dest+i, src[i]}, pe);
        }
        //hclib::Actor<P>::send({loc, val}, pe);
    }
};

template<typename T, typename P=AtomicIncPkt<T>>
class AtomicInc : public hclib::Actor<P> {

  void process(P pkt, int sender_rank) {
      *(pkt.loc) = *(pkt.loc) + 1;
  }

  public:
    AtomicInc() {
        hclib::Actor<P>::mb[0].process = [this](P pkt, int sender_rank) { this->process(pkt, sender_rank);};
        //hclib::Actor<P>::start();
    }

    void operator()(T *dest, int pe) {
        hclib::Actor<P>::send({dest}, pe);
    }
};

template<typename T, typename P=AtomicAddPkt<T>>
class AtomicAdd : public hclib::Actor<P> {

  void process(P pkt, int sender_rank) {
      *(pkt.loc) = *(pkt.loc) + pkt.val;;
  }

  public:
    AtomicAdd() {
        hclib::Actor<P>::mb[0].process = [this](P pkt, int sender_rank) { this->process(pkt, sender_rank);};
        //hclib::Actor<P>::start();
    }

    void operator()(T *dest, T val,  int pe) {
        hclib::Actor<P>::send({dest, val}, pe);
    }
};

template<typename T, typename P=GetPkt<T>>
class Get_nbi: public hclib::Selector<2, P>{

  void req_process(P pkt, int sender_rank) {
      pkt.val = *(pkt.src);
      hclib::Selector<2, P>::send(RESPONSE, pkt, sender_rank);
  }

  void resp_process(P pkt, int sender_rank) {
      *(pkt.dest) = pkt.val;
  }

  public:

    Get_nbi() {
        hclib::Selector<2, P>::mb[REQUEST].process = [this](P pkt, int sender_rank) { this->req_process(pkt, sender_rank); };
        hclib::Selector<2, P>::mb[RESPONSE].process = [this](P pkt, int sender_rank) { this->resp_process(pkt, sender_rank); };
        //hclib::Selector<2, P>::start();
    }

    void operator()(T *dest, T *src, size_t nelems, int pe) {
        for(int i=0;i<nelems;i++) {
            hclib::Selector<2, P>::send(REQUEST,  {dest+i, src+i}, pe);
        }
        //hclib::Selector<2, P>::send(REQUEST, {dest, src}, pe);
    }

    void done() {
        hclib::Selector<2, P>::done(REQUEST);
    }
};

#endif

