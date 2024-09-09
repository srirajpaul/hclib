#include <iostream>
#include <shmem.h>
#include "selector.h"
#include <map>
#include <string>

// mock format for HTTP message 
// struct MockHTTPMessage {
    // std::string method;
    // std::string path;
    // std::string version;
    // std::map<std::string, std::string> headers;
    // std::string body;
// };
struct packet {
    int64_t id;
    int64_t dest;
    int64_t src;
};

enum MailBoxType{REQUEST, REPLY};

// Dual Mapping from PE role to PE number
std::map<std::string, int> role_to_pe = {
  {"client", 0},
  {"gateway", 1},
  {"service_A", 2},
  {"service_B", 3}
};

class GatewayActor: public hclib::Selector<2, packet> {
    //std::map<int64_t, packet> service_queue;

    void process_request(packet pkt, int sender_rank) {
        //
        // forward to correct PE
        std::cout<<"PE"<<shmem_my_pe()<<": Gateway received request from PE"<<sender_rank<<std::endl;
        //service_queue[pkt.id] = pkt;
        auto rank = role_to_pe["service_A"];
        send(REQUEST, pkt, 2);
        std::cout<<"PE"<<shmem_my_pe()<<": Gateway sent request to Service@PE"<<rank<<std::endl;
    }
    void process_response(packet pkt, int sender_rank) {
        // send to client
        std::cout<<"Gateway received response from PE"<<sender_rank<<std::endl;
        // fetch request from queue
        //service_queue.erase(pkt.id);
        send(REPLY, pkt, role_to_pe["client"]);
        initiate_global_done();
    }
    public:
    void start() {
        std::cout<<"Gateway started"<<std::endl;
        hclib::Selector<2, packet>::start();
    }
    GatewayActor() {
        std::cout<<"PE"<<shmem_my_pe()<<": GatewayActor initialized"<<std::endl;
        mb[REQUEST].process = [this](packet pkt, int sender_rank) { this->process_request(pkt, sender_rank);};
        mb[REPLY].process = [this](packet pkt, int sender_rank) { this->process_response(pkt, sender_rank);};
    }
};

class ClientActor: public hclib::Selector<2, packet> {

    void process_request(packet pkt, int sender_rank) {
        send(REQUEST, pkt, role_to_pe["gateway"]);
    }
    void process_response(packet pkt, int sender_rank) {
        std::cout<<"Client recieved response from PE"<<sender_rank<<std::endl;
        std::cout<<"...should terminate here..."<<std::endl;
        initiate_global_done();
    }
    public:
    void start() {
        std::cout<<"Client started"<<std::endl;
        hclib::Selector<2, packet>::start();
    }
    ClientActor() {
        std::cout<<"PE"<<shmem_my_pe()<<": ClientActor initialized"<<std::endl;
        mb[REQUEST].process = [this](packet pkt, int sender_rank) { this->process_request(pkt, sender_rank);};
        mb[REPLY].process = [this](packet pkt, int sender_rank) { this->process_response(pkt, sender_rank);};
    }
};

class ServiceActor: public hclib::Selector<2, packet> {
    void process_request(packet pkt, int sender_rank) {
        std::cout<<"(request) Service received request from PE"<<sender_rank<<std::endl;
        send(REPLY, pkt, shmem_my_pe());
    }
    void process_reply(packet pkt, int sender_rank) {
        std::cout<<"(reply) Service received request from PE"<<role_to_pe["gateway"]<<std::endl;
        send(REPLY, pkt, role_to_pe["gateway"]);
        initiate_global_done();  
    }
    public:
    void start() {
        std::cout<<"Service started"<<std::endl;
        hclib::Selector<2, packet>::start();
    }
    ServiceActor() {
        std::cout<<"PE"<<shmem_my_pe()<<": ServiceActor initialized"<<std::endl;
        mb[REQUEST].process = [this](packet pkt, int sender_rank) { this->process_request(pkt, sender_rank);};
        mb[REPLY].process = [this](packet pkt, int sender_rank) { this->process_reply(pkt, sender_rank);};
    }
};

// SPMD
int main(int argc, char * argv[]) {
  // Initialize SHMEM
  // shmem_init();

  const char *deps[] = { "system", "bale_actor" };
  hclib::launch(deps, 2, [=] {
    std::cout<<"Starting with "<<shmem_n_pes()<<" PEs"<<std::endl;
    hclib::Selector<2, packet>* current_actor;
    switch(shmem_my_pe()) {
        case 0:
            current_actor = new ClientActor();
            break;
        case 1:
            current_actor = new GatewayActor();
            break;
        case 2:
            current_actor = new ServiceActor();
            break;
    }
    shmem_barrier_all();
    hclib::finish([=]() {
        
        current_actor->start();
        if (shmem_my_pe() == role_to_pe["client"]) {
            packet pkt;
            pkt.id = 1;
            pkt.dest = role_to_pe["service_A"];
            pkt.src = role_to_pe["client"];
            current_actor->send(REQUEST, pkt, shmem_my_pe()); //dummy message so can have (outside) -> PE0 mb[REQUEST] 
        }
        //if (shmem_my_pe() == 0) {
        //    current_actor->done(REQUEST);
        //}
    });
    shmem_barrier_all();
  });
  // Finalize SHMEM
  // shmem_finalize();
  return 0;
}
