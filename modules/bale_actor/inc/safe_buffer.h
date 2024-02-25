
#ifndef SAFE_BUFFER_H
#define SAFE_BUFFER_H

#include "boost/circular_buffer.hpp"
#include <mutex>

namespace hclib {
namespace conveyor {

template<typename T>
class safe_buffer {
    
    boost::circular_buffer<T> cb;
    std::mutex mtx;

  public:

    safe_buffer() {
#ifdef SELECTOR_DEBUG
        printf("creating safe buffer\n");
#endif
    }

    ~safe_buffer() {
#ifdef SELECTOR_DEBUG
        printf("deleting safe buffer\n");
#endif        
    }

    safe_buffer(size_t capacity) {
#ifdef SELECTOR_DEBUG
        printf("creating safe buffer with capacity\n");
#endif
        cb.set_capacity(capacity);
    }

    void set_capacity(size_t capacity) {
        cb.set_capacity(capacity);
    }

    std::mutex& get_mutex() {
        return mtx;
    }

    void push_back(const T& val) {
#ifdef USE_LOCK
        std::lock_guard<std::mutex> lg(mtx); 
#endif
        cb.push_back(val);
    }

    T& operator [](size_t index) {
        return cb[index];
    }

    T& at(size_t index) {
        return cb.at(index);
    }

    void erase_begin(size_t index) {
        cb.erase_begin(index);
    }

    size_t size() const noexcept {
        return cb.size();
    }

    bool full() const noexcept {
        return cb.full();
    }
}; // class safe_buffer

}; // namespace conveyor
}; // namespace hclib
#endif // SAFE_BUFFER_H

