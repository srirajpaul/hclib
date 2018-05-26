#ifndef HCLIB_RESILIENCE_H
#define HCLIB_RESILIENCE_H

#include "hclib_cpp.h"
#include "hclib_ref-count.h"
#define USE_STD_VEC
//#define USE_C_ARRAY
#ifdef USE_STD_VEC
#include <vector>
#include <mutex>
#else
#include "hclib-atomics.h"
#endif

namespace hclib {
namespace resilience {

/*
Resilient object.
Users should extend this class to provide their own
equals method to be used for comparison at the end
of resilient task
*/
class obj {
  public:
    virtual ~obj() {}

    //return true if objs are equal else false
    virtual bool equals(obj* obj2) { assert(false); return false; }
};

namespace util {
/*
vector that allows concurrent insertion
Read operations are not concurrent with insertion
*/
template<typename T>
class safe_vector {
#ifdef USE_STD_VEC
    std::vector<T> vec;
    std::mutex mtx;

  public:
    //TODO: can this take both lvalue and rvalue?
    void push_back(T&& value) {
        std::lock_guard<std::mutex> lkg(mtx);
        vec.push_back(std::forward<T>(value));
    }

    size_t size() const noexcept {
        return vec.size();
    }

    T* data() noexcept {
        return vec.data();
    }

    void clear() noexcept {
        vec.clear();
    }

    //auto begin() noexcept -> decltype(vec.begin()) {
    //    return vec.begin();
    //}

    //auto end() noexcept -> decltype(vec.end()) {
    //    return vec.end();
    //}
#elif defined USE_C_ARRAY
#ifndef SAFE_VECTOR_CAPACITY
#define SAFE_VECTOR_CAPACITY 128
#endif
    T vec[SAFE_VECTOR_CAPACITY];
    volatile int pos = -1;
#ifdef USE_MUTEX_LOCK
    std::mutex mtx;
#endif

  public:
    //TODO: can this take both lvalue and rvalue?
    void push_back(T&& value) {
#ifdef USE_MUTEX_LOCK
	std::lock_guard<std::mutex> lkg(mtx);
        pos++;
        assert(pos < SAFE_VECTOR_CAPACITY);
        vec[pos] = value;
#else
	hc_atomic_inc(&pos);
	assert(pos < SAFE_VECTOR_CAPACITY);
	vec[pos] = std::forward<T>(value);
#endif
    }

    size_t size() const noexcept {
        return pos + 1;
    }

    T* data() noexcept {
        return &vec[0];
    }

    void clear() noexcept {
	// TODO: mutual exclusion?
        pos = -1;
    }
#endif
};

template<typename T>
class safe_promise_vector : public safe_vector<T> {
  public:
    void do_puts(int index) {
        //perform the actual put
        auto data_var = safe_vector<T>::data();
        auto size_var = safe_vector<T>::size();
        for(int i=0; i<size_var; i++)
            data_var[i]->put_actual(index);
    }
};

template<typename T>
class safe_future_vector : public safe_vector<T> {
  public:
    void do_releases() {
        //perform the release
        auto data_var = safe_vector<T>::data();
        auto size_var = safe_vector<T>::size();
        for(int i=0; i<size_var; i++)
            data_var[i]->release();
    }
};

} // namespace util

} // namespace resilience
} // namespace hclib

#include "hclib_resilience_promise.h"
#include "hclib_resilience_diamond.h"
#include "hclib_resilience_replay.h"
#include "hclib_resilience_abft.h"

#endif

