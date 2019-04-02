#include "hclib_node.js_openshmem.h"
#include "addon_hclib_node.js_openshmem.h"
#include <shmem.h>

#ifdef USE_OFFLOAD
    #define START_IS_OFFLOAD \
        hclib::finish([=] () {\
        hclib::async_nb_at([=] () {
#else
    #define START_IS_OFFLOAD
#endif

#ifdef USE_OFFLOAD
    #define END_IS_OFFLOAD \
        }, nic);\
        });
#else
    #define END_IS_OFFLOAD
#endif

namespace hclib {

Napi::Number shmem_malloc_sync_fn(const Napi::CallbackInfo& info) {
    Napi::Env env = info.Env();
    int64_t size = info[0].As<Napi::Number>().Int64Value();
    void **out_alloc = (void **)malloc(sizeof(void *));
    START_IS_OFFLOAD
        *out_alloc = shmem_malloc(size);
    END_IS_OFFLOAD
    void *ptr = *out_alloc;
    free(out_alloc);
    return Napi::Number::New(env, (uintptr_t)ptr);
}

Napi::ArrayBuffer shmem_malloc_ab_sync_fn(const Napi::CallbackInfo& info) {
    Napi::Env env = info.Env();
    int64_t size = info[0].As<Napi::Number>().Int64Value();
    void **out_alloc = (void **)malloc(sizeof(void *));
    START_IS_OFFLOAD
        *out_alloc = shmem_malloc(size);
    END_IS_OFFLOAD
    Napi::ArrayBuffer buffer = Napi::ArrayBuffer::New(info.Env(), *out_alloc, size);
    free(out_alloc);
    return buffer;
}

Napi::Number shmem_free_sync_fn(const Napi::CallbackInfo& info) {
    Napi::Env env = info.Env();
    void *ptr = (void*) info[0].As<Napi::Number>().Int64Value();
    START_IS_OFFLOAD
        shmem_free(ptr);
    END_IS_OFFLOAD
    return Napi::Number::New(env, 1);
}

Napi::Number shmem_free_ab_sync_fn(const Napi::CallbackInfo& info) {
    Napi::Env env = info.Env();
    Napi::ArrayBuffer buffer = info[0].As<Napi::ArrayBuffer>();
    void* ptr = static_cast<void*>(buffer.Data());
    START_IS_OFFLOAD
        shmem_free(ptr);
    END_IS_OFFLOAD
    return Napi::Number::New(env, 1);
}

Napi::Number shmem_long_g_sync_fn(const Napi::CallbackInfo& info) {
    Napi::Env env = info.Env();
    long *src = (long*) info[0].As<Napi::Number>().Int64Value();
    int64_t pe = info[1].As<Napi::Number>().Int64Value();
    long *val = (long*)malloc(sizeof(long));
    START_IS_OFFLOAD
        *val = shmem_long_g(src, pe);
    END_IS_OFFLOAD
    long val_tmp = *val;
    free(val);
    return Napi::Number::New(env, val_tmp);
}

Napi::Number shmem_long_p_sync_fn(const Napi::CallbackInfo& info) {
    Napi::Env env = info.Env();
    long *dest = (long*) info[0].As<Napi::Number>().Int64Value();
    long val = info[1].As<Napi::Number>().Int64Value();
    int64_t pe = info[2].As<Napi::Number>().Int64Value();
    START_IS_OFFLOAD
        shmem_long_p(dest, val, pe);
    END_IS_OFFLOAD
    return Napi::Number::New(env, 1);
}

Napi::Number shmem_double_g_sync_fn(const Napi::CallbackInfo& info) {
    Napi::Env env = info.Env();
    double *src = (double*) info[0].As<Napi::Number>().Int64Value();
    int64_t pe = info[1].As<Napi::Number>().Int64Value();
    double *val = (double*)malloc(sizeof(long));
    START_IS_OFFLOAD
        *val = shmem_double_g(src, pe);
    END_IS_OFFLOAD
    double val_tmp = *val;
    free(val);
    return Napi::Number::New(env, val_tmp);
}

Napi::Number shmem_double_p_sync_fn(const Napi::CallbackInfo& info) {
    Napi::Env env = info.Env();
    double *dest = (double*) info[0].As<Napi::Number>().Int64Value();
    double val = info[1].As<Napi::Number>().DoubleValue();
    int64_t pe = info[2].As<Napi::Number>().Int64Value();
    START_IS_OFFLOAD
        shmem_double_p(dest, val, pe);
    END_IS_OFFLOAD
    return Napi::Number::New(env, 1);
}

Napi::Number shmem_ab_double_g_sync_fn(const Napi::CallbackInfo& info) {
    Napi::Env env = info.Env();
    Napi::ArrayBuffer remoteBuffer = info[0].As<Napi::ArrayBuffer>();
    int64_t index = info[1].As<Napi::Number>().Int64Value();
    int64_t remotePE = info[2].As<Napi::Number>().Int64Value();
    double* dataPtr = static_cast<double*>(remoteBuffer.Data());
    double* valPtr = (double*)malloc(sizeof(double));

    START_IS_OFFLOAD
	*valPtr = shmem_double_g(dataPtr + index, remotePE);
    END_IS_OFFLOAD

    double remoteVal = *valPtr;
    free(valPtr);

    return Napi::Number::New(env, remoteVal);
}

Napi::Number shmem_ab_double_p_sync_fn(const Napi::CallbackInfo& info) {
    Napi::Env env = info.Env();
    Napi::ArrayBuffer remoteBuffer = info[0].As<Napi::ArrayBuffer>();
    int64_t index = info[1].As<Napi::Number>().Int64Value();
    double val = info[2].As<Napi::Number>().DoubleValue();
    int64_t remotePE = info[3].As<Napi::Number>().Int64Value();
    double* dataPtr = static_cast<double*>(remoteBuffer.Data());

    START_IS_OFFLOAD
	shmem_double_p(dataPtr + index, val, remotePE);
    END_IS_OFFLOAD

    return Napi::Number::New(env, 1);
}

Napi::Number shmem_double_g_nbi_sync_fn(const Napi::CallbackInfo& info) {
    Napi::Env env = info.Env();
    Napi::ArrayBuffer dest_buffer = info[0].As<Napi::ArrayBuffer>();
    double* dest_ptr = static_cast<double*>(dest_buffer.Data());
    int64_t dest_index = info[1].As<Napi::Number>().Int64Value();
    Napi::ArrayBuffer src_buffer = info[2].As<Napi::ArrayBuffer>();
    double* src_ptr = static_cast<double*>(src_buffer.Data());
    int64_t src_index = info[3].As<Napi::Number>().Int64Value();
    int64_t pe = info[4].As<Napi::Number>().Int64Value();

    START_IS_OFFLOAD
        shmem_double_get_nbi(dest_ptr+dest_index, src_ptr+src_index, 1, pe);
    END_IS_OFFLOAD
    return Napi::Number::New(env, 1);
}

Napi::Number shmem_double_p_nbi_sync_fn(const Napi::CallbackInfo& info) {
    Napi::Env env = info.Env();
    Napi::ArrayBuffer dest_buffer = info[0].As<Napi::ArrayBuffer>();
    double* dest_ptr = static_cast<double*>(dest_buffer.Data());
    int64_t dest_index = info[1].As<Napi::Number>().Int64Value();
    Napi::ArrayBuffer src_buffer = info[2].As<Napi::ArrayBuffer>();
    double* src_ptr = static_cast<double*>(src_buffer.Data());
    int64_t src_index = info[3].As<Napi::Number>().Int64Value();
    int64_t pe = info[4].As<Napi::Number>().Int64Value();

    START_IS_OFFLOAD
        shmem_double_put_nbi(dest_ptr+dest_index, src_ptr+src_index, 1, pe);
    END_IS_OFFLOAD
    return Napi::Number::New(env, 1);
}

Napi::Number shmem_quiet_sync_fn(const Napi::CallbackInfo& info) {
    Napi::Env env = info.Env();
    START_IS_OFFLOAD
        shmem_quiet();
    END_IS_OFFLOAD
    return Napi::Number::New(env, 1);
}

Napi::Number shmem_barrier_all_sync_fn(const Napi::CallbackInfo& info) {
    Napi::Env env = info.Env();
    START_IS_OFFLOAD
        shmem_barrier_all();
    END_IS_OFFLOAD
    return Napi::Number::New(env, 1);
}

Napi::Number shmem_set_lock_sync_fn(const Napi::CallbackInfo& info) {
    Napi::Env env = info.Env();

    //TODO: only one task per rank can lock as of now. otherwise it will cause deadlock.
    //Use promise/future list as done in openshmem module.

    void *ptr = (void*) info[0].As<Napi::Number>().Int64Value();
    START_IS_OFFLOAD
        shmem_set_lock((long*)ptr);
    END_IS_OFFLOAD
    return Napi::Number::New(env, 1);
}

Napi::Number shmem_clear_lock_sync_fn(const Napi::CallbackInfo& info) {
    Napi::Env env = info.Env();
    void *ptr = (void*) info[0].As<Napi::Number>().Int64Value();
    START_IS_OFFLOAD
        shmem_clear_lock((long*)ptr);
    END_IS_OFFLOAD
    return Napi::Number::New(env, 1);
}

Napi::Number shmem_broadcast64_sync_fn(const Napi::CallbackInfo& info) {
    Napi::Env env = info.Env();
    void *dest = (void*) info[0].As<Napi::Number>().Int64Value();
    void *src = (void*) info[1].As<Napi::Number>().Int64Value();
    size_t nelems = info[2].As<Napi::Number>().Int64Value();
    int PE_root = info[3].As<Napi::Number>().Int64Value();
    int PE_start = info[4].As<Napi::Number>().Int64Value();
    int logPE_stride = info[5].As<Napi::Number>().Int64Value();
    int PE_size = info[6].As<Napi::Number>().Int64Value();
    long *pSync = (long*) info[7].As<Napi::Number>().Int64Value();

    START_IS_OFFLOAD
        shmem_broadcast64(dest, src, nelems, PE_root, PE_start, logPE_stride, PE_size, pSync);
    END_IS_OFFLOAD
    return Napi::Number::New(env, 1);
}

Napi::Number shmem_long_sum_to_all_sync_fn(const Napi::CallbackInfo& info) {
    Napi::Env env = info.Env();
    long *dest = (long*) info[0].As<Napi::Number>().Int64Value();
    long *src  = (long*) info[1].As<Napi::Number>().Int64Value();
    int nreduce = info[2].As<Napi::Number>().Int64Value();
    int PE_start = info[3].As<Napi::Number>().Int64Value();
    int logPE_stride = info[4].As<Napi::Number>().Int64Value();
    int PE_size = info[5].As<Napi::Number>().Int64Value();
    long *pWrk = (long*) info[6].As<Napi::Number>().Int64Value();
    long *pSync = (long*) info[7].As<Napi::Number>().Int64Value();

    START_IS_OFFLOAD
        shmem_long_sum_to_all(dest, src, nreduce, PE_start, logPE_stride, PE_size, pWrk, pSync);
    END_IS_OFFLOAD
    return Napi::Number::New(env, 1);
}

Napi::Number shmem_double_sum_to_all_sync_fn(const Napi::CallbackInfo& info) {
    Napi::Env env = info.Env();
    double *dest = (double*) info[0].As<Napi::Number>().Int64Value();
    double *src  = (double*) info[1].As<Napi::Number>().Int64Value();
    int nreduce = info[2].As<Napi::Number>().Int64Value();
    int PE_start = info[3].As<Napi::Number>().Int64Value();
    int logPE_stride = info[4].As<Napi::Number>().Int64Value();
    int PE_size = info[5].As<Napi::Number>().Int64Value();
    double *pWrk = (double*) info[6].As<Napi::Number>().Int64Value();
    long *pSync = (long*) info[7].As<Napi::Number>().Int64Value();

    START_IS_OFFLOAD
        shmem_double_sum_to_all(dest, src, nreduce, PE_start, logPE_stride, PE_size, pWrk, pSync);
    END_IS_OFFLOAD
    return Napi::Number::New(env, 1);
}

Napi::Number shmem_int64_atomic_xor_sync_fn(const Napi::CallbackInfo& info) {
    Napi::Env env = info.Env();
    int64_t *dest = (int64_t*) info[0].As<Napi::Number>().Int64Value();
    int64_t val  = (int64_t) info[1].As<Napi::Number>().Int64Value();
    int64_t pe = info[2].As<Napi::Number>().Int64Value();

    START_IS_OFFLOAD
        shmem_int64_atomic_xor(dest, val, pe);
    END_IS_OFFLOAD
    return Napi::Number::New(env, 1);
}

Napi::Number shmem_int64_atomic_add_sync_fn(const Napi::CallbackInfo& info) {
    Napi::Env env = info.Env();
    int64_t *dest = (int64_t*) info[0].As<Napi::Number>().Int64Value();
    int64_t val  = (int64_t) info[1].As<Napi::Number>().Int64Value();
    int64_t pe = info[2].As<Napi::Number>().Int64Value();

    START_IS_OFFLOAD
        shmem_int64_atomic_add(dest, val, pe);
    END_IS_OFFLOAD
    return Napi::Number::New(env, 1);
}

Napi::Value Init_sync(Napi::Env env, Napi::Object exports) {

  exports.Set(Napi::String::New(env, "malloc_sync"),
              Napi::Function::New(env, shmem_malloc_sync_fn));
  exports.Set(Napi::String::New(env, "malloc_ab_sync"),
              Napi::Function::New(env, shmem_malloc_ab_sync_fn));
  exports.Set(Napi::String::New(env, "free_sync"),
              Napi::Function::New(env, shmem_free_sync_fn));
  exports.Set(Napi::String::New(env, "free_ab_sync"),
              Napi::Function::New(env, shmem_free_ab_sync_fn));
  exports.Set(Napi::String::New(env, "long_g_sync"),
              Napi::Function::New(env, shmem_long_g_sync_fn));
  exports.Set(Napi::String::New(env, "long_p_sync"),
              Napi::Function::New(env, shmem_long_p_sync_fn));
  exports.Set(Napi::String::New(env, "double_g_sync"),
              Napi::Function::New(env, shmem_double_g_sync_fn));
  exports.Set(Napi::String::New(env, "double_p_sync"),
              Napi::Function::New(env, shmem_double_p_sync_fn));
  exports.Set(Napi::String::New(env, "ab_double_g_sync"),
	      Napi::Function::New(env, shmem_ab_double_g_sync_fn));
  exports.Set(Napi::String::New(env, "ab_double_p_sync"),
	      Napi::Function::New(env, shmem_ab_double_p_sync_fn));
  exports.Set(Napi::String::New(env, "double_g_nbi_sync"),
              Napi::Function::New(env, shmem_double_g_nbi_sync_fn));
  exports.Set(Napi::String::New(env, "double_p_nbi_sync"),
              Napi::Function::New(env, shmem_double_p_nbi_sync_fn));
  exports.Set(Napi::String::New(env, "quiet_sync"),
              Napi::Function::New(env, shmem_quiet_sync_fn));
  exports.Set(Napi::String::New(env, "barrier_all_sync"),
              Napi::Function::New(env, shmem_barrier_all_sync_fn));
  exports.Set(Napi::String::New(env, "set_lock_sync"),
              Napi::Function::New(env, shmem_set_lock_sync_fn));
  exports.Set(Napi::String::New(env, "clear_lock_sync"),
              Napi::Function::New(env, shmem_clear_lock_sync_fn));
  exports.Set(Napi::String::New(env, "broadcast64_sync"),
              Napi::Function::New(env, shmem_broadcast64_sync_fn));
  exports.Set(Napi::String::New(env, "long_sum_to_all_sync"),
              Napi::Function::New(env, shmem_long_sum_to_all_sync_fn));
  exports.Set(Napi::String::New(env, "double_sum_to_all_sync"),
              Napi::Function::New(env, shmem_double_sum_to_all_sync_fn));
  exports.Set(Napi::String::New(env, "int64_atomic_xor_sync"),
              Napi::Function::New(env, shmem_int64_atomic_xor_sync_fn));
  exports.Set(Napi::String::New(env, "int64_atomic_add_sync"),
              Napi::Function::New(env, shmem_int64_atomic_add_sync_fn));
  return Napi::Number::New(env, 1);
}

}
