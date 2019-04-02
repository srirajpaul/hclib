
//const oshmem = require(process.env.HCLIB_HOME+'/modules/node.js_openshmem/lib/node_hclib_openshmem.node');
const addon = require('./hclib_node.js_openshmem.node');
//const yaml = require('js-yaml');
//const fs = require('fs');

var promFlush;
var ptr_array = new Object();

module.exports = {

init_hc : function (filename) {

  //const type_arr = yaml.safeLoad(fs.readFileSync(filename, 'utf8'));
  //const indentedJson = JSON.stringify(type_arr, null, 4);
  //console.log(indentedJson);

  //for (const one_type of type_arr.arrays) {
    //console.log(one_type);
    //var size = one_type.size * one_type.count;
    //console.log("size " +size);
    //ptr_array[one_type.name] = addon.shmem_malloc(size);
  //}

  addon.init(function(){});
  promFlush = setInterval(() => {
   // addon.clear_put_promises();
  }, 100);
},

finalize_hc : function () {
  addon.finalize();
  clearInterval(promFlush);
},

get_ptr_array : function(name) {
    return ptr_array[name];
},

my_rank: function(){
    return addon.my_rank();
},

num_ranks: function(){
    return addon.num_ranks();
},

get_long_value_at: function(src, index) {
    return addon.get_long_value(src, index);
},

put_long_value_at:function(dest, val, index) {
    return addon.put_long_value(dest, val, index);
},

get_long_value: function(src) {
    return addon.get_long_value(src, 0);
},

put_long_value:function(dest, val) {
    return addon.put_long_value(dest, val, 0);
},

get_double_value_at: function(src, index) {
    return addon.get_double_value(src, index);
},

put_double_value_at:function(dest, val, index) {
    return addon.put_double_value(dest, val, index);
},

get_double_value: function(src) {
    return addon.get_double_value(src, 0);
},

put_double_value:function(dest, val) {
    return addon.put_double_value(dest, val, 0);
},

malloc_sync: function(size) {
    return addon.malloc_sync(size);
},

malloc_ab_sync: function(size) {
    return addon.malloc_ab_sync(size);
},

free_sync: function(ptr) {
    return addon.free_sync(ptr);
},

free_ab_sync: function(array_buffer) {
    return addon.free_ab_sync(array_buffer);
},

long_g_sync: function(src, pe) {
    return addon.long_g_sync(src, pe);
},

long_p_sync: function(dest, val, pe) {
    return addon.long_p_sync(dest, val, pe);
},

double_g_sync: function(src, pe) {
    return addon.double_g_sync(src, pe);
},

double_p_sync: function(dest, val, pe) {
    return addon.double_p_sync(dest, val, pe);
},

ab_double_g_sync: function(dest_buf, index, pe) {
    return addon.ab_double_g_sync(dest_buf, index, pe);
},

ab_double_p_sync: function(dest_buf, index, val, pe) {
    return addon.ab_double_p_sync(dest_buf, index, val, pe);
},

double_g_nbi_sync: function(dest_buf, dest_idx, src_buf, src_idx, pe) {
    return addon.double_g_nbi_sync(dest_buf, dest_idx, src_buf, src_idx, pe);
},

double_p_nbi_sync: function(dest_buf, dest_idx, src_buf, src_idx, pe) {
    return addon.double_p_nbi_sync(dest_buf, dest_idx, src_buf, src_idx, pe);
},

quiet_sync: function() {
    return addon.quiet_sync();
},

barrier_all_sync: function() {
    return addon.barrier_all_sync();
},

set_lock_sync: function(ptr) {
    return addon.set_lock_sync(ptr);
},

clear_lock_sync: function(ptr) {
    return addon.clear_lock_sync(ptr);
},

broadcast64_sync: function(dest, src, nelems, PE_root, PE_start, logPE_stride, PE_size, pSync) {
    return addon.broadcast64_sync(dest, src, nelems, PE_root, PE_start, logPE_stride, PE_size, pSync);
},

long_sum_to_all_sync: function(dest, src, nreduce, PE_start, logPE_stride, PE_size, pWrk, pSync) {
    return addon.long_sum_to_all_sync(dest, src, nreduce, PE_start, logPE_stride, PE_size, pWrk, pSync);
},

double_sum_to_all_sync: function(dest, src, nreduce, PE_start, logPE_stride, PE_size, pWrk, pSync) {
    return addon.double_sum_to_all_sync(dest, src, nreduce, PE_start, logPE_stride, PE_size, pWrk, pSync);
},

int64_atomic_xor_sync: function(dest, val, pe) {
    return addon.int64_atomic_xor_sync(dest, val, pe);
},

int64_atomic_add_sync: function(dest, val, pe) {
    return addon.int64_atomic_add_sync(dest, val, pe);
},

long_g_async: function(src, pe) {
    return addon.long_g_async(src, pe);
},

long_g_async_callback: function(src, pe, cb) {
    return addon.long_g_async_callback(src, pe, cb);
},

double_g_async: function(src, pe) {
    return addon.double_g_async(src, pe);
},

double_g_async_callback: function(src, pe) {
    return addon.double_g_async_callback(src, pe);
},

int64_atomic_add_async: function(dest, val, pe) {
    return addon.int64_atomic_add_async(dest, val, pe);
},

set_lock_async: function(ptr) {
    return addon.set_lock_async(ptr);
},

set_lock_async_callback: function(ptr, cb) {
    return addon.set_lock_async_callback(ptr, cb);
},

clear_lock_async: function(ptr) {
    return addon.clear_lock_async(ptr);
},

clear_lock_async_noreturn: function(ptr) {
    return addon.clear_lock_async_noreturn(ptr);
},

clear_lock_async_finish: function(ptr) {
    return addon.clear_lock_async_finish(ptr);
},

clear_lock_async_callback: function(ptr) {
    return addon.clear_lock_async_callback(ptr);
},

quiet_async: function() {
    return addon.quiet_async();
},

finish: function(func) {
    //TODO: include context to isolate this finish scope from others
    func();
    return addon.quiet_async();
}

};

module.exports.SHMEM_BCAST_SYNC_SIZE = addon.get_SHMEM_BCAST_SYNC_SIZE();
module.exports.SHMEM_BARRIER_SYNC_SIZE = addon.get_SHMEM_BARRIER_SYNC_SIZE();
module.exports.SHMEM_REDUCE_SYNC_SIZE = addon.get_SHMEM_REDUCE_SYNC_SIZE();
module.exports.SHMEM_REDUCE_MIN_WRKDATA_SIZE = addon.get_SHMEM_REDUCE_MIN_WRKDATA_SIZE();
module.exports.SHMEM_SYNC_VALUE = addon.get_SHMEM_SYNC_VALUE();

module.exports.SIZEOF_64_t = 8;
module.exports.at = function(ptr, i) {
    return ptr + module.exports.SIZEOF_64_t*i;
};

