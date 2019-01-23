{
  "targets": [
    {
      "target_name": "hclib_openshmem",
      "cflags!": [ "-fno-exceptions" ],
      "cflags_cc!": [ "-fno-exceptions" ],
      #"cflags": [ "-O3", "-fPIC", "-DUSE_OFFLOAD", "-std=c++11" ],
      "cflags": [ "-O3", "-fPIC", "-std=c++11" ],
      "sources": [
        "../src/addon_hclib_node.js_openshmem_util.cpp",
        "../src/addon_hclib_node.js_openshmem_sync.cpp",
        "../src/addon_hclib_node.js_openshmem_async.cpp",
        "../src/addon_hclib_node.js_openshmem.cpp"
      ],
      "include_dirs": [
        "<!@(node -p \"require('node-addon-api').include\")",
        "<!(echo $HCLIB_ROOT)/include",
        "<!(echo $HCLIB_HOME)/modules/system/inc",
        "<!(echo $HCLIB_HOME)/modules/node.js_openshmem/inc",
#        "<!(echo $OPENSHMEM_INSTALL)/include",
      ],
      "link_settings": {
        "libraries": [
          "-lhclib -ldl -lhclib_system -lhclib_node.js_openshmem",
        ],
        "library_dirs": [
          "<!(echo $HCLIB_ROOT)/lib",
          "<!(echo $HCLIB_HOME)/modules/system/lib",
          "<!(echo $HCLIB_HOME)/modules/node.js_openshmem/lib",
#          "<!(echo $OPENSHMEM_INSTALL)/lib",
        ],
      },
      'defines': [ 'NAPI_DISABLE_CPP_EXCEPTIONS' ],
    }
  ]
}
