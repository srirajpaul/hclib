## Description
This translator is a Clang AST pass, which converts the lambda version of the selector to the class-based version. Currently, it supports 1 mailbox and we will plan to support multiple mailboxes.

## Build & Installation instructions

Step 1: Clone the LLVM repository

```
$ git clone https://github.com/llvm/llvm-project.git
```

Step 2: Create a symbolic link to this directory (`translator`) in `llvm-project/clang-tools-extra`
```
$ cd llvm-project/clang-tools-extra
$ ln -s path_to/translator .
```

Step 3: Update CMakeLists.txt 
```
// in llvm-project/clang-tools-extra
$ emacs CMakeLists.txt
```

Add `add_subdirectory(translator)` after `add_subdirectory(tool-template)`.

Step 4: Build 

Please also see LLVM's official doc: [here](https://llvm.org/docs/GettingStarted.html)

```
$ cd path_to/llvm-project
$ mkdir build
$ cd build
$ cmake -G Ninja -DLLVM_ENABLE_PROJECTS="clang;libcxx;libcxxabi;clang-tools-extra;openmp" -DLLVM_TARGETS_TO_BUILD="X86" -DLLVM_BUILD_TOOLS=OFF -DCMAKE_INSTALL_PREFIX=$PWD/../install -DLLVM_PARALLEL_LINK_JOBS=1 ../llvm
$ ninja install
```

If you want to save your disk, it'd be worth adding `-DCMAKE_BUILD_TYPE="MinSizeRel"`, which produces smallest possible binaries.

Step 5: Test
```
$ cd path_to/hclib/modules/bale_actor/test
$ ../llvm-project/install/bin/selector-trans histo_selector_lambda.cpp -- -I../inc -I$HCLIB_ROOT/include -I$BALE_INSTALL/include
```
