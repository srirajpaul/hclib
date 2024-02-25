#!/bin/bash -x

$LLVM_HOME/bin/selector-trans $1_selector_lambda$2.cpp -- -I/root/local/include -I$BALE_INSTALL/include -I$HCLIB_ROOT/include -I$PWD/../inc | sed 's/#define USE_LAMBDA/#undef USE_LAMBDA/' | $LLVM_HOME/bin/clang-format > $1_selector_lambda$2.trans.cpp
