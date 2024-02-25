set (LLVM_LINK_COMPONENTS Support)

add_clang_tool(selector-trans SelectorTrans.cpp)

target_link_libraries(selector-trans
  PRIVATE
  clangAST
  clangBasic
  clangDriver
  clangFrontend
  clangRewriteFrontend
  clangSerialization
  clangTooling
  )

install(TARGETS selector-trans RUNTIME DESTINATION bin)
