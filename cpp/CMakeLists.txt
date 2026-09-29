# cpp/nmc.cpp is a single self-contained file with its own embedded BAND solver (added in T8).
if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/nmc.cpp)
  add_executable(nmc_cpp nmc.cpp)
  # No fused multiply-add contraction: faithful mode reproduces the original program exactly
  target_compile_options(nmc_cpp PRIVATE
    $<$<CXX_COMPILER_ID:GNU,Clang,AppleClang>:-ffp-contract=off -Wall -Wextra>)
endif()
