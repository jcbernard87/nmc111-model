cmake_minimum_required(VERSION 3.20)
project(nmc_model VERSION 0.1.0 LANGUAGES C CXX Fortran)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# The Fortran program uses bandsolver. By default it is fetched from GitHub;
# use a local checkout with -DFETCHCONTENT_SOURCE_DIR_BANDSOLVER=/path/to/bandsolver
# Build everything, bandsolver included, without fused multiply-add contraction: faithful mode
# reproduces the original programs (built at -O0) bit for bit only then (docs/validation.md).
if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
  string(APPEND CMAKE_Fortran_FLAGS " -ffp-contract=off")
endif()
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
  string(APPEND CMAKE_CXX_FLAGS " -ffp-contract=off")
endif()
if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang")
  string(APPEND CMAKE_C_FLAGS " -ffp-contract=off")
endif()

include(FetchContent)
FetchContent_Declare(bandsolver
  GIT_REPOSITORY https://github.com/jcbernard87/bandsolver.git
  GIT_TAG v0.1.2)
set(BANDSOLVER_BUILD_TESTS OFF CACHE BOOL "" FORCE)
set(BANDSOLVER_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(bandsolver)

enable_testing()
add_subdirectory(fortran)
add_subdirectory(cpp)
