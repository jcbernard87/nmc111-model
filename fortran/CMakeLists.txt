# fortran/nmc.f90: single-file program (added in T7).
if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/nmc.f90)
  add_executable(nmc_f nmc.f90)
  target_link_libraries(nmc_f PRIVATE bandsolver_f)
  target_include_directories(nmc_f PRIVATE ${CMAKE_BINARY_DIR}/modules)
  # No fused multiply-add contraction: faithful mode reproduces the original program exactly
  target_compile_options(nmc_f PRIVATE $<$<Fortran_COMPILER_ID:GNU>:-ffp-contract=off -std=f2008 -Wall>)
endif()
