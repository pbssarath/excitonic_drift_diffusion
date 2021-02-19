cmake_minimum_required(VERSION 3.1.0 FATAL_ERROR)

project(drift_diffusion)

set(DD_SOURCES
  # header files
  include/DDEquation.h
  include/DDExciton.h
  include/DDGridField.h
  include/DDInputData.h
  include/DDNodeData.h
  include/initial_guess.h
  include/lsqr.h
  include/rhs.h
  include/DT.hh
  include/gen_diss.hh
  include/interpolator.h
  # source files to compile
  src/gen_diss.cpp
  src/main.cpp
)

add_executable(DD ${DD_SOURCES})

# add the ./include directory to the include search path so we can #include
# the header files we put in include/
target_include_directories(DD PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)


find_package(talyfem REQUIRED)

# Add reproducibility information
include (TalyReproducibility)
add_taly_reproducibility(TARGET DD DEPENDS ${DD_SOURCES})

include (TalyDoxygen)
add_doxygen(TARGET DD INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/include)

# this will also copy all the compiler definitions, etc. that TalyFEM was
# compiled with (that were marked as PUBLIC) - e.g. this will automatically
# add "-DENABLE_HDF5" if TalyFEM was compiled with HDF5 support
target_link_libraries(DD talyfem)

# automatically generate documentation with Doxygen when running 'make docs'
# requires us to have the TalyDoxygen CMake module
# add_doxygen(TARGET DD INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/include)
