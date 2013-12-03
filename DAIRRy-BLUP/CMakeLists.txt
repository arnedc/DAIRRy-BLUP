project(parreml)

cmake_minimum_required(VERSION 2.6) 

enable_language(Fortran)

find_package(MPI REQUIRED)
find_package(HDF5 REQUIRED)
include_directories(${HDF5_INCLUDE_DIR} ${MPI_INCLUDE_PATH} ${PROJECT_SOURCE_DIR}/src)
include(CMakeForceCompiler)
  CMAKE_FORCE_C_COMPILER(mpicc ${CMAKE_C_COMPILER_ID})
  CMAKE_FORCE_CXX_COMPILER(mpicxx ${CMAKE_CXX_COMPILER_ID})

aux_source_directory(src helpfunctions)
add_executable(parreml ${helpfunctions} main.cpp)

if(MPI_COMPILE_FLAGS)
  set_target_properties(parreml PROPERTIES
    COMPILE_FLAGS "-mpe=mpilog ${MPI_COMPILE_FLAGS}")
endif()

if(MPI_LINK_FLAGS)
  set_target_properties(parreml PROPERTIES
    LINK_FLAGS "-L/opt/intel/mkl/lib/intel64 -L/opt/intel/lib/intel64 ${MPI_LINK_FLAGS}")
endif()

find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)

target_link_libraries(parreml 
                      ${LAPACK_LIBRARIES}
                      ${BLAS_LIBRARIES}
                      ${MPI_LIBRARIES}
		      ${HDF5_LIBRARIES}
		      -lmkl_scalapack_lp64          # dit is de MKL variant van ScaLapack
		      -lmkl_intel_lp64
		      -lmkl_sequential
		      -lmkl_core
                      -lmkl_blacs_intelmpi_lp64     # dit is de MKL variant van BLACS
		      #-liomp5
		      -lm -lpthread)

