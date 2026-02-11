cmake_minimum_required(VERSION 3.20)
# loss
project(linear2018 LANGUAGES CXX CUDA)
project(nonlinear2024 LANGUAGES CXX CUDA)
project(whamweight LANGUAGES CXX CUDA)
project(impcons LANGUAGES CXX CUDA)
# dca
project(Filter LANGUAGES C)
project(Moment LANGUAGES C)
project(LM LANGUAGES C)
project(PLM LANGUAGES CXX CUDA)
# io
project(GetLambda LANGUAGES Fortran)
project(GetSteps LANGUAGES Fortran)

# Note: added CUDA to all project commands, and changed cuda_add_library to add_library and cuda_add_executable to add_executable for cmake 4.0 compatibility

find_package(CUDA 9.0 REQUIRED)

find_package(OpenMP REQUIRED)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")

find_package(MPI REQUIRED)
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${MPI_C_LIBRARIES} ${MPI_C_LINK_FLAGS}")

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif(NOT CMAKE_BUILD_TYPE)

# loss
set(linear2018_sources
  loss/linear2018.cu
  loss/linear2018.h
)
set(nonlinear2024_sources
  loss/nonlinear2024.cu
  loss/nonlinear2024.h
)
set(whamweight_sources
  loss/whamweight.cu
  loss/whamweight.h
)
set(impcons_sources
  loss/impcons.cu
  loss/impcons.h
)
# dca
set(filter_sources
  dca/Filter.c
  dca/Filter.h
)
set(moment_sources
  dca/ObservedMoments.c
  dca/ObservedMoments.h
)
set(lm_sources
  dca/LM.c
  dca/LM.h
)
set(plm_sources
  dca/PLM.cu
  dca/PLM.h
)
# io
set(GetLambda_sources
  io/GetLambda.F90
)
set(GetSteps_sources
  io/GetSteps.F90
)

# dca
add_executable(Filter ${filter_sources})
add_executable(Moment ${moment_sources})
add_executable(LM ${lm_sources})
target_link_libraries(Filter mpi)
target_link_libraries(Moment mpi)
target_link_libraries(LM m mpi)
target_include_directories(Filter PUBLIC ${MPI_C_INCLUDE_PATH} .)
target_include_directories(Moment PUBLIC ${MPI_C_INCLUDE_PATH} .)
target_include_directories(LM PUBLIC ${MPI_C_INCLUDE_PATH} .)
# io
add_executable(GetLambda ${GetLambda_sources})
add_executable(GetSteps ${GetSteps_sources})

set(CUDA_NVCC_FLAGS_DEBUG ${CUDA_NVCC_FLAGS_DEBUG} "-G")
set(CUDA_NVCC_FLAGS_RELEASE ${CUDA_NVCC_FLAGS_RELEASE} "--use_fast_math")

# loss
add_library(linear2018 SHARED ${linear2018_sources})
# add_executable(linear2018 ${linear2018_sources})
# add_library(nonlinear2024 SHARED ${nonlinear2024_sources})
add_executable(nonlinear2024 ${nonlinear2024_sources})
# add_library(whamweight SHARED ${whamweight_sources})
add_executable(whamweight ${whamweight_sources})
# add_library(impcons SHARED ${impcons_sources})
add_executable(impcons ${impcons_sources})
# dca
# add_library(PLM SHARED ${plm_sources})
add_executable(PLM ${plm_sources})
