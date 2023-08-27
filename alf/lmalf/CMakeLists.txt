cmake_minimum_required(VERSION 3.0)
project(whamweight LANGUAGES CXX)
project(lmalf LANGUAGES CXX)

find_package(CUDA 9.0)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif(NOT CMAKE_BUILD_TYPE)

set(whamweight_sources
  src/whamweight.cu
  src/whamweight.h
)
set(lmalf_sources
  src/lmalf.cu
  src/lmalf.h
)

if (CUDA_FOUND)

set(CUDA_NVCC_FLAGS_DEBUG ${CUDA_NVCC_FLAGS_DEBUG} "-G")
set(CUDA_NVCC_FLAGS_RELEASE ${CUDA_NVCC_FLAGS_RELEASE} "--use_fast_math")

# cuda_add_library(whamweight SHARED ${whamweight_sources})
cuda_add_executable(whamweight ${whamweight_sources})
# cuda_add_library(lmalf SHARED ${lmalf_sources})
cuda_add_executable(lmalf ${lmalf_sources})

else()

message(FATAL_ERROR "No CUDA found. LMALF requires CUDA")

endif()
