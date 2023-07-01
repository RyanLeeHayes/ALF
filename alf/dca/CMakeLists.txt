cmake_minimum_required(VERSION 3.0)
project(Filter LANGUAGES C)
project(Moment LANGUAGES C)
project(LM LANGUAGES C)
project(PLM LANGUAGES CXX)

find_package(CUDA 9.0)

find_package(OpenMP REQUIRED)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")

find_package(MPI REQUIRED)
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${MPI_C_LIBRARIES} ${MPI_C_LINK_FLAGS}")

set(filter_sources
  src/Filter.c
  src/Filter.h
)
set(moment_sources
  src/ObservedMoments.c
  src/ObservedMoments.h
)
set(lm_sources
  src/LM.c
  src/LM.h
)
set(plm_sources
  src/PLM.cu
  src/PLM.h
)

# add_library(Filter SHARED ${filter_sources})
# add_library(Moment SHARED ${moment_sources})
# add_library(LM SHARED ${lm_sources})
add_executable(Filter ${filter_sources})
add_executable(Moment ${moment_sources})
add_executable(LM ${lm_sources})
target_link_libraries(LM m)
target_include_directories(Filter PUBLIC ${MPI_C_INCLUDE_PATH} .)
target_include_directories(Moment PUBLIC ${MPI_C_INCLUDE_PATH} .)
target_include_directories(LM PUBLIC ${MPI_C_INCLUDE_PATH} .)

if (CUDA_FOUND)

set(CUDA_NVCC_FLAGS_DEBUG ${CUDA_NVCC_FLAGS_DEBUG} "-G")
set(CUDA_NVCC_FLAGS_RELEASE ${CUDA_NVCC_FLAGS_RELEASE} "--use_fast_math")

# add_library(PLM SHARED ${plm_sources})
cuda_add_executable(PLM ${plm_sources})

else()

message(FATAL_ERROR "No CUDA found. PLMGPU is hardcoded on in src/PLM.cu, which requires CUDA")

# add_library(PLM SHARED ${plm_sources})
add_executable(PLM ${plm_sources})

endif()
