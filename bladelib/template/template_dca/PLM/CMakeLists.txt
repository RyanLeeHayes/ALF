cmake_minimum_required(VERSION 3.0)
project(PLM)

find_package(CUDA 9.0)

find_package(OpenMP REQUIRED)
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")

if (CUDA_FOUND)

set(CUDA_NVCC_FLAGS_DEBUG ${CUDA_NVCC_FLAGS_DEBUG} "-G")
set(CUDA_NVCC_FLAGS_RELEASE ${CUDA_NVCC_FLAGS_RELEASE} "--use_fast_math")

file(GLOB_RECURSE sources src/*.cu src/*.h)
# add_definitions(-DNOMPI)
cuda_add_executable(PLM ${sources})
# target_compile_options(md PUBLIC -DNOMPI -g -O3 -arch=sm_20)

else()

file(GLOB_RECURSE sources src/*.c src/*.h)
# add_definitions(-DNOMPI)
add_executable(PLM ${sources})
# target_compile_options(md PUBLIC -DNOMPI -g -O3 -arch=sm_20)
target_link_libraries(PLM m)

endif()
