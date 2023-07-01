cmake_minimum_required(VERSION 3.0)
# project(wham LANGUAGES CUDA CXX)
# project(wham LANGUAGES CUDA)
project(wham CXX)

find_package(CUDA 9.0)

if (CUDA_FOUND)

file(GLOB_RECURSE sources src/*.cu src/*.h)
# add_definitions(-DNOMPI)
# cuda_add_executable(wham ${sources})
cuda_add_library(wham SHARED ${sources})
# target_compile_options(md PUBLIC -DNOMPI -g -O3 -arch=sm_20)

else()

file(GLOB_RECURSE sources src/*.c src/*.h)
# add_definitions(-DNOMPI)
# add_executable(wham ${sources})
add_library(wham SHARED ${sources})
# target_compile_options(md PUBLIC -DNOMPI -g -O3 -arch=sm_20)
target_link_libraries(wham m)

endif()
