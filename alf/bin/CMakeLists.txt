cmake_minimum_required(VERSION 3.20)
project(ALF LANGUAGES C CXX CUDA Fortran)

find_package(OpenMP REQUIRED)
find_package(MPI REQUIRED)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

# --- loss (CUDA) ---

add_library(linear2018 SHARED loss/linear2018.cu loss/linear2018.h)
target_compile_options(linear2018 PRIVATE
  $<$<CONFIG:Release>:--use_fast_math>
  $<$<CONFIG:Debug>:-G>)

add_executable(nonlinear2024 loss/nonlinear2024.cu loss/nonlinear2024.h)
target_compile_options(nonlinear2024 PRIVATE
  $<$<CONFIG:Release>:--use_fast_math>
  $<$<CONFIG:Debug>:-G>)

add_executable(whamweight loss/whamweight.cu loss/whamweight.h)
target_compile_options(whamweight PRIVATE
  $<$<CONFIG:Release>:--use_fast_math>
  $<$<CONFIG:Debug>:-G>)

add_executable(impcons loss/impcons.cu loss/impcons.h)
target_compile_options(impcons PRIVATE
  $<$<CONFIG:Release>:--use_fast_math>
  $<$<CONFIG:Debug>:-G>)

# --- dca (C with MPI, CUDA) ---

add_executable(Filter dca/Filter.c dca/Filter.h)
target_link_libraries(Filter PRIVATE OpenMP::OpenMP_C MPI::MPI_C)

add_executable(Moment dca/ObservedMoments.c dca/ObservedMoments.h)
target_link_libraries(Moment PRIVATE OpenMP::OpenMP_C MPI::MPI_C)

add_executable(LM dca/LM.c dca/LM.h)
target_link_libraries(LM PRIVATE m OpenMP::OpenMP_C MPI::MPI_C)

add_executable(PLM dca/PLM.cu dca/PLM.h)
target_link_libraries(PLM PRIVATE OpenMP::OpenMP_CXX)
target_compile_options(PLM PRIVATE
  $<$<CONFIG:Release>:--use_fast_math>
  $<$<CONFIG:Debug>:-G>)

# --- io (Fortran) ---

add_executable(GetLambda io/GetLambda.F90)
add_executable(GetSteps io/GetSteps.F90)
