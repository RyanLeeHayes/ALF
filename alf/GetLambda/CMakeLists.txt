cmake_minimum_required(VERSION 3.0)
project(GetLambda LANGUAGES Fortran)
project(GetSteps LANGUAGES Fortran)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif(NOT CMAKE_BUILD_TYPE)

set(GetLambda_sources
  src/GetLambda.F90
)
set(GetSteps_sources
  src/GetSteps.F90
)
add_executable(GetLambda ${GetLambda_sources})
add_executable(GetSteps ${GetSteps_sources})
