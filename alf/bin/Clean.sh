#! /bin/bash
# Remove files created by cmake to begin a clean compilation attempt

rm -r CMakeCache.txt CMakeFiles cmake_install.cmake Makefile
# from loss
rm liblinear2018.so nonlinear2024 whamweight impcons
# from dca
rm PLM LM Filter Moment
# from io
rm GetLambda GetSteps
