#! /bin/bash
# Remove files created by cmake to begin a clean compilation attempt

rm -r CMakeCache.txt CMakeFiles cmake_install.cmake Makefile
# from loss
rm linear2018 nonlinear2024 whamweight impcons
# from dca
rm PLM LM Filter Moment
# from io
rm GetLambda GetSteps
