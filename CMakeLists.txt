cmake_minimum_required(VERSION 3.9)
project(WaveMarching)

set(CMAKE_CXX_STANDARD 11)

add_executable(WaveMarching
        no_opencv_version/2D/wave_mm_2D.cpp
        no_opencv_version/TYPES/WMMStructs.h)
