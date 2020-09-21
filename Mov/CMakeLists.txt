cmake_minimum_required(VERSION 3.17)
project(Mov)

set(CMAKE_CXX_STANDARD 14)

add_executable(Mov main.cpp utility_functions.h onepass.h offline.h fantom.h samplegreedy.h multipass.h)