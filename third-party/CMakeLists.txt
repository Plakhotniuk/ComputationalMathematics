cmake_minimum_required(VERSION 3.17)
project(ComputationalMathematics)
set(CMAKE_CXX_STANDARD 17)


add_library(eigen-3.4.0 INTERFACE)
target_include_directories(eigen-3.4.0 INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})
