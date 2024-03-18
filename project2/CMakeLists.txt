cmake_minimum_required(VERSION 3.25)
project(project2)

set(CMAKE_CXX_STANDARD 17)

add_executable(project2 main.cpp Indexing.h Flux.cpp Flux.h Solver_main.cpp Solver_main.h Chemistry.cpp Chemistry.h VariableTransform.cpp VariableTransform.h)
