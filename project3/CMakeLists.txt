cmake_minimum_required(VERSION 3.25)
project(project3)

set(CMAKE_CXX_STANDARD 17)

add_executable(project3 main.cpp Indexing.h Chemistry.cpp Chemistry.h StateVariables.h Solver_main.h LinearSolve.cpp LinearSolve.h Solver_main.cpp)
