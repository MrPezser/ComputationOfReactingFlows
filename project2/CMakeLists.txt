cmake_minimum_required(VERSION 3.25)
project(project2)

set(CMAKE_CXX_STANDARD 17)

add_executable(project2 main.cpp Indexing.h RHS.cpp RHS.h Solver_main.cpp Solver_main.h Chemistry.cpp Chemistry.h Jacobian.cpp Jacobian.h LinearSolve.cpp LinearSolve.h StateVariables.h ChemSource.cpp ChemSource.h)
