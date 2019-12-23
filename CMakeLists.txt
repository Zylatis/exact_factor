cmake_minimum_required(VERSION 2.8)
set(CMAKE_CXX_STANDARD 11)
set(
    CMAKE_RUNTIME_OUTPUT_DIRECTORY
    ${CMAKE_HOME_DIRECTORY}
    )
project(ef_investigation C CXX)
FIND_PACKAGE ( OpenMP  REQUIRED )
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -fopenmp ")
add_executable(solver src/solver.cpp)
add_executable(unit_tests tests/unit_tests.cpp)
#target_link_libraries(solver fftw3)
#target_link_libraries(unit_tests fftw3)
	
