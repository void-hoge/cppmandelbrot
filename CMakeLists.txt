cmake_minimum_required(VERSION 3.10)

project(cppmandelbrot CXX)
set(CMAKE_CXX_STANDARD 20)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -mtune=native -march=native -DNDEBUG")

find_package(OpenMP)
if (OpenMP_FOUND)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

add_executable(mandel test/main_mandelbrot.cpp mandelbrot/mandelbrot.cpp)
target_link_libraries(mandel gmp gmpxx)

add_executable(buddha test/main_buddhabrot.cpp mandelbrot/mandelbrot.cpp buddhabrot/buddhabrot.cpp)
target_link_libraries(buddha gmp gmpxx)
