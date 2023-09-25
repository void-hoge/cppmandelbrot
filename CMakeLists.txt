cmake_minimum_required(VERSION 3.10)

project(cppmandelbrot CXX)
set(CMAKE_CXX_STANDARD 20)

add_compile_options(-O3 -mtune=native -march=native -DNDEBUG -Wall)

find_package(OpenMP)
if (OpenMP_FOUND)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

add_subdirectory(mandelbrot)
add_subdirectory(buddhabrot)
add_executable(mandel test/main_mandelbrot.cpp)
target_link_libraries(mandel mandelbrot)

add_executable(buddha test/main_buddhabrot.cpp)
target_link_libraries(buddha buddhabrot)
