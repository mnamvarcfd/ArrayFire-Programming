cmake_minimum_required(VERSION 3.8)

project(MAIN_MONTECARLO)

set(CMAKE_CXX_STANDARD 17)

find_package(ArrayFire REQUIRED)

add_executable(MAIN_MONTECARLO MAIN_MONTECARLO.cpp
IO/IO.h
IO/console.h
IO/console.cpp
IO/file.h
IO/file.cpp
IO/vtk.h
IO/vtk.cpp
stackoverflow/stackoverflow.h
stackoverflow/stackoverflow.cpp
type/type.h
type/type.cpp)

# To use Unified backend, do the following.
# Unified backend lets you choose the backend at runtime
target_link_libraries(MAIN_MONTECARLO ArrayFire::af)

