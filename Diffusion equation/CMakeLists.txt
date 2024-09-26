cmake_minimum_required(VERSION 3.27)

project(Diffusion)

enable_language(Fortran)

if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
  set (CMAKE_Fortran_FLAGS "-Wall -Wextra -std=f2008")
endif()

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/../run)

add_subdirectory(source)
