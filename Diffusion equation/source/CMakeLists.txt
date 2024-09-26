include_directories(${CMAKE_BINARY_DIR}/source/)

add_executable(program.out
  kinds.f08
  math.f08
  schemes.f08
  input_output.f08
  main.f08
  )
