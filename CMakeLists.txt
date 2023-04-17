cmake_minimum_required(VERSION 3.0.0)
project(launch VERSION 0.1.0)

add_executable(
  launch
  main.cpp

  sources/LocationData.cpp
)

target_include_directories(
  launch PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/includes
)

