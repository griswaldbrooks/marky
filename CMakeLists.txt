cmake_minimum_required(VERSION 3.22)
project(marky CXX)

add_subdirectory(third_party SYSTEM)

if(CMAKE_CXX_COMPILER_ID MATCHES "(GNU|Clang)")
  add_compile_options(-Werror -Wall -Wextra -Wpedantic -Wshadow -Wconversion -Wsign-conversion)
endif()

add_library(marky INTERFACE)
add_library(marky::marky ALIAS marky)
target_include_directories(marky INTERFACE $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>)
target_compile_features(marky INTERFACE cxx_std_20)
target_link_libraries(marky INTERFACE mp-units::mp-units)

add_executable(marky_example
  src/main.cpp
)
target_link_libraries(marky_example PRIVATE mdspan expected)
target_include_directories(marky_example INTERFACE $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>)
target_compile_features(marky_example INTERFACE cxx_std_20)

if(NOT PROJECT_IS_TOP_LEVEL)
  return()
endif()

include(CTest)
if(BUILD_TESTING)
  add_subdirectory(test)
endif()

add_subdirectory(benchmarks)
