cmake_minimum_required(VERSION 3.16)
project(landy CXX)

add_subdirectory(third_party SYSTEM)

if(CMAKE_CXX_COMPILER_ID MATCHES "(GNU|Clang)")
  add_compile_options(-Werror -Wall -Wextra -Wpedantic -Wshadow -Wconversion -Wsign-conversion)
endif()

add_library(landy INTERFACE)
add_library(landy::landy ALIAS landy)
target_include_directories(landy INTERFACE $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>)
target_compile_features(landy INTERFACE cxx_std_20)

add_executable(landy_example
  src/main.cpp
)
target_link_libraries(landy_example PRIVATE mdspan expected)
target_include_directories(landy_example INTERFACE $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>)
target_compile_features(landy_example INTERFACE cxx_std_20)

if(NOT PROJECT_IS_TOP_LEVEL)
  return()
endif()

include(CTest)
if(BUILD_TESTING)
  add_subdirectory(test)
endif()
