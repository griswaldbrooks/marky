find_package(benchmark REQUIRED)

add_executable(benchmark_geometry benchmark_geometry.cpp)
target_link_libraries(benchmark_geometry PRIVATE landy::landy benchmark::benchmark)

install(
  TARGETS benchmark_geometry
  LIBRARY DESTINATION lib
  RUNTIME DESTINATION bin)
