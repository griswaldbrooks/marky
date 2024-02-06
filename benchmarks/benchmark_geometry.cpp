#include "marky/pga.hpp"
#include <iostream>
#include <benchmark/benchmark.h>
#include <random>
#include <tuple>

namespace pga {

struct random_line_generator {
  std::default_random_engine rng{2024};
  std::uniform_real_distribution<double> dist{-1e6, 1e6};
  std::tuple<multivector, multivector> operator()() {
    auto const p1 = point(dist(rng), dist(rng), dist(rng));
    auto const p2 = point(dist(rng), dist(rng), dist(rng));
    auto const p3 = point(dist(rng), dist(rng), dist(rng));
    auto const p4 = point(dist(rng), dist(rng), dist(rng));
    auto const l1 = p1 & p2;
    auto const l2 = p3 & p4;
    return std::make_tuple(l1, l2);   
  }
};

multivector midpoint(multivector const& a, multivector const& b) {
  // compute the meet of the two lines 
  auto const m = a ^ b; 
  // project the meet onto the first line 
  auto const m1 = (m | a) ^ a;
  // project the meet onto the second line 
  auto const m2 = (m | b) ^ a;
  // find the midpoint of the meetss
  return (m1 + m2) * 0.5;
}
static void midpoint_benchmark(benchmark::State& state) {
  auto make_lines = random_line_generator{};
  for (auto _ : state) {
    auto const& [a, b] = make_lines();
    benchmark::DoNotOptimize(midpoint(a, b));
  }
}
// Register the function as a benchmark
BENCHMARK(midpoint_benchmark);
}  // namespace marky::pga
// Run the benchmark
BENCHMARK_MAIN();
