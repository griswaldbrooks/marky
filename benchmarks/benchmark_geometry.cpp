#include "landy/landy.hpp"
#include <iostream>
#include <benchmark/benchmark.h>
#include <random>
#include <tuple>

namespace landy::geometry {

struct random_line_generator {
  std::default_random_engine rng{2024};
  std::uniform_real_distribution<double> dist{-1e6, 1e6};
  std::tuple<line_t, line_t> operator()() {
    auto const l1 = line_t{{dist(rng), dist(rng), dist(rng)}, {dist(rng), dist(rng), dist(rng)}};
    auto const l2 = line_t{{dist(rng), dist(rng), dist(rng)}, {dist(rng), dist(rng), dist(rng)}};
    return std::make_tuple(l1, l2);
  }
};

static void midpoint_benchmark(benchmark::State& state) {
  auto make_lines = random_line_generator{};
  for (auto _ : state) {
    auto const& [a, b] = make_lines();
    benchmark::DoNotOptimize(midpoint(a, b));
  }
}
// Register the function as a benchmark
BENCHMARK(midpoint_benchmark);
}  // namespace landy::geometry
// Run the benchmark
BENCHMARK_MAIN();
