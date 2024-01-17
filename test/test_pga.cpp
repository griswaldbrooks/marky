#include <numeric>
#include <gtest/gtest.h>
#include "landy/pga.hpp"
namespace pga {

multivector make_zero_multivector() {
  auto m = multivector{};
  for (auto& element : m) {
    element = 0.;
  }
  return m;
}

multivector brute_force_add(multivector const& a, multivector const& b) {
  auto result = make_zero_multivector();
  for (std::size_t i = 0; i < 16; ++i) {
    result[i] = a[i] + b[i];
  }
  return result;
}

void fill_increasing(multivector& m) {
  std::iota(m.begin(), m.end(), 0.);
}

void fill_decreasing(multivector& m) {
  std::iota(m.rbegin(), m.rend(), 0.);
}

TEST(Multivector, ZeroInitialized) {
 auto const a = multivector{};
  auto const expected = make_zero_multivector();
  EXPECT_EQ(a, expected);
}

TEST(Multivector, Addition) {
  auto const a = e0;
  auto const b = e1;
  auto const result = a + b;
  auto const expected = brute_force_add(a, b);
  EXPECT_EQ(result, expected);
}

// TEST(Multivector, Reverse) {
//   auto const a = make_zero_multivector();
//   fill_increasing(a);
// }
} // namespace landy::pga
