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

TEST(Multivector, Reverse) {
  auto const a = multivector{1., // scalar
    2., 3., 4., 5., // vector
    6., 7., 8., 9., 10., 11.,// bivector
    12., 13., 14., 15.,// trivector
    16.}; // pseudoscalar
  auto const expected = multivector{1., 2., 3., 4., 5., -6., -7., -8.,
                                    -9., -10., -11., -12., -13., -14., -15., 16.};
  EXPECT_EQ(~a, expected);
}

TEST(Multivector, Dual) {
  auto const a = multivector{1., // scalar
    2., 3., 4., 5., // vector
    6., 7., 8., 9., 10., 11.,// bivector
    12., 13., 14., 15.,// trivector
    16.}; // pseudoscalar
  auto const expected = multivector{16., 15., 14., 13., 12., 11., 10., 9.,
                                    8., 7., 6., 5., 4., 3., 2., 1.};
  EXPECT_EQ(!a, expected);
}

TEST(Multivector, Conjugate) {
  auto const a = multivector{1., // scalar
    2., 3., 4., 5., // vector
    6., 7., 8., 9., 10., 11.,// bivector
    12., 13., 14., 15.,// trivector
    16.}; // pseudoscalar
  auto const expected = multivector{1.,//
    -2., -3., -4., -5., //
    -6., -7., -8., -9., -10., -11., //
    12., 13., 14., 15., //
    16.};
  EXPECT_EQ(conjugate(a), expected);
}

TEST(Multivector, Involute){
  auto const a = multivector{1., // scalar
    2., 3., 4., 5., // vector
    6., 7., 8., 9., 10., 11.,// bivector
    12., 13., 14., 15.,// trivector
    16.}; // pseudoscalar
  auto const expected = multivector{1., // scalar
    -2., -3., -4., -5., // vector
    6., 7., 8., 9., 10., 11.,// bivector
    -12., -13., -14., -15.,// trivector
    16.}; // pseudoscalar
  EXPECT_EQ(involute(a), expected);
}
} // namespace landy::pga
