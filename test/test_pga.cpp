#include <string>
#include <numeric>
#include <ostream>
#include <gtest/gtest.h>
#include "landy/pga.hpp"
namespace pga {

TEST(Multivector, ZeroInitialized) {
 auto const a = multivector{};
  auto const expected = multivector{0.,
    0., 0., 0., 0.,
    0., 0., 0., 0., 0., 0.,
    0., 0., 0., 0.,
    0.};
  EXPECT_EQ(a, expected);
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

TEST(Multivector, Addition){
  auto const a = multivector{1., // scalar
    2., 3., 4., 5., // vector
    6., 7., 8., 9., 10., 11.,// bivector
    12., 13., 14., 15.,// trivector
    16.}; // pseudoscalar
  auto const b = multivector{1., // scalar
    1., 1., 1., 1., // vector
    1., 1., 1., 1., 1., 1.,// bivector
    1., 1., 1., 1.,// trivector
    1.}; // pseudoscalar
  auto const expected = multivector{2., // scalar
    3., 4., 5., 6., // vector
    7., 8., 9., 10., 11., 12.,// bivector
    13., 14., 15., 16.,// trivector
    17.}; // pseudoscalar
  EXPECT_EQ(a + b, expected);
}

TEST(Multivector, Subtraction){
  auto const a = multivector{1., // scalar
    2., 3., 4., 5., // vector
    6., 7., 8., 9., 10., 11.,// bivector
    12., 13., 14., 15.,// trivector
    16.}; // pseudoscalar
  auto const b = multivector{1., // scalar
    1., 1., 1., 1., // vector
    1., 1., 1., 1., 1., 1.,// bivector
    1., 1., 1., 1.,// trivector
    1.}; // pseudoscalar
  auto const expected = multivector{0., // scalar
    1., 2., 3., 4., // vector
    5., 6., 7., 8., 9., 10.,// bivector
    11., 12., 13., 14.,// trivector
    15.}; // pseudoscalar
  EXPECT_EQ(a - b, expected);
}

TEST(Multivector, ScalarAddition) {
  auto const a = multivector{1., // scalar
    2., 3., 4., 5., // vector
    6., 7., 8., 9., 10., 11.,// bivector
    12., 13., 14., 15.,// trivector
    16.}; // pseudoscalar

  auto const expected = multivector{2., // scalar 
    2., 3., 4., 5., // vector
    6., 7., 8., 9., 10., 11.,// bivector
    12., 13., 14., 15.,// trivector
    16.}; // pseudoscalar
  EXPECT_EQ(a + 1., expected);
}

TEST(Multivector, ScalarSubtraction){
  auto const a = multivector{1., // scalar
    2., 3., 4., 5., // vector
    6., 7., 8., 9., 10., 11.,// bivector
    12., 13., 14., 15.,// trivector
    16.}; // pseudoscalar

  auto const expected = multivector{0., // scalar 
    2., 3., 4., 5., // vector
    6., 7., 8., 9., 10., 11.,// bivector
    12., 13., 14., 15.,// trivector
    16.}; // pseudoscalar
  EXPECT_EQ(a - 1., expected);
}

TEST(Multivector, ScalarMultiplication){
  auto const a = multivector{1., // scalar
    2., 3., 4., 5., // vector
    6., 7., 8., 9., 10., 11.,// bivector
    12., 13., 14., 15.,// trivector
    16.}; // pseudoscalar

  auto const expected = multivector{2., // scalar 
    4., 6., 8., 10., // vector
    12., 14., 16., 18., 20., 22.,// bivector
    24., 26., 28., 30.,// trivector
    32.}; // pseudoscalar
  EXPECT_EQ(a * 2., expected);
}

namespace geometric_product{

struct scenario {
    std::string display;
    multivector lhs;
    multivector rhs;
    multivector expected;
};

std::ostream& operator<<(std::ostream& os, scenario const& s) {
    return os << s.display;
}

auto const scenarios = ::testing::Values(
    scenario{"1 * 1", one, one, one},
    scenario{"e0123 * e0123", e0123, e0123, zero});

struct GeometricProductFixture : public ::testing::TestWithParam<scenario> {};

TEST_P(GeometricProductFixture, GeometricProduct) {
    // GIVEN two multivectors 
    auto const [_, lhs, rhs, expected] = GetParam();

    // WHEN we take their geometric product
    auto const result = lhs * rhs; 

    // THEN it will be the product we expect
    EXPECT_EQ(result, expected);
}

// associates test and scenario as each can be reused
// values can also be automatically generated or combined
// http://google.github.io/googletest/reference/testing.html#INSTANTIATE_TEST_SUITE_P
INSTANTIATE_TEST_SUITE_P(GeometricProductGroup, GeometricProductFixture, scenarios);

}  // namespace test_transitions
} // namespace landy::pga
