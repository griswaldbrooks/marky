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

TEST(Multivector, ScalarProduct){
  EXPECT_EQ(one * one, one);
}

struct scenario {
    std::string display;
    multivector lhs;
    multivector rhs;
    multivector expected;
};

std::ostream& operator<<(std::ostream& os, scenario const& s) {
    return os << s.display;
}

auto const vector_product = ::testing::Values(
    scenario{"e0 * e0", e0, e0, zero},
    scenario{"e0 * e1", e0, e1, e01},
    scenario{"e0 * e2", e0, e2, e02},
    scenario{"e0 * e3", e0, e3, e03},
    scenario{"e1 * e0", e1, e0, -e01},
    scenario{"e1 * e1", e1, e1, one},
    scenario{"e1 * e2", e1, e2, e12},
    scenario{"e1 * e3", e1, e3, -e31},
    scenario{"e2 * e0", e2, e0, -e02},
    scenario{"e2 * e1", e2, e1, -e12},
    scenario{"e2 * e2", e2, e2, one},
    scenario{"e2 * e3", e2, e3, e23},
    scenario{"e3 * e0", e3, e0, -e03},
    scenario{"e3 * e1", e3, e1, e31},
    scenario{"e3 * e2", e3, e2, -e23},
    scenario{"e3 * e3", e3, e3, one});

auto const bivector_product = ::testing::Values(
    scenario{"e01 * e01", e01, e01, zero},
    scenario{"e01 * e02", e01, e02, zero},
    scenario{"e01 * e03", e01, e03, zero},
    scenario{"e01 * e12", e01, e12, e02},
    scenario{"e01 * e23", e01, e23, e0123},
    scenario{"e01 * e31", e01, e31, -e03},

    scenario{"e02 * e01", e02, e01, zero},
    scenario{"e02 * e02", e02, e02, zero},
    scenario{"e02 * e03", e02, e03, zero},
    scenario{"e02 * e12", e02, e12, -e01},
    scenario{"e02 * e23", e02, e23, e03},
    scenario{"e02 * e31", e02, e31, e0123},

    scenario{"e03 * e01", e03, e01, zero},
    scenario{"e03 * e02", e03, e02, zero},
    scenario{"e03 * e03", e03, e03, zero},
    scenario{"e03 * e12", e03, e12, e0123},
    scenario{"e03 * e23", e03, e23, -e02},
    scenario{"e03 * e31", e03, e31, e01},

    scenario{"e12 * e01", e12, e01, -e02},
    scenario{"e12 * e02", e12, e02, e01},
    scenario{"e12 * e03", e12, e03, e0123},
    scenario{"e12 * e12", e12, e12, -one},
    scenario{"e12 * e23", e12, e23, -e31},
    scenario{"e12 * e31", e12, e31, e23},

    scenario{"e23 * e01", e23, e01, e0123},
    scenario{"e23 * e02", e23, e02, -e03},
    scenario{"e23 * e03", e23, e03, e02},
    scenario{"e23 * e12", e23, e12, e31},
    scenario{"e23 * e23", e23, e23, -one},
    scenario{"e23 * e31", e23, e31, -e12},

    scenario{"e31 * e01", e31, e01, e03},
    scenario{"e31 * e02", e31, e02, e0123},
    scenario{"e31 * e03", e31, e03, -e01},
    scenario{"e31 * e12", e31, e12, -e23},
    scenario{"e31 * e23", e31, e23, e12},
    scenario{"e31 * e31", e31, e31, -one});

auto const bivector_vector_product = ::testing::Values(
    scenario{"e01 * e0", e01, e0, zero},
    scenario{"e01 * e1", e01, e1, e0},
    scenario{"e01 * e2", e01, e2, -e021},
    scenario{"e01 * e3", e01, e3, e013},

    scenario{"e02 * e0", e02, e0, zero},
    scenario{"e02 * e1", e02, e1, e021},
    scenario{"e02 * e2", e02, e2, e0},
    scenario{"e02 * e3", e02, e3, -e032},

    scenario{"e03 * e0", e03, e0, zero},
    scenario{"e03 * e1", e03, e1, -e013},
    scenario{"e03 * e2", e03, e2, e032},
    scenario{"e03 * e3", e03, e3, e0},

    scenario{"e12 * e0", e12, e0, -e021},
    scenario{"e12 * e1", e12, e1, -e2},
    scenario{"e12 * e2", e12, e2, e1},
    scenario{"e12 * e3", e12, e3, e123},

    scenario{"e31 * e0", e31, e0, -e013},
    scenario{"e31 * e1", e31, e1, e3},
    scenario{"e31 * e2", e31, e2, e123},
    scenario{"e31 * e3", e31, e3, -e1},

    scenario{"e23 * e0", e23, e0, -e032},
    scenario{"e23 * e1", e23, e1, e123},
    scenario{"e23 * e2", e23, e2, -e3},
    scenario{"e23 * e3", e23, e3, e2});

auto const vector_bivector_product = ::testing::Values(
    scenario{"e0 * e01", e0, e01, zero},
    scenario{"e0 * e02", e0, e02, zero},
    scenario{"e0 * e03", e0, e03, zero},
    scenario{"e0 * e12", e0, e12, -e021},
    scenario{"e0 * e31", e0, e31, -e013},
    scenario{"e0 * e23", e0, e23, -e032},

    scenario{"e1 * e01", e1, e01, -e0},
    scenario{"e1 * e02", e1, e02, e021},
    scenario{"e1 * e03", e1, e03, -e013},
    scenario{"e1 * e12", e1, e12, e2},
    scenario{"e1 * e31", e1, e31, -e3},
    scenario{"e1 * e23", e1, e23, e123},

    scenario{"e2 * e01", e2, e01, -e021},
    scenario{"e2 * e02", e2, e02, -e0},
    scenario{"e2 * e03", e2, e03, e032},
    scenario{"e2 * e12", e2, e12, -e1},
    scenario{"e2 * e31", e2, e31, e123},
    scenario{"e2 * e23", e2, e23, e3},

    scenario{"e3 * e01", e3, e01, e013},
    scenario{"e3 * e02", e3, e02, -e032},
    scenario{"e3 * e03", e3, e03, -e0},
    scenario{"e3 * e12", e3, e12, e123},
    scenario{"e3 * e31", e3, e31, e1},
    scenario{"e3 * e23", e3, e23, -e2});

auto const trivector_product = ::testing::Values(
    scenario("e021 * e021", e021, e021, zero),
    scenario("e021 * e013", e021, e013, zero),
    scenario("e021 * e032", e021, e032, zero),
    scenario("e021 * e123", e021, e123, e03),

    scenario("e013 * e021", e013, e021, zero),
    scenario("e013 * e013", e013, e013, zero),
    scenario("e013 * e032", e013, e032, zero),
    scenario("e013 * e123", e013, e123, e02),

    scenario("e032 * e021", e032, e021, zero),
    scenario("e032 * e013", e032, e013, zero),
    scenario("e032 * e032", e032, e032, zero),
    scenario("e032 * e123", e032, e123, e01),

    scenario("e123 * e021", e123, e021, -e03),
    scenario("e123 * e013", e123, e013, -e02),
    scenario("e123 * e032", e123, e032, -e01),
    scenario("e123 * e123", e123, e123, -one)
);

auto const trivector_vector_product = ::testing::Values(
    scenario("e021 * e0", e021, e0, zero),
    scenario("e021 * e1", e021, e1, e02),
    scenario("e021 * e2", e021, e2, -e01),
    scenario("e021 * e3", e021, e3, -e0123),

    scenario("e013 * e0", e013, e0, zero),
    scenario("e013 * e1", e013, e1, -e03),
    scenario("e013 * e2", e013, e2, -e0123),
    scenario("e013 * e3", e013, e3, e01),

    scenario("e032 * e0", e032, e0, zero),
    scenario("e032 * e1", e032, e1, -e0123),
    scenario("e032 * e2", e032, e2, e03),
    scenario("e032 * e3", e032, e3, -e02),

    scenario("e123 * e0", e123, e0, -e0123),
    scenario("e123 * e1", e123, e1, e23),
    scenario("e123 * e2", e123, e2, e31),
    scenario("e123 * e3", e123, e3, e12)
);

auto const vector_trivector_product = ::testing::Values(
  scenario("e0 * e021", e0, e021, zero),
  scenario("e0 * e013", e0, e013, zero),
  scenario("e0 * e032", e0, e032, zero),
  scenario("e0 * e123", e0, e123, e0123),

  scenario("e1 * e021", e1, e021, e02),
  scenario("e1 * e013", e1, e013, -e03),
  scenario("e1 * e032", e1, e032, e0123),
  scenario("e1 * e123", e1, e123, e23),

  scenario("e2 * e021", e2, e021, -e01),
  scenario("e2 * e013", e2, e013, e0123),
  scenario("e2 * e032", e2, e032, e03),
  scenario("e2 * e123", e2, e123, e31),

  scenario("e3 * e021", e3, e021, e0123),
  scenario("e3 * e013", e3, e013, e01),
  scenario("e3 * e032", e3, e032, -e02),
  scenario("e3 * e123", e3, e123, e12)
);

auto const trivector_bivector_product = ::testing::Values(
    scenario("e021 * e01", e021, e01, zero),
    scenario("e021 * e02", e021, e02, zero),
    scenario("e021 * e03", e021, e03, zero),
    scenario("e021 * e12", e021, e12, e0),
    scenario("e021 * e31", e021, e31, e032),
    scenario("e021 * e23", e021, e23, -e013),

    scenario("e013 * e01", e013, e01, zero),
    scenario("e013 * e02", e013, e02, zero),
    scenario("e013 * e03", e013, e03, zero),
    scenario("e013 * e12", e013, e12, -e032),
    scenario("e013 * e31", e013, e31, e0),
    scenario("e013 * e23", e013, e23, e021),

    scenario("e032 * e01", e032, e01, zero),
    scenario("e032 * e02", e032, e02, zero),
    scenario("e032 * e03", e032, e03, zero),
    scenario("e032 * e12", e032, e12, e013),
    scenario("e032 * e31", e032, e31, -e021),
    scenario("e032 * e23", e032, e23, e0),

    scenario("e123 * e01", e123, e01, e032),
    scenario("e123 * e02", e123, e02, e013),
    scenario("e123 * e03", e123, e03, e021),
    scenario("e123 * e12", e123, e12, -e3),
    scenario("e123 * e31", e123, e31, -e2),
    scenario("e123 * e23", e123, e23, -e1)
);

auto const bivector_trivector_product = ::testing::Values(
    scenario("e01 * e021", e01, e021, zero),
    scenario("e01 * e013", e01, e013, zero),
    scenario("e01 * e032", e01, e032, zero),
    scenario("e01 * e123", e01, e123, -e032),

    scenario("e02 * e021", e02, e021, zero),
    scenario("e02 * e013", e02, e013, zero),
    scenario("e02 * e032", e02, e032, zero),
    scenario("e02 * e123", e02, e123, -e013),

    scenario("e03 * e021", e03, e021, zero),
    scenario("e03 * e013", e03, e013, zero),
    scenario("e03 * e032", e03, e032, zero),
    scenario("e03 * e123", e03, e123, -e021),

    scenario("e12 * e021", e12, e021, e0),
    scenario("e12 * e013", e12, e013, e032),
    scenario("e12 * e032", e12, e032, -e013),
    scenario("e12 * e123", e12, e123, -e3),

    scenario("e31 * e021", e31, e021, -e032),
    scenario("e31 * e013", e31, e013, e0),
    scenario("e31 * e032", e31, e032, e021),
    scenario("e31 * e123", e31, e123, -e2),

    scenario("e23 * e021", e23, e021, e013),
    scenario("e23 * e013", e23, e013, -e021),
    scenario("e23 * e032", e23, e032, e0),
    scenario("e23 * e123", e23, e123, -e1)
);

auto const antiscalar_product = ::testing::Values(
    scenario("e0123 * e0", e0123, e0, zero),
    scenario("e0123 * e1", e0123, e1, -e032),
    scenario("e0123 * e2", e0123, e2, -e013),
    scenario("e0123 * e3", e0123, e3, -e021),
    scenario("e0 * e0123", e0, e0123, zero),
    scenario("e1 * e0123", e1, e0123, e032),
    scenario("e2 * e0123", e2, e0123, e013),
    scenario("e3 * e0123", e3, e0123, e021),

    scenario("e0123 * e01", e0123, e01, zero),
    scenario("e0123 * e02", e0123, e02, zero),
    scenario("e0123 * e03", e0123, e03, zero),
    scenario("e0123 * e12", e0123, e12, -e03),
    scenario("e0123 * e31", e0123, e31, -e02),
    scenario("e0123 * e23", e0123, e23, -e01),
    scenario("e01 * e0123", e01, e0123, zero),
    scenario("e02 * e0123", e02, e0123, zero),
    scenario("e03 * e0123", e03, e0123, zero),
    scenario("e12 * e0123", e12, e0123, -e03),
    scenario("e31 * e0123", e31, e0123, -e02),
    scenario("e23 * e0123", e23, e0123, -e01),

    scenario("e0123 * e021", e0123, e021, zero),
    scenario("e0123 * e013", e0123, e013, zero),
    scenario("e0123 * e032", e0123, e032, zero),
    scenario("e0123 * e123", e0123, e123, -e0),
    scenario("e021 * e0123", e021, e0123, zero),
    scenario("e013 * e0123", e013, e0123, zero),
    scenario("e032 * e0123", e032, e0123, zero),
    scenario("e123 * e0123", e123, e0123, e0),

    scenario("e0123 * e0123", e0123, e0123, zero)
);

struct GeometricProductFixture : public ::testing::TestWithParam<scenario> {};

TEST_P(GeometricProductFixture, GeometricProduct) {
    // GIVEN two multivectors
    auto const [_, lhs, rhs, expected] = GetParam();

    // WHEN we take their geometric product
    auto const result = lhs * rhs;

    // THEN it will be the product we expect
    EXPECT_EQ(result, expected);
}

INSTANTIATE_TEST_SUITE_P(GeometricProductVector, GeometricProductFixture, vector_product);
INSTANTIATE_TEST_SUITE_P(GeometricProductBivector, GeometricProductFixture, bivector_product);
INSTANTIATE_TEST_SUITE_P(GeometricProductBivectorVector, GeometricProductFixture, bivector_vector_product);
INSTANTIATE_TEST_SUITE_P(GeometricProductVectorBivector, GeometricProductFixture, vector_bivector_product);
INSTANTIATE_TEST_SUITE_P(GeometricProductTrivector, GeometricProductFixture, trivector_product);
INSTANTIATE_TEST_SUITE_P(GeometricProductTrivectorVector, GeometricProductFixture, trivector_vector_product);
INSTANTIATE_TEST_SUITE_P(GeometricProductVectorTrivector, GeometricProductFixture, vector_trivector_product);
INSTANTIATE_TEST_SUITE_P(GeometricProductTrivectorBivector, GeometricProductFixture, trivector_bivector_product);
INSTANTIATE_TEST_SUITE_P(GeometricProductBivectorTrivector, GeometricProductFixture, bivector_trivector_product);
INSTANTIATE_TEST_SUITE_P(GeometricProductAntiscalar, GeometricProductFixture, antiscalar_product);

}  // namespace geometric_product
} // namespace landy::pga
