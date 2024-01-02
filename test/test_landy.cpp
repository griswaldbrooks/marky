#include <gtest/gtest.h>
#include <span>
#include "landy/landy.hpp"
#include "mp-units/systems/si/si.h"
#include "mp-units/systems/isq/isq.h"
#include <mp-units/customization_points.h>

using namespace mp_units;
using namespace::si::unit_symbols;

using point_t = std::array<double, 3>;

// template<>
// inline constexpr bool mp_units::treat_as_floating_point<point_t> = mp_units::treat_as_floating_point<double>;
//
// template<>
// inline constexpr bool mp_units::is_vector<point_t> = true;
//
//
// template<Reference R>
// Quantity auto operator*(point_t rep, R) {
//   return make_quantity<R>(rep);
// }
//
namespace landy::geometry {

// TEST(Norm, Norm){
//   auto const point = point_t{0., 1.1, 2.2} * isq::position_vector[m];
//   EXPECT_EQ(norm(point), 1.);
// }

// TEST(Dot, Dot) {
//   auto const a = point_t{1, 2, 3};
//   auto const b = point_t{4, 5, 6};
//   EXPECT_EQ(dot(a, b), 32);
// }

// TEST(Cross, Cross) {
//   auto const a = point_t{1, 2, 3};
//   auto const b = point_t{4, 5, 6};
//   auto const expected = point_t{-3, 6, -3};
//   EXPECT_EQ(cross(a, b), expected);
// }
//
// TEST(IsNear, IsNear) {
//   auto const a = point_t{1, 2, 3};
//   auto const b = point_t{1, 2, 3};
//   EXPECT_TRUE(is_near(a, b));
// }
//
// TEST(IsNear, IsNearTolerance) {
//   auto const a = point_t{1, 2, 3};
//   auto const b = point_t{1, 2, 3.0000001};
//   EXPECT_TRUE(is_near(a, b));
// }
//
// TEST(IsNear, IsNearToleranceFalse) {
//   auto const a = point_t{1, 2, 3};
//   auto const b = point_t{1, 2, 3.0000001};
//   EXPECT_FALSE(is_near(a, b, 1e-8));
// }
//
// TEST(Norm, Norm) {
//   auto const a = point_t{1, 2, 3};
//   EXPECT_TRUE(is_near(norm(a), std::sqrt(14)));
// }
//
// TEST(IsSkew, IsSkew) {
//   auto const a = line_t{point_t{0, 0, 0}, point_t{1, 0, 0}};
//   auto const b = line_t{point_t{0, 0, 1}, point_t{0, 1, 0}};
//   EXPECT_TRUE(is_skew(a, b));
// }
//
// TEST(IsSkew, SameLine) {
//   auto const a = line_t{point_t{0, 0, 0}, point_t{1, 0, 0}};
//   auto const b = line_t{point_t{0, 0, 0}, point_t{1, 0, 0}};
//   EXPECT_FALSE(is_skew(a, b));
// }
//
// TEST(IsSkew, Parallel) {
//   auto const a = line_t{point_t{0, 0, 0}, point_t{1, 0, 0}};
//   auto const b = line_t{point_t{0, 0, 1}, point_t{1, 0, 1}};
//   EXPECT_FALSE(is_skew(a, b));
// }
//
// TEST(IsSkew, Intersect) {
//   auto const a = line_t{point_t{0, 0, 0}, point_t{1, 0, 0}};
//   auto const b = line_t{point_t{0, 0, 0}, point_t{0, 0, 1}};
//   EXPECT_FALSE(is_skew(a, b));
// }
//
// TEST(MidPoint, Intersect) {
//   auto const a = line_t{point_t{0, 0, 0}, point_t{1, 0, 0}};
//   auto const b = line_t{point_t{0, 0, 0}, point_t{0, 0, 1}};
//   auto const mp = midpoint(a, b).value();
//   auto const expected = point_t{0, 0, 0};
//   EXPECT_TRUE(is_near(mp, expected));
// }
//
// TEST(MidPoint, Skew) {
//   auto const a = line_t{point_t{0, 0, 0}, point_t{1, 0, 0}};
//   auto const b = line_t{point_t{0, 0, 1}, point_t{0, 1, 0}};
//   auto const mp = midpoint(a, b).value();
//   auto const expected = point_t{0., 0., 0.5};
//   EXPECT_TRUE(is_near(mp, expected));
// }
//
// TEST(MidPoint, Parallel) {
//   auto const a = line_t{point_t{0, 0, 0}, point_t{1, 0, 0}};
//   auto const b = line_t{point_t{0, 0, 1}, point_t{1, 0, 1}};
//   auto const mp = midpoint(a, b);
//   EXPECT_FALSE(mp.has_value());
// }
} // namespace landy::geometry
