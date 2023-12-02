#include <gtest/gtest.h>
#include <span>
#include "landy/landy.hpp"

namespace landy::geometry {


// TEST(ComputePose, NoOp) {
//   auto const observations = std::array{//
//     point_t{0, 0, 0}, //
//     point_t{0, 0, 0}, //
//     point_t{0, 0, 0}, //
//   };
//   auto const constellation = std::array{
//     point_t{0, 0, 0}, //
//     point_t{0, 0, 0}, //
//     point_t{0, 0, 0}, //
//   };
//   auto const pose = compute_pose(std::span{observations}, std::span{constellation});
//   EXPECT_TRUE(pose.x == 0.);
// }

TEST(Dot, Dot) {
  auto const a = point_t{1, 2, 3};
  auto const b = point_t{4, 5, 6};
  EXPECT_EQ(dot(a, b), 32);
}

TEST(Cross, Cross) {
  auto const a = point_t{1, 2, 3};
  auto const b = point_t{4, 5, 6};
  auto const expected = point_t{-3, 6, -3};
  EXPECT_EQ(cross(a, b), expected);
}

TEST(IsNear, IsNear) {
  auto const a = point_t{1, 2, 3};
  auto const b = point_t{1, 2, 3};
  EXPECT_TRUE(is_near(a, b));
}

TEST(IsNear, IsNearTolerance) {
  auto const a = point_t{1, 2, 3};
  auto const b = point_t{1, 2, 3.0000001};
  EXPECT_TRUE(is_near(a, b));
}

TEST(IsNear, IsNearToleranceFalse) {
  auto const a = point_t{1, 2, 3};
  auto const b = point_t{1, 2, 3.0000001};
  EXPECT_FALSE(is_near(a, b, 1e-8));
}

TEST(Norm, Norm) {
  auto const a = point_t{1, 2, 3};
  EXPECT_TRUE(is_near(norm(a), std::sqrt(14)));
}

TEST(IsSkew, IsSkew) {
  auto const a = line_t{point_t{0, 0, 0}, point_t{1, 0, 0}};
  auto const b = line_t{point_t{0, 0, 1}, point_t{0, 1, 0}};
  EXPECT_TRUE(is_skew(a, b));
}

TEST(IsSkew, SameLine) {
  auto const a = line_t{point_t{0, 0, 0}, point_t{1, 0, 0}};
  auto const b = line_t{point_t{0, 0, 0}, point_t{1, 0, 0}};
  EXPECT_FALSE(is_skew(a, b));
}

TEST(IsSkew, Parallel) {
  auto const a = line_t{point_t{0, 0, 0}, point_t{1, 0, 0}};
  auto const b = line_t{point_t{0, 0, 1}, point_t{1, 0, 1}};
  EXPECT_FALSE(is_skew(a, b));
}

TEST(IsSkew, Intersect) {
  auto const a = line_t{point_t{0, 0, 0}, point_t{1, 0, 0}};
  auto const b = line_t{point_t{0, 0, 0}, point_t{0, 0, 1}};
  EXPECT_FALSE(is_skew(a, b));
}

TEST(MidPoint, Intersect) {
  auto const a = line_t{point_t{0, 0, 0}, point_t{1, 0, 0}};
  auto const b = line_t{point_t{0, 0, 0}, point_t{0, 0, 1}};
  auto const mp = midpoint(a, b).value();
  auto const expected = point_t{0, 0, 0};
  EXPECT_TRUE(is_near(mp, expected));
}

TEST(MidPoint, Skew) {
  auto const a = line_t{point_t{0, 0, 0}, point_t{1, 0, 0}};
  auto const b = line_t{point_t{0, 0, 1}, point_t{0, 1, 0}};
  auto const mp = midpoint(a, b).value();
  auto const expected = point_t{0., 0., 0.5};
  EXPECT_TRUE(is_near(mp, expected));
}

TEST(MidPoint, Parallel) {
  auto const a = line_t{point_t{0, 0, 0}, point_t{1, 0, 0}};
  auto const b = line_t{point_t{0, 0, 1}, point_t{1, 0, 1}};
  auto const mp = midpoint(a, b);
  EXPECT_FALSE(mp.has_value());
} 
} // namespace landy::geometry
