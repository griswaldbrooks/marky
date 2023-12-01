#include <gtest/gtest.h>
#include <span>
#include "landy/landy.hpp"

TEST(ComputePose, NoOp) {
  using landy::geometry::point_t;
  using landy::geometry::compute_pose;
  auto const observations = std::array{//
    point_t{0, 0, 0}, //
    point_t{0, 0, 0}, //
    point_t{0, 0, 0}, //
  };
  auto const constellation = std::array{
    point_t{0, 0, 0}, //
    point_t{0, 0, 0}, //
    point_t{0, 0, 0}, //
  };
  auto const pose = compute_pose(std::span{observations}, std::span{constellation});
  EXPECT_TRUE(pose.x == 0.);
}

TEST(Dot, Dot) {
  using landy::geometry::point_t;
  using landy::geometry::dot;
  auto const a = point_t{1, 2, 3};
  auto const b = point_t{4, 5, 6};
  EXPECT_EQ(dot(a, b), 32);
}

TEST(Cross, Cross) {
  using landy::geometry::point_t;
  using landy::geometry::cross;
  auto const a = point_t{1, 2, 3};
  auto const b = point_t{4, 5, 6};
  auto const expected = point_t{-3, 6, -3};
  EXPECT_EQ(cross(a, b), expected);
}
