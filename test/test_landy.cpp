#include <gtest/gtest.h>
#include <span>
#include "landy/landy.hpp"

TEST(ComputePose, NoOp) {
  using landy::geometry::pose_t;
  using landy::geometry::compute_pose;
  auto const observations = std::array{//
    pose_t{0, 0, 0}, //
    pose_t{0, 0, 0}, //
    pose_t{0, 0, 0}, //
  };
  auto const constellation = std::array{
    pose_t{0, 0, 0}, //
    pose_t{0, 0, 0}, //
    pose_t{0, 0, 0}, //
  };
  auto const pose = compute_pose(std::span{observations}, std::span{constellation});
  EXPECT_TRUE(pose.x == 0.);
}
