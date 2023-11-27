#pragma once
#include <span>
namespace landy::geometry {

struct pose_t { double x, y, z; };

/**
* @brief Given a set of observed landmarks, produce the observer's pose
* @param observations of the landmarks in the scene
* @param constellation of known landmarks in environment
* @note constellation will become it's own type, but observations will probably stay a list of poses.
*       what i mean to say is that any iterable sequence i will take and i don't care if it's a 
*       vector, array, or generator.
* @returns pose of observer
* @note something about this is wrong, since there is a required number of observations for this to be solved
*       exactly. too few and it's ambiguous, too many and a fitting function will be required.
*       this is meant to be the precursor to whatever estimator tracks the pose and i'm sure there's a normal
*       name for this.
*/
template <size_t N, size_t M>
[[nodiscard]] pose_t compute_pose([[maybe_unused]] std::span<pose_t const, N> observations, [[maybe_unused]] std::span<pose_t const, M> constellation) {
  return {};
}

} // geometry
