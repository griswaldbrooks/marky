#pragma once
#include <optional>
#include <cmath>
namespace landy::geometry {

template <typename T>
[[nodiscard]] bool is_near(T const& lhs, T const& rhs, T const& tolerance = 1e-6) {
  return std::abs(lhs - rhs) < tolerance;
}

struct point_t {
  double x = 0.;
  double y = 0.;
  double z = 0.;
};

[[nodiscard]] double norm(point_t const& lhs) {
  return std::hypot(lhs.x, lhs.y, lhs.z);
}

point_t operator+(point_t const& lhs, point_t const& rhs) {
  return {lhs.x + rhs.x, lhs.y + rhs.y, lhs.z + rhs.z};
}

point_t operator-(point_t const& lhs, point_t const& rhs) {
  return {lhs.x - rhs.x, lhs.y - rhs.y, lhs.z - rhs.z};
}

point_t operator*(point_t const& lhs, double const& rhs) {
  return {lhs.x * rhs, lhs.y * rhs, lhs.z * rhs};
}

point_t operator*(double const& lhs, point_t const& rhs) {
  return rhs * lhs;
}

point_t operator/(point_t const& lhs, double const& rhs) {
  return {lhs.x / rhs, lhs.y / rhs, lhs.z / rhs};
}

[[nodiscard]] bool is_near(point_t const& lhs, point_t const& rhs, double tolerance = 1e-6) {
  return norm(lhs - rhs) < tolerance;
}

bool operator==(point_t const& lhs, point_t const& rhs) {
  return is_near(lhs, rhs);
}

struct observation_t {
  point_t point;
  // bearing or range? can it be both or is this mutally exclusive?
  // std::variant, tenplate, tuple?
};

struct line_t {
  point_t p;
  point_t q;
};

[[nodiscard]] point_t cross(point_t const& p, point_t const& q) {
  return {p.y * q.z - p.z * q.y, p.z * q.x - p.x * q.z, p.x * q.y - p.y * q.x};
}

[[nodiscard]] double dot(point_t const& p, point_t const& q) {
  return p.x * q.x + p.y * q.y + p.z * q.z;
}

[[nodiscard]] bool is_skew(line_t const& a, line_t const& b) {
  // https://en.wikipedia.org/wiki/Skew_lines#Nearest_points
  // https://mathforyou.net/en/online/vectors/volume/tetrahedron/
  // Compute the triple product to get the volume of the tetrahedron
  auto const volume = dot(a.p - b.p, cross(a.q - a.p, b.q - b.p));
  return !is_near(volume, 0.);
}

[[nodiscard]] bool is_parallel(line_t const& a, line_t const& b){
  return is_near(norm(cross(a.q - a.p, b.q - b.p)), 0.);
}

// 3d lines that are neither parallel nor intersecting are called skew lines
// https://en.wikipedia.org/wiki/Skew_lines#Nearest_points
// change to std::expected
std::optional<point_t> midpoint(line_t const& a, line_t const& b) {
  if (is_parallel(a, b)) {
    return std::nullopt;
  }
  auto const n = cross(a.q, b.q); // normal to both lines
  auto const n1 = cross(a.q, n); // normal to line a and n
  auto const n2 = cross(b.q, n); // normal to line b and n

  auto const c1 = a.p + a.q * dot(b.p - a.p, n2) / dot(a.q, n2);
  auto const c2 = b.p + b.q * dot(a.p - b.p, n1) / dot(b.q, n1);

  return 0.5 * (c1 + c2);
}

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
[[nodiscard]] point_t compute_pose([[maybe_unused]] std::span<point_t const, N> observations, [[maybe_unused]] std::span<point_t const, M> constellation) {
  // given a known constellation, just knowing which landmarks are in view isn't enough.
  // for example, if the landmarks allow for unique identification, and nothing else,
  // then just having a subset of the constellation tells you nothing.
  // you need to know
  // - distance to landmarks
  // - bearing to landmarks
  // - field of view or something else about the sensor telling you
  //   something about the robot by being able to observe these landmarks
  //
  // if you know the range to a landmark, then your robot exists on a circle relative to that point
  // if you know the yaw and pitch (azimuth/elevation?) then you are on a line passing through that point
  //
  // two range observations reduces your pose set to two points (circles intersect)
  // two bearing observations reduces you pose set to two points (i think your pose mirrors over a plane?)
  //
  // so, three observations uniquely identify your pose.
  // since these measurements will be noisy and therefore three points will likely not have a solution, will
  // an optimization algorithm prefer to work with the raw lines and circles, or will it be ok with points?
  // feels like working on the lines and points is preferable
  // observed lines might not even intersect,
  //
  // need a name for those sets (the lines or circles, etc) for when there are less than 3 observations
  // actually in 3d these are spheres for ranges, so two observations give a ellipsoid, not two points
  // i bet you can give the ground plane constrain for known 2d cases
  //
  //
  //
  return {};
}

} // geometry
