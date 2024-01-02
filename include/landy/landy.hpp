#pragma once
#include <optional>
#include <cmath>
#include <utility>
#include <numbers>
#include "mp-units/systems/si/si.h"
#include "mp-units/systems/isq/isq.h"
#include "mp-units/math.h"
#include <klein/klein.hpp>

using namespace mp_units;
using namespace::si::unit_symbols;
using namespace std::numbers;

kln::point testklein() {
  // Create a rotor representing a pi/2 rotation about the z-axis
// Normalization is done automatically
  auto const r = kln::rotor{pi_v<float> * 0.5f, 0.f, 0.f, 1.f};

// Create a translator that represents a translation of 1 unit
// in the yz-direction. Normalization is done automatically.
  auto const t = kln::translator{1.f, 0.f, 1.f, 1.f};

// Create a motor that combines the action of the rotation and
// translation above.
  auto const m = r * t;

// Construct a point at position (1, 0, 0)
  auto const p1 = kln::point{1, 0, 0};

// Apply the motor to the point. This is equivalent to the conjugation
// operator m * p1 * ~m where * is the geometric product and ~ is the
// reverse operation.
  auto const p2 = m(p1);
  return p2;
}

namespace landy::geometry {

template <typename T>
[[nodiscard]] bool is_near(T const& lhs, T const& rhs, T const& tolerance = 1e-6) {
  return std::abs(lhs - rhs) < tolerance;
}

[[nodiscard]] constexpr auto norm(mp_units::QuantityOf<isq::position_vector> auto const& lhs) {
  using mp_units::hypot;
  return hypot(lhs.x, lhs.y, lhs.z);
}

// point_t operator+(point_t const& lhs, point_t const& rhs) {
//   return {lhs.x + rhs.x, lhs.y + rhs.y, lhs.z + rhs.z};
// }
//
// point_t operator-(point_t const& lhs, point_t const& rhs) {
//   return {lhs.x - rhs.x, lhs.y - rhs.y, lhs.z - rhs.z};
// }
//
// point_t operator*(point_t const& lhs, double const& rhs) {
//   return {lhs.x * rhs, lhs.y * rhs, lhs.z * rhs};
// }
//
// point_t operator*(double const& lhs, point_t const& rhs) {
//   return rhs * lhs;
// }
//
// point_t operator/(point_t const& lhs, double const& rhs) {
//   return {lhs.x / rhs, lhs.y / rhs, lhs.z / rhs};
// }
//
// [[nodiscard]] bool is_near(point_t const& lhs, point_t const& rhs, double tolerance = 1e-6) {
//   return norm(lhs - rhs) < tolerance;
// }
//
// bool operator==(point_t const& lhs, point_t const& rhs) {
//   return is_near(lhs, rhs);
// }
//
// struct observation_t {
//   point_t point;
//   // bearing or range? can it be both or is this mutally exclusive?
//   // std::variant, tenplate, tuple?
// };
//
// struct line_t {
//   point_t p;
//   point_t q;
// };
//
// [[nodiscard]] point_t cross(point_t const& p, point_t const& q) {
//   return {p.y * q.z - p.z * q.y, p.z * q.x - p.x * q.z, p.x * q.y - p.y * q.x};
// }
//
// [[nodiscard]] double dot(point_t const& p, point_t const& q) {
//   return p.x * q.x + p.y * q.y + p.z * q.z;
// }
//
// [[nodiscard]] bool is_skew(line_t const& a, line_t const& b) {
//   // https://en.wikipedia.org/wiki/Skew_lines#Nearest_points
//   // https://mathforyou.net/en/online/vectors/volume/tetrahedron/
//   // Compute the triple product to get the volume of the tetrahedron
//   auto const volume = dot(a.p - b.p, cross(a.q - a.p, b.q - b.p));
//   return !is_near(volume, 0.);
// }
//
// [[nodiscard]] bool is_parallel(line_t const& a, line_t const& b){
//   return is_near(norm(cross(a.q - a.p, b.q - b.p)), 0.);
// }
//
// // 3d lines that are neither parallel nor intersecting are called skew lines
// // https://en.wikipedia.org/wiki/Skew_lines#Nearest_points
// // change to std::expected
// std::optional<point_t> midpoint(line_t const& a, line_t const& b) {
//   if (is_parallel(a, b)) {
//     return std::nullopt;
//   }
//   auto const n = cross(a.q, b.q); // normal to both lines
//   auto const n1 = cross(a.q, n); // normal to line a and n
//   auto const n2 = cross(b.q, n); // normal to line b and n
//
//   auto const c1 = a.p + a.q * dot(b.p - a.p, n2) / dot(a.q, n2);
//   auto const c2 = b.p + b.q * dot(a.p - b.p, n1) / dot(b.q, n1);
//
//   return 0.5 * (c1 + c2);
// }
} // geometry
