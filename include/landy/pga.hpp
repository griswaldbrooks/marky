#pragma once
#include <array>
#include <cmath>
#include <numbers>

namespace pga {

using namespace std::numbers;

struct multivector {
 template <std::size_t Index>
    constexpr void set(double value) {
    static_assert(Index < 16, "Multivectors only have 16 elements");
        std::get<Index>(elements_) = value;
    } 
  constexpr double& operator[](std::size_t idx) { return elements_[idx]; }
  constexpr double operator[](std::size_t idx) const { return elements_[idx]; }
  // implements const iterator for the elements 
  constexpr auto cbegin() const { return elements_.cbegin(); }
  constexpr auto cend() const { return elements_.cend(); }
  constexpr auto begin() { return elements_.begin(); }
  constexpr auto end() { return elements_.end(); }

private:
  std::array<double, 16> elements_ = {};
};

template <std::size_t Index>
constexpr multivector make_vector(double value) {
  multivector result;
  result.set<Index>(value);
  return result;
}

constexpr bool operator==(multivector const& a, multivector const& b) {
  return std::equal(a.cbegin(), a.cend(), b.cbegin(), b.cend());
}
//***********************
// multivector.Reverse : res = ~a
// Reverse the order of the basis blades.
//***********************
constexpr multivector operator~(multivector const& a) {
  multivector res;
  res[0] = a[0];
  res[1] = a[1];
  res[2] = a[2];
  res[3] = a[3];
  res[4] = a[4];
  res[5] = -a[5];
  res[6] = -a[6];
  res[7] = -a[7];
  res[8] = -a[8];
  res[9] = -a[9];
  res[10] = -a[10];
  res[11] = -a[11];
  res[12] = -a[12];
  res[13] = -a[13];
  res[14] = -a[14];
  res[15] = a[15];
  return res;
};

//***********************
// multivector.Dual : res = !a
// Poincare duality operator.
//***********************
constexpr multivector operator!(multivector const&a) {
  multivector res;
  res[0] = a[15];
  res[1] = a[14];
  res[2] = a[13];
  res[3] = a[12];
  res[4] = a[11];
  res[5] = a[10];
  res[6] = a[9];
  res[7] = a[8];
  res[8] = a[7];
  res[9] = a[6];
  res[10] = a[5];
  res[11] = a[4];
  res[12] = a[3];
  res[13] = a[2];
  res[14] = a[1];
  res[15] = a[0];
  return res;
};

//***********************
// multivector.Conjugate : res = a.Conjugate()
// Clifford Conjugation
//***********************
constexpr multivector conjugate(multivector const& a) {
  multivector res;
  res[0] = a[0];
  res[1] = -a[1];
  res[2] = -a[2];
  res[3] = -a[3];
  res[4] = -a[4];
  res[5] = -a[5];
  res[6] = -a[6];
  res[7] = -a[7];
  res[8] = -a[8];
  res[9] = -a[9];
  res[10] = -a[10];
  res[11] = a[11];
  res[12] = a[12];
  res[13] = a[13];
  res[14] = a[14];
  res[15] = a[15];
  return res;
};

//***********************
// multivector.Involute : res = a.Involute()
// Main involution
//***********************
constexpr multivector involute(multivector const& a) {
  multivector res;
  res[0] = a[0];
  res[1] = -a[1];
  res[2] = -a[2];
  res[3] = -a[3];
  res[4] = -a[4];
  res[5] = a[5];
  res[6] = a[6];
  res[7] = a[7];
  res[8] = a[8];
  res[9] = a[9];
  res[10] = a[10];
  res[11] = -a[11];
  res[12] = -a[12];
  res[13] = -a[13];
  res[14] = -a[14];
  res[15] = a[15];
  return res;
};

//***********************
// multivector.Mul : res = a * b
// The geometric product.
//***********************
constexpr multivector operator*(multivector const& a, multivector const& b) {
  multivector res;
  res[0] = b[0] * a[0] + b[2] * a[2] + b[3] * a[3] + b[4] * a[4] - b[8] * a[8] -
           b[9] * a[9] - b[10] * a[10] - b[14] * a[14];
  res[1] = b[1] * a[0] + b[0] * a[1] - b[5] * a[2] - b[6] * a[3] - b[7] * a[4] +
           b[2] * a[5] + b[3] * a[6] + b[4] * a[7] + b[11] * a[8] +
           b[12] * a[9] + b[13] * a[10] + b[8] * a[11] + b[9] * a[12] +
           b[10] * a[13] + b[15] * a[14] - b[14] * a[15];
  res[2] = b[2] * a[0] + b[0] * a[2] - b[8] * a[3] + b[9] * a[4] + b[3] * a[8] -
           b[4] * a[9] - b[14] * a[10] - b[10] * a[14];
  res[3] = b[3] * a[0] + b[8] * a[2] + b[0] * a[3] - b[10] * a[4] -
           b[2] * a[8] - b[14] * a[9] + b[4] * a[10] - b[9] * a[14];
  res[4] = b[4] * a[0] - b[9] * a[2] + b[10] * a[3] + b[0] * a[4] -
           b[14] * a[8] + b[2] * a[9] - b[3] * a[10] - b[8] * a[14];
  res[5] = b[5] * a[0] + b[2] * a[1] - b[1] * a[2] - b[11] * a[3] +
           b[12] * a[4] + b[0] * a[5] - b[8] * a[6] + b[9] * a[7] +
           b[6] * a[8] - b[7] * a[9] - b[15] * a[10] - b[3] * a[11] +
           b[4] * a[12] + b[14] * a[13] - b[13] * a[14] - b[10] * a[15];
  res[6] = b[6] * a[0] + b[3] * a[1] + b[11] * a[2] - b[1] * a[3] -
           b[13] * a[4] + b[8] * a[5] + b[0] * a[6] - b[10] * a[7] -
           b[5] * a[8] - b[15] * a[9] + b[7] * a[10] + b[2] * a[11] +
           b[14] * a[12] - b[4] * a[13] - b[12] * a[14] - b[9] * a[15];
  res[7] = b[7] * a[0] + b[4] * a[1] - b[12] * a[2] + b[13] * a[3] -
           b[1] * a[4] - b[9] * a[5] + b[10] * a[6] + b[0] * a[7] -
           b[15] * a[8] + b[5] * a[9] - b[6] * a[10] + b[14] * a[11] -
           b[2] * a[12] + b[3] * a[13] - b[11] * a[14] - b[8] * a[15];
  res[8] = b[8] * a[0] + b[3] * a[2] - b[2] * a[3] + b[14] * a[4] +
           b[0] * a[8] + b[10] * a[9] - b[9] * a[10] + b[4] * a[14];
  res[9] = b[9] * a[0] - b[4] * a[2] + b[14] * a[3] + b[2] * a[4] -
           b[10] * a[8] + b[0] * a[9] + b[8] * a[10] + b[3] * a[14];
  res[10] = b[10] * a[0] + b[14] * a[2] + b[4] * a[3] - b[3] * a[4] +
            b[9] * a[8] - b[8] * a[9] + b[0] * a[10] + b[2] * a[14];
  res[11] = b[11] * a[0] - b[8] * a[1] + b[6] * a[2] - b[5] * a[3] +
            b[15] * a[4] - b[3] * a[5] + b[2] * a[6] - b[14] * a[7] -
            b[1] * a[8] + b[13] * a[9] - b[12] * a[10] + b[0] * a[11] +
            b[10] * a[12] - b[9] * a[13] + b[7] * a[14] - b[4] * a[15];
  res[12] = b[12] * a[0] - b[9] * a[1] - b[7] * a[2] + b[15] * a[3] +
            b[5] * a[4] + b[4] * a[5] - b[14] * a[6] - b[2] * a[7] -
            b[13] * a[8] - b[1] * a[9] + b[11] * a[10] - b[10] * a[11] +
            b[0] * a[12] + b[8] * a[13] + b[6] * a[14] - b[3] * a[15];
  res[13] = b[13] * a[0] - b[10] * a[1] + b[15] * a[2] + b[7] * a[3] -
            b[6] * a[4] - b[14] * a[5] - b[4] * a[6] + b[3] * a[7] +
            b[12] * a[8] - b[11] * a[9] - b[1] * a[10] + b[9] * a[11] -
            b[8] * a[12] + b[0] * a[13] + b[5] * a[14] - b[2] * a[15];
  res[14] = b[14] * a[0] + b[10] * a[2] + b[9] * a[3] + b[8] * a[4] +
            b[4] * a[8] + b[3] * a[9] + b[2] * a[10] + b[0] * a[14];
  res[15] = b[15] * a[0] + b[14] * a[1] + b[13] * a[2] + b[12] * a[3] +
            b[11] * a[4] + b[10] * a[5] + b[9] * a[6] + b[8] * a[7] +
            b[7] * a[8] + b[6] * a[9] + b[5] * a[10] - b[4] * a[11] -
            b[3] * a[12] - b[2] * a[13] - b[1] * a[14] + b[0] * a[15];
  return res;
};

//***********************
// multivector.Wedge : res = a ^ b
// The outer product. (MEET)
//***********************
 constexpr multivector operator^(multivector const& a, multivector const& b) {
  multivector res;
  res[0] = b[0] * a[0];
  res[1] = b[1] * a[0] + b[0] * a[1];
  res[2] = b[2] * a[0] + b[0] * a[2];
  res[3] = b[3] * a[0] + b[0] * a[3];
  res[4] = b[4] * a[0] + b[0] * a[4];
  res[5] = b[5] * a[0] + b[2] * a[1] - b[1] * a[2] + b[0] * a[5];
  res[6] = b[6] * a[0] + b[3] * a[1] - b[1] * a[3] + b[0] * a[6];
  res[7] = b[7] * a[0] + b[4] * a[1] - b[1] * a[4] + b[0] * a[7];
  res[8] = b[8] * a[0] + b[3] * a[2] - b[2] * a[3] + b[0] * a[8];
  res[9] = b[9] * a[0] - b[4] * a[2] + b[2] * a[4] + b[0] * a[9];
  res[10] = b[10] * a[0] + b[4] * a[3] - b[3] * a[4] + b[0] * a[10];
  res[11] = b[11] * a[0] - b[8] * a[1] + b[6] * a[2] - b[5] * a[3] -
            b[3] * a[5] + b[2] * a[6] - b[1] * a[8] + b[0] * a[11];
  res[12] = b[12] * a[0] - b[9] * a[1] - b[7] * a[2] + b[5] * a[4] +
            b[4] * a[5] - b[2] * a[7] - b[1] * a[9] + b[0] * a[12];
  res[13] = b[13] * a[0] - b[10] * a[1] + b[7] * a[3] - b[6] * a[4] -
            b[4] * a[6] + b[3] * a[7] - b[1] * a[10] + b[0] * a[13];
  res[14] = b[14] * a[0] + b[10] * a[2] + b[9] * a[3] + b[8] * a[4] +
            b[4] * a[8] + b[3] * a[9] + b[2] * a[10] + b[0] * a[14];
  res[15] = b[15] * a[0] + b[14] * a[1] + b[13] * a[2] + b[12] * a[3] +
            b[11] * a[4] + b[10] * a[5] + b[9] * a[6] + b[8] * a[7] +
            b[7] * a[8] + b[6] * a[9] + b[5] * a[10] - b[4] * a[11] -
            b[3] * a[12] - b[2] * a[13] - b[1] * a[14] + b[0] * a[15];
  return res;
};

//***********************
// multivector.Vee : res = a & b
// The regressive product. (JOIN)
//***********************
 constexpr multivector operator&(multivector const& a, multivector const& b) {
  multivector res;
  res[15] = 1 * (a[15] * b[15]);
  res[14] = -1 * (a[14] * -1 * b[15] + a[15] * b[14] * -1);
  res[13] = -1 * (a[13] * -1 * b[15] + a[15] * b[13] * -1);
  res[12] = -1 * (a[12] * -1 * b[15] + a[15] * b[12] * -1);
  res[11] = -1 * (a[11] * -1 * b[15] + a[15] * b[11] * -1);
  res[10] = 1 * (a[10] * b[15] + a[13] * -1 * b[14] * -1 -
                 a[14] * -1 * b[13] * -1 + a[15] * b[10]);
  res[9] = 1 * (a[9] * b[15] + a[12] * -1 * b[14] * -1 -
                a[14] * -1 * b[12] * -1 + a[15] * b[9]);
  res[8] = 1 * (a[8] * b[15] + a[11] * -1 * b[14] * -1 -
                a[14] * -1 * b[11] * -1 + a[15] * b[8]);
  res[7] = 1 * (a[7] * b[15] + a[12] * -1 * b[13] * -1 -
                a[13] * -1 * b[12] * -1 + a[15] * b[7]);
  res[6] = 1 * (a[6] * b[15] - a[11] * -1 * b[13] * -1 +
                a[13] * -1 * b[11] * -1 + a[15] * b[6]);
  res[5] = 1 * (a[5] * b[15] + a[11] * -1 * b[12] * -1 -
                a[12] * -1 * b[11] * -1 + a[15] * b[5]);
  res[4] = 1 * (a[4] * b[15] - a[7] * b[14] * -1 + a[9] * b[13] * -1 -
                a[10] * b[12] * -1 - a[12] * -1 * b[10] + a[13] * -1 * b[9] -
                a[14] * -1 * b[7] + a[15] * b[4]);
  res[3] = 1 * (a[3] * b[15] - a[6] * b[14] * -1 - a[8] * b[13] * -1 +
                a[10] * b[11] * -1 + a[11] * -1 * b[10] - a[13] * -1 * b[8] -
                a[14] * -1 * b[6] + a[15] * b[3]);
  res[2] = 1 * (a[2] * b[15] - a[5] * b[14] * -1 + a[8] * b[12] * -1 -
                a[9] * b[11] * -1 - a[11] * -1 * b[9] + a[12] * -1 * b[8] -
                a[14] * -1 * b[5] + a[15] * b[2]);
  res[1] = 1 * (a[1] * b[15] + a[5] * b[13] * -1 + a[6] * b[12] * -1 +
                a[7] * b[11] * -1 + a[11] * -1 * b[7] + a[12] * -1 * b[6] +
                a[13] * -1 * b[5] + a[15] * b[1]);
  res[0] = 1 * (a[0] * b[15] + a[1] * b[14] * -1 + a[2] * b[13] * -1 +
                a[3] * b[12] * -1 + a[4] * b[11] * -1 + a[5] * b[10] +
                a[6] * b[9] + a[7] * b[8] + a[8] * b[7] + a[9] * b[6] +
                a[10] * b[5] - a[11] * -1 * b[4] - a[12] * -1 * b[3] -
                a[13] * -1 * b[2] - a[14] * -1 * b[1] + a[15] * b[0]);
  return res;
};

//***********************
// multivector.Dot : res = a | b
// The inner product.
//***********************
 constexpr multivector operator|(multivector const& a, multivector const& b) {
  multivector res;
  res[0] = b[0] * a[0] + b[2] * a[2] + b[3] * a[3] + b[4] * a[4] - b[8] * a[8] -
           b[9] * a[9] - b[10] * a[10] - b[14] * a[14];
  res[1] = b[1] * a[0] + b[0] * a[1] - b[5] * a[2] - b[6] * a[3] - b[7] * a[4] +
           b[2] * a[5] + b[3] * a[6] + b[4] * a[7] + b[11] * a[8] +
           b[12] * a[9] + b[13] * a[10] + b[8] * a[11] + b[9] * a[12] +
           b[10] * a[13] + b[15] * a[14] - b[14] * a[15];
  res[2] = b[2] * a[0] + b[0] * a[2] - b[8] * a[3] + b[9] * a[4] + b[3] * a[8] -
           b[4] * a[9] - b[14] * a[10] - b[10] * a[14];
  res[3] = b[3] * a[0] + b[8] * a[2] + b[0] * a[3] - b[10] * a[4] -
           b[2] * a[8] - b[14] * a[9] + b[4] * a[10] - b[9] * a[14];
  res[4] = b[4] * a[0] - b[9] * a[2] + b[10] * a[3] + b[0] * a[4] -
           b[14] * a[8] + b[2] * a[9] - b[3] * a[10] - b[8] * a[14];
  res[5] = b[5] * a[0] - b[11] * a[3] + b[12] * a[4] + b[0] * a[5] -
           b[15] * a[10] - b[3] * a[11] + b[4] * a[12] - b[10] * a[15];
  res[6] = b[6] * a[0] + b[11] * a[2] - b[13] * a[4] + b[0] * a[6] -
           b[15] * a[9] + b[2] * a[11] - b[4] * a[13] - b[9] * a[15];
  res[7] = b[7] * a[0] - b[12] * a[2] + b[13] * a[3] + b[0] * a[7] -
           b[15] * a[8] - b[2] * a[12] + b[3] * a[13] - b[8] * a[15];
  res[8] = b[8] * a[0] + b[14] * a[4] + b[0] * a[8] + b[4] * a[14];
  res[9] = b[9] * a[0] + b[14] * a[3] + b[0] * a[9] + b[3] * a[14];
  res[10] = b[10] * a[0] + b[14] * a[2] + b[0] * a[10] + b[2] * a[14];
  res[11] = b[11] * a[0] + b[15] * a[4] + b[0] * a[11] - b[4] * a[15];
  res[12] = b[12] * a[0] + b[15] * a[3] + b[0] * a[12] - b[3] * a[15];
  res[13] = b[13] * a[0] + b[15] * a[2] + b[0] * a[13] - b[2] * a[15];
  res[14] = b[14] * a[0] + b[0] * a[14];
  res[15] = b[15] * a[0] + b[0] * a[15];
  return res;
};

//***********************
// multivector.Add : res = a + b
// Multivector addition
//***********************
 constexpr multivector operator+(multivector const& a, multivector const& b) {
  multivector res;
  res[0] = a[0] + b[0];
  res[1] = a[1] + b[1];
  res[2] = a[2] + b[2];
  res[3] = a[3] + b[3];
  res[4] = a[4] + b[4];
  res[5] = a[5] + b[5];
  res[6] = a[6] + b[6];
  res[7] = a[7] + b[7];
  res[8] = a[8] + b[8];
  res[9] = a[9] + b[9];
  res[10] = a[10] + b[10];
  res[11] = a[11] + b[11];
  res[12] = a[12] + b[12];
  res[13] = a[13] + b[13];
  res[14] = a[14] + b[14];
  res[15] = a[15] + b[15];
  return res;
};

//***********************
// multivector.Sub : res = a - b
// Multivector subtraction
//***********************
 constexpr multivector operator-(multivector const& a, multivector const& b) {
  multivector res;
  res[0] = a[0] - b[0];
  res[1] = a[1] - b[1];
  res[2] = a[2] - b[2];
  res[3] = a[3] - b[3];
  res[4] = a[4] - b[4];
  res[5] = a[5] - b[5];
  res[6] = a[6] - b[6];
  res[7] = a[7] - b[7];
  res[8] = a[8] - b[8];
  res[9] = a[9] - b[9];
  res[10] = a[10] - b[10];
  res[11] = a[11] - b[11];
  res[12] = a[12] - b[12];
  res[13] = a[13] - b[13];
  res[14] = a[14] - b[14];
  res[15] = a[15] - b[15];
  return res;
};

//***********************
// multivector.smul : res = a * b
// scalar/multivector multiplication
//***********************
 constexpr multivector operator*(double a, multivector const& b) {
  multivector res;
  res[0] = a * b[0];
  res[1] = a * b[1];
  res[2] = a * b[2];
  res[3] = a * b[3];
  res[4] = a * b[4];
  res[5] = a * b[5];
  res[6] = a * b[6];
  res[7] = a * b[7];
  res[8] = a * b[8];
  res[9] = a * b[9];
  res[10] = a * b[10];
  res[11] = a * b[11];
  res[12] = a * b[12];
  res[13] = a * b[13];
  res[14] = a * b[14];
  res[15] = a * b[15];
  return res;
};

//***********************
// multivector.muls : res = a * b
// multivector/scalar multiplication
//***********************
 constexpr multivector operator*(multivector const& a, double b) {
  multivector res;
  res[0] = a[0] * b;
  res[1] = a[1] * b;
  res[2] = a[2] * b;
  res[3] = a[3] * b;
  res[4] = a[4] * b;
  res[5] = a[5] * b;
  res[6] = a[6] * b;
  res[7] = a[7] * b;
  res[8] = a[8] * b;
  res[9] = a[9] * b;
  res[10] = a[10] * b;
  res[11] = a[11] * b;
  res[12] = a[12] * b;
  res[13] = a[13] * b;
  res[14] = a[14] * b;
  res[15] = a[15] * b;
  return res;
};

//***********************
// multivector.sadd : res = a + b
// scalar/multivector addition
//***********************
 constexpr multivector operator+(double a, multivector const& b) {
  multivector res;
  res[0] = a + b[0];
  res[1] = b[1];
  res[2] = b[2];
  res[3] = b[3];
  res[4] = b[4];
  res[5] = b[5];
  res[6] = b[6];
  res[7] = b[7];
  res[8] = b[8];
  res[9] = b[9];
  res[10] = b[10];
  res[11] = b[11];
  res[12] = b[12];
  res[13] = b[13];
  res[14] = b[14];
  res[15] = b[15];
  return res;
};

//***********************
// multivector.adds : res = a + b
// multivector/scalar addition
//***********************
 constexpr multivector operator+(multivector const& a, double b) {
  multivector res;
  res[0] = a[0] + b;
  res[1] = a[1];
  res[2] = a[2];
  res[3] = a[3];
  res[4] = a[4];
  res[5] = a[5];
  res[6] = a[6];
  res[7] = a[7];
  res[8] = a[8];
  res[9] = a[9];
  res[10] = a[10];
  res[11] = a[11];
  res[12] = a[12];
  res[13] = a[13];
  res[14] = a[14];
  res[15] = a[15];
  return res;
};

//***********************
// multivector.ssub : res = a - b
// scalar/multivector subtraction
//***********************
 constexpr multivector operator-(double a, multivector const& b) {
  multivector res;
  res[0] = a - b[0];
  res[1] = -b[1];
  res[2] = -b[2];
  res[3] = -b[3];
  res[4] = -b[4];
  res[5] = -b[5];
  res[6] = -b[6];
  res[7] = -b[7];
  res[8] = -b[8];
  res[9] = -b[9];
  res[10] = -b[10];
  res[11] = -b[11];
  res[12] = -b[12];
  res[13] = -b[13];
  res[14] = -b[14];
  res[15] = -b[15];
  return res;
};

//***********************
// multivector.subs : res = a - b
// multivector/scalar subtraction
//***********************
 constexpr multivector operator-(multivector const& a, double b) {
  multivector res;
  res[0] = a[0] - b;
  res[1] = a[1];
  res[2] = a[2];
  res[3] = a[3];
  res[4] = a[4];
  res[5] = a[5];
  res[6] = a[6];
  res[7] = a[7];
  res[8] = a[8];
  res[9] = a[9];
  res[10] = a[10];
  res[11] = a[11];
  res[12] = a[12];
  res[13] = a[13];
  res[14] = a[14];
  res[15] = a[15];
  return res;
};

// Define the basis blades
static constexpr auto e0 = make_vector<1>(1.);
static constexpr auto e1 = make_vector<2>(1.);
static constexpr auto e2 = make_vector<3>(1.);
static constexpr auto e3 = make_vector<4>(1.);
static constexpr auto e01 = e0 ^ e1;
static constexpr auto e02 = e0 ^ e2;
static constexpr auto e03 = e0 ^ e3;
static constexpr auto e12 = e1 ^ e2;
static constexpr auto e31 = e3 ^ e1;
static constexpr auto e23 = e2 ^ e3;
static constexpr auto e021 = e0 ^ e2 ^ e1;
static constexpr auto e013 = e0 ^ e1 ^ e3;
static constexpr auto e032 = e0 ^ e3 ^ e2;
static constexpr auto e123 = e1 ^ e2 ^ e3;
static constexpr auto e0123 = e0 ^ e1 ^ e2 ^ e3;

 constexpr double norm(multivector const& a) {
  return std::sqrt(std::abs((a * conjugate(a))[0]));
}
 constexpr double inorm(multivector const& a) {
  return norm(!a);
}
 constexpr multivector normalized(multivector const& a) {
  return a * (1. / norm(a));
}

// A rotor (Euclidean line) and translator (Ideal line)
constexpr multivector rotor(double angle, multivector line) {
  return std::cos(angle / 2.) + std::sin(angle / 2.) * normalized(line);
}
constexpr multivector translator(double dist, multivector line) {
  return 1. + dist / 2. * line;
}


// A plane is defined using its homogeneous equation ax + by + cz + d = 0
constexpr multivector plane(double a, double b, double c, double d) {
  return a * e1 + b * e2 + c * e3 + d * e0;
}


// A point is just a homogeneous point, euclidean coordinates plus the origin
constexpr multivector point(double x, double y, double z) {
  return e123 + x * e032 + y * e013 + z * e021;
}

// for our toy problem (generate points on the surface of a torus)
// we start with a function that generates motors.
// circle(t) with t going from 0 to 1.
constexpr multivector circle(double t, double radius, multivector const& line) {
  return rotor(t * 2. * pi, line) * translator(radius, e1 * e0);
}

// a torus is now the product of two circles.
constexpr multivector torus(double s, double t, double r1, multivector const& l1, double r2, multivector const& l2) {
  return circle(s, r2, l2) * circle(t, r1, l1);
}

// and to sample its points we simply sandwich the origin ..
constexpr multivector point_on_torus(double s, double t) {
  multivector to = torus(s, t, 0.25f, e1 * e2, 0.6f, e1 * e3);
  return to * e123 * ~to;
}

void test() {

  // Elements of the even subalgebra (scalar + bivector + pss) of unit length
  // are motors
  [[maybe_unused]] auto const rot = rotor(pi / 2., e1 * e2);

  // The outer product ^ is the MEET. Here we intersect the yz (x=0) and xz
  // (y=0) planes.
  [[maybe_unused]] auto const ax_z = e1 ^ e2;

  // line and plane meet in point. We intersect the line along the z-axis
  // (x=0,y=0) with the xy (z=0) plane.
  [[maybe_unused]] auto const orig = ax_z ^ e3;

  // We can also easily create points and join them into a line using the
  // regressive (vee, &) product.
  [[maybe_unused]] auto const px = point(1.0, 0.0, 0.0);
  [[maybe_unused]] auto const line = orig & px;

  // Lets also create the plane with equation 2x + z - 3 = 0
  [[maybe_unused]] auto const p = plane(2, 0, 1, -3);

  // rotations work on all elements
  [[maybe_unused]] auto const rotated_plane = rot * p * ~rot;
  [[maybe_unused]] auto const rotated_line = rot * line * ~rot;
  [[maybe_unused]] auto const rotated_point = rot * px * ~rot;

  // See the 3D PGA Cheat sheet for a huge collection of useful formulas
  [[maybe_unused]] auto const point_on_plane = (p | px) * p;

  // Some output.
  // printf("a point       : ");
  // px.log();
  // printf("a line        : ");
  // line.log();
  // printf("a plane       : ");
  // p.log();
  // printf("a rotor       : ");
  // rot.log();
  // printf("rotated line  : ");
  // rotated_line.log();
  // printf("rotated point : ");
  // rotated_point.log();
  // printf("rotated plane : ");
  // rotated_plane.log();
  // printf("point on plane: ");
  // point_on_plane.normalized().log();
  // printf("point on torus: ");
  // point_on_torus(0., 0.).log();
  // (e0 - 1.).log();
  // (1. - e0).log();
}
} // namespace pga
