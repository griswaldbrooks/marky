#include <gtest/gtest.h>
#include "landy/pga.hpp"
namespace pga {

multivector make_zero_multivector() {
  auto m = multivector{};
  for (auto& element : m) {
    element = 0.;
  }
  return m;
}

TEST(Multivector, ZeroInitialized) {
 auto const a = multivector{};
  auto const result = make_zero_multivector();
  EXPECT_EQ(a, result);
}
} // namespace landy::pga
