#include "math_util.h"
#include "backend.h"
#include "gtest/gtest.h"
#include "unit/test_util.h"

namespace bdm {

TEST(MathUtilTest, Norm) {
  std::array<double, 3> vector = {1.1, 2.2, 3.3};
  auto result = Math::Norm(vector);

  EXPECT_NEAR(16.94, result, abs_error<double>::value);
}

TEST(MathUtilTest, NormZero) {
  std::array<double, 3> vector = {0, 0, 0};
  auto result = Math::Norm(vector);

  EXPECT_NEAR(1, result, abs_error<double>::value);
}

TEST(MathUtilTest, NormalizeZero) {
  std::array<double, 3> vector = {0, 0, 0};
  auto result = Math::Normalize(vector);

  EXPECT_NEAR(0, result[0], abs_error<double>::value);
  EXPECT_NEAR(0, result[1], abs_error<double>::value);
  EXPECT_NEAR(0, result[2], abs_error<double>::value);
}

TEST(MathUtilTest, Normalize) {
  std::array<double, 3> vector = {1.1, 2.2, 3.3};
  auto result = Math::Normalize(vector);

  EXPECT_NEAR(0.0649350649350649351, result[0], abs_error<double>::value);
  EXPECT_NEAR(0.1298701298701298701, result[1], abs_error<double>::value);
  EXPECT_NEAR(0.1948051948051948052, result[2], abs_error<double>::value);
}

TEST(MathUtilTest, L2Distance) {
  std::array<double, 3> vector1 = {0, 0, 0};
  std::array<double, 3> vector2 = {1, 2, 3};
  auto result1 = Math::GetL2Distance(vector1, vector2);
  auto result2 = Math::GetL2Distance(vector2, vector1);

  EXPECT_NEAR(3.7416573867739413855, result1, abs_error<double>::value);
  EXPECT_NEAR(3.7416573867739413855, result2, abs_error<double>::value);
}

}  // namespace bdm
