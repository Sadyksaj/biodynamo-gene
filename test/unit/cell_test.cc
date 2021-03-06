#include "cell.h"
#include "gtest/gtest.h"
#include "unit/cell_test.h"
#include "unit/test_util.h"

namespace bdm {
namespace cell_test_internal {

TEST(CellTest, TransformCoordinatesGlobalToPolar) {
  TestCell cell;
  cell.TestTransformCoordinatesGlobalToPolar();
}

TEST(CellTest, DivideVolumeRatioPhiTheta) {
  TestCell mother;
  mother.SetPosition({5, 6, 7});
  mother.SetMassLocation({5, 6, 7});
  mother.SetTractorForce({0, 0, 0});
  mother.SetDiameter(10);
  mother.UpdateVolume();
  mother.SetAdherence(1.1);
  mother.SetMass(5);
  mother.SetXAxis({1, 2, 3});
  mother.SetYAxis({4, 5, 6});
  mother.SetZAxis({7, 8, 9});
  mother.AddBiologyModule(GrowthModule());
  mother.AddBiologyModule(MovementModule({1, 2, 3}));
  mother.SetBoxIdx(123);

  TestCell daughter;
  mother.Divide(&daughter, 0.75, 0.12, 0.34);

  const double kEpsilon = abs_error<double>::value;

  // verify mother data members
  EXPECT_NEAR(0.16369893089539111, mother.GetPosition()[0], kEpsilon);
  EXPECT_NEAR(0.39656253332927616, mother.GetPosition()[1], kEpsilon);
  EXPECT_NEAR(0.6294261357631612, mother.GetPosition()[2], kEpsilon);

  EXPECT_NEAR(0.16369893089539111, mother.GetMassLocation()[0], kEpsilon);
  EXPECT_NEAR(0.39656253332927616, mother.GetMassLocation()[1], kEpsilon);
  EXPECT_NEAR(0.6294261357631612, mother.GetMassLocation()[2], kEpsilon);

  EXPECT_NEAR(0, mother.GetTractorForce()[0], kEpsilon);
  EXPECT_NEAR(0, mother.GetTractorForce()[1], kEpsilon);
  EXPECT_NEAR(0, mother.GetTractorForce()[2], kEpsilon);

  EXPECT_NEAR(8.2982653336624335, mother.GetDiameter(), kEpsilon);
  // differs slightly from the value in branch validation due to more precise
  // value of PI
  EXPECT_NEAR(299.19930034188491, mother.GetVolume(), kEpsilon);
  EXPECT_NEAR(1.1, mother.GetAdherence(), kEpsilon);
  EXPECT_NEAR(2.8571428571428563, mother.GetMass(), kEpsilon);

  EXPECT_NEAR(1, mother.GetXAxis()[0], kEpsilon);
  EXPECT_NEAR(2, mother.GetXAxis()[1], kEpsilon);
  EXPECT_NEAR(3, mother.GetXAxis()[2], kEpsilon);

  EXPECT_NEAR(4, mother.GetYAxis()[0], kEpsilon);
  EXPECT_NEAR(2, mother.GetXAxis()[1], kEpsilon);
  EXPECT_NEAR(3, mother.GetXAxis()[2], kEpsilon);

  EXPECT_NEAR(7, mother.GetZAxis()[0], kEpsilon);
  EXPECT_NEAR(8, mother.GetZAxis()[1], kEpsilon);
  EXPECT_NEAR(9, mother.GetZAxis()[2], kEpsilon);

  EXPECT_EQ(123u, mother.GetBoxIdx());

  // verify daughter data members
  EXPECT_NEAR(11.448401425472813, daughter.GetPosition()[0], kEpsilon);
  EXPECT_NEAR(13.471249955560966, daughter.GetPosition()[1], kEpsilon);
  EXPECT_NEAR(15.494098485649118, daughter.GetPosition()[2], kEpsilon);

  EXPECT_NEAR(11.448401425472813, daughter.GetMassLocation()[0], kEpsilon);
  EXPECT_NEAR(13.471249955560966, daughter.GetMassLocation()[1], kEpsilon);
  EXPECT_NEAR(15.494098485649118, daughter.GetMassLocation()[2], kEpsilon);

  EXPECT_NEAR(0, daughter.GetTractorForce()[0], kEpsilon);
  EXPECT_NEAR(0, daughter.GetTractorForce()[1], kEpsilon);
  EXPECT_NEAR(0, daughter.GetTractorForce()[2], kEpsilon);

  EXPECT_NEAR(7.5394744112915388, daughter.GetDiameter(), kEpsilon);
  // differs slightly from the value in branch validation due to more precise
  // value of PI
  EXPECT_NEAR(224.39947525641387, daughter.GetVolume(), kEpsilon);
  EXPECT_NEAR(1.1, daughter.GetAdherence(), kEpsilon);
  EXPECT_NEAR(2.1428571428571437, daughter.GetMass(), kEpsilon);

  EXPECT_NEAR(1, daughter.GetXAxis()[0], kEpsilon);
  EXPECT_NEAR(2, daughter.GetXAxis()[1], kEpsilon);
  EXPECT_NEAR(3, daughter.GetXAxis()[2], kEpsilon);

  EXPECT_NEAR(4, daughter.GetYAxis()[0], kEpsilon);
  EXPECT_NEAR(2, daughter.GetXAxis()[1], kEpsilon);
  EXPECT_NEAR(3, daughter.GetXAxis()[2], kEpsilon);

  EXPECT_NEAR(7, daughter.GetZAxis()[0], kEpsilon);
  EXPECT_NEAR(8, daughter.GetZAxis()[1], kEpsilon);
  EXPECT_NEAR(9, daughter.GetZAxis()[2], kEpsilon);

  // biology modules mother
  EXPECT_EQ(2u, mother.GetBiologyModules().size());
  EXPECT_EQ(1u, daughter.GetBiologyModules().size());
  if (get_if<GrowthModule>(&(daughter.GetBiologyModules()[0])) == nullptr) {
    FAIL() << "Variant type at position 0 is not a GrowthModule";
  }

  EXPECT_EQ(123u, daughter.GetBoxIdx());

  // additional check
  EXPECT_NEAR(5, mother.GetMass() + daughter.GetMass(), kEpsilon);
}

TEST(CellTest, Divide) {
  TestCell cell;
  gRandom.SetSeed(42);

  cell.check_input_parameters_ = true;
  cell.expected_volume_ratio_ = 1.0455127360065737;
  cell.expected_phi_ = 1.9633629889829609;
  cell.expected_theta_ = 4.2928196812086608;

  TestCell daughter;
  cell.Divide(&daughter);
}

TEST(CellTest, DivideVolumeRatio) {
  TestCell cell;
  gRandom.SetSeed(42);

  cell.check_input_parameters_ = true;
  cell.expected_volume_ratio_ = 0.59;
  cell.expected_phi_ = 1.1956088797871529;
  cell.expected_theta_ = 4.5714174264720571;

  TestCell daughter;
  cell.Divide(&daughter, 0.59);
}

TEST(CellTest, DivideAxis) {
  TestCell cell;
  cell.SetMassLocation({1, 2, 3});
  gRandom.SetSeed(42);

  cell.check_input_parameters_ = true;
  cell.expected_volume_ratio_ = 1.0455127360065737;
  cell.expected_phi_ = 1.0442265974045177;
  cell.expected_theta_ = 0.72664234068172562;

  TestCell daughter;
  cell.Divide(&daughter, {9, 8, 7});
}

TEST(CellTest, DivideVolumeRatioAxis) {
  TestCell cell;
  cell.SetMassLocation({1, 2, 3});
  gRandom.SetSeed(42);

  cell.check_input_parameters_ = true;
  cell.expected_volume_ratio_ = 0.456;
  cell.expected_phi_ = 1.0442265974045177;
  cell.expected_theta_ = 0.72664234068172562;

  TestCell daughter;
  cell.Divide(&daughter, 0.456, {9, 8, 7});
}

TEST(CellTest, BiologyModule) {
  TestCell cell;
  double diameter = cell.GetDiameter();
  auto position = cell.GetPosition();

  cell.AddBiologyModule(MovementModule({1, 2, 3}));
  cell.AddBiologyModule(GrowthModule());

  cell.RunBiologyModules<TestCell>();

  EXPECT_NEAR(diameter + 0.5, cell.GetDiameter(), abs_error<double>::value);
  EXPECT_NEAR(position[0] + 1, cell.GetPosition()[0], abs_error<double>::value);
  EXPECT_NEAR(position[1] + 2, cell.GetPosition()[1], abs_error<double>::value);
  EXPECT_NEAR(position[2] + 3, cell.GetPosition()[2], abs_error<double>::value);
}

TEST(CellTest, IO) { RunIOTest(); }

}  // namespace cell_test_internal
}  // namespace bdm
