#ifndef UNIT_SIMULATION_OBJECT_TEST_H_
#define UNIT_SIMULATION_OBJECT_TEST_H_

#include "gtest/gtest.h"
#include "simulation_object.h"

namespace bdm {
namespace simulation_object_test_internal {

template <typename TBackend = Soa>
struct CTParam {
  template <typename TTBackend>
  using Self = CTParam<TTBackend>;
  using Backend = TBackend;
};

inline void RunPushBackAndClearTest() {
  SimulationObject<CTParam<>> soa;
  // call clear, because creating a SOA object with default constructor will
  // already have one element inside
  soa.clear();
  EXPECT_EQ(0u, soa.size());
  SimulationObject<CTParam<Scalar>> so;
  soa.push_back(so);
  soa.push_back(so);
  EXPECT_EQ(2u, soa.size());
  soa.clear();
  EXPECT_EQ(0u, soa.size());
}

}  // namespace simulation_object_test_internal
}  // namespace bdm

#endif  // UNIT_SIMULATION_OBJECT_TEST_H_
