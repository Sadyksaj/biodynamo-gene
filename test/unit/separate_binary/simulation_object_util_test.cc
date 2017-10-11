#include "gtest/gtest.h"

#include "unit/separate_binary/simulation_object_util_test.h"

namespace bdm {
// namespace simulation_object_util_test_internal {

// The following tests check if code insertion in new classes works as intended
// Therefore SimulationObject is extended in two stages: first by CellExt and
// then by NeuronExt

TEST(SimulationObjectUtilTest, ContainerFunctionality) {
  using Scalar = ContainerTestClass;
  auto vector = Scalar::NewEmptySoa();

  EXPECT_EQ(0u, vector.size());
  EXPECT_EQ(0u, vector.GetTotalSize());

  EXPECT_EQ(0u, vector.DelayedPushBack(Scalar(1, 1.0)));
  EXPECT_EQ(1u, vector.DelayedPushBack(Scalar(2, 2.0)));
  EXPECT_EQ(2u, vector.DelayedPushBack(Scalar(3, 3.0)));

  // changes have not been commited yet
  EXPECT_EQ(0u, vector.size());
  //   but data is already in vector
  EXPECT_EQ(3u, vector.GetVecDm1().size());
  EXPECT_EQ(3u, vector.GetVecDm2().size());
  EXPECT_EQ(3u, vector.GetTotalSize());

  vector.Commit();

  EXPECT_EQ(3u, vector.size());
  EXPECT_EQ(3u, vector.GetVecDm1().size());
  EXPECT_EQ(3u, vector.GetVecDm2().size());
  EXPECT_EQ(3u, vector.GetTotalSize());

  EXPECT_EQ(1, vector[0].GetDm1());
  EXPECT_EQ(1.0, vector[0].GetDm2());
  EXPECT_EQ(2, vector[1].GetDm1());
  EXPECT_EQ(2.0, vector[1].GetDm2());
  EXPECT_EQ(3, vector[2].GetDm1());
  EXPECT_EQ(3.0, vector[2].GetDm2());

  vector.DelayedRemove(0);

  // changes have not been commited yet
  EXPECT_EQ(3u, vector.size());
  EXPECT_EQ(3u, vector.GetVecDm1().size());
  EXPECT_EQ(3u, vector.GetVecDm2().size());
  EXPECT_EQ(3u, vector.GetTotalSize());

  vector.Commit();

  EXPECT_EQ(2u, vector.size());
  EXPECT_EQ(2u, vector.GetVecDm1().size());
  EXPECT_EQ(2u, vector.GetVecDm2().size());
  EXPECT_EQ(2u, vector.GetTotalSize());

  EXPECT_EQ(3, vector[0].GetDm1());
  EXPECT_EQ(3.0, vector[0].GetDm2());
  EXPECT_EQ(2, vector[1].GetDm1());
  EXPECT_EQ(2.0, vector[1].GetDm2());

  vector.DelayedRemove(0);
  vector.DelayedRemove(1);
  vector.Commit();
  EXPECT_EQ(0u, vector.size());
  EXPECT_EQ(0u, vector.GetTotalSize());

  // push_back
  vector.push_back(Scalar(9, 9.0));
  EXPECT_EQ(1u, vector.size());
  EXPECT_EQ(1u, vector.GetTotalSize());
  EXPECT_EQ(9, vector[0].GetDm1());

  vector.DelayedPushBack(Scalar(10, 10.0));
  EXPECT_EQ(1u, vector.size());
  EXPECT_EQ(2u, vector.GetTotalSize());
  try {
    vector.push_back(Scalar(11, 11.0));
    FAIL() << "Should have thrown a logic_error";
  } catch(std::logic_error& e) {}
}

template <typename T>
void RunDefaultConstructorTest(const T& neuron) {
  EXPECT_EQ(1u, neuron.size());

  EXPECT_EQ(6.28, neuron.GetDiameter());
  auto& position = neuron.GetPosition();
  EXPECT_EQ(1, position[0]);
  EXPECT_EQ(2, position[1]);
  EXPECT_EQ(3, position[2]);
  auto neurites = neuron.GetNeurites();
  EXPECT_EQ(0u, neurites.size());
}

TEST(SimulationObjectUtilTest, DefaultConstructor) {
  // are data members in all extensions correctly initialized?
  Neuron neuron;
  RunDefaultConstructorTest(Neuron());
  RunDefaultConstructorTest(SoaNeuron());
}

TEST(SimulationObjectUtilTest, NewEmptySoa) {
  // Create an empty SOA container
  // Creating it using e.g. `SoaNeuron soa;` will already have one element
  // inside - (the one with default parameters)
  auto neurons = Neuron::NewEmptySoa();
  EXPECT_EQ(0u, neurons.size());
  EXPECT_EQ(0u, neurons.neurites_.size());
  EXPECT_EQ(0u, neurons.diameter_.size());
  EXPECT_EQ(0u, neurons.position_.size());
}

TEST(SimulationObjectUtilTest, NonDefaultConstructor) {
  std::vector<Neurite> neurites;
  neurites.push_back(Neurite(2));
  neurites.push_back(Neurite(3));

  Neuron neuron(neurites, std::array<double, 3>{4, 5, 6});

  EXPECT_EQ(6.28, neuron.GetDiameter());
  auto& position = neuron.GetPosition();
  EXPECT_EQ(4, position[0]);
  EXPECT_EQ(5, position[1]);
  EXPECT_EQ(6, position[2]);
  EXPECT_EQ(2u, neurites.size());
}

TEST(SimulationObjectUtilTest, SoaRef) {
  SoaNeuron neurons;

  auto neurons_ref = neurons[0];

  EXPECT_EQ(1u, neurons_ref.size());

  // check if changes are visible for the referenced object
  neurons_ref.SetDiameter(12.34);
  EXPECT_EQ(12.34, neurons.GetDiameter());
  neurons_ref.push_back(Neuron());
  EXPECT_EQ(2u, neurons.size());
}

TEST(SimulationObjectUtilTest, Soa_push_back_AndSubscriptOperator) {
  std::vector<Neurite> neurites;
  neurites.push_back(Neurite(2));
  neurites.push_back(Neurite(3));

  Neuron neuron1(neurites, std::array<double, 3>{4, 5, 6});

  neurites.push_back(Neurite(4));
  Neuron neuron2(neurites, std::array<double, 3>{9, 8, 7});

  auto neurons = Neuron::NewEmptySoa();
  neurons.push_back(neuron1);
  neurons.push_back(neuron2);

  EXPECT_EQ(2u, neurons.size());

  EXPECT_EQ(6.28, neurons[0].GetDiameter());
  auto& position1 = neurons[0].GetPosition();
  EXPECT_EQ(4, position1[0]);
  EXPECT_EQ(5, position1[1]);
  EXPECT_EQ(6, position1[2]);
  EXPECT_EQ(2u, neurons[0].GetNeurites().size());

  // test if return type of subscript operator has SoaRef backend
  auto element1 = neurons[1];
  Neuron::Self<SoaRef>* cast_result =
      dynamic_cast<Neuron::Self<SoaRef>*>(&element1);
  EXPECT_TRUE(cast_result != nullptr);

  EXPECT_EQ(6.28, neurons[1].GetDiameter());
  auto& position2 = neurons[1].GetPosition();
  EXPECT_EQ(9, position2[0]);
  EXPECT_EQ(8, position2[1]);
  EXPECT_EQ(7, position2[2]);
  EXPECT_EQ(3u, neurons[1].GetNeurites().size());
}

TEST(SimulationObjectUtilTest, Soa_clear) {
  SoaNeuron neurons;
  EXPECT_EQ(1u, neurons.size());
  neurons.clear();
  EXPECT_EQ(0u, neurons.size());
  EXPECT_EQ(0u, neurons.neurites_.size());
  EXPECT_EQ(0u, neurons.diameter_.size());
  EXPECT_EQ(0u, neurons.position_.size());
}

TEST(SimulationObjectUtilTest, Soa_reserve) {
  SoaNeuron neurons;
  neurons.reserve(10);
  EXPECT_EQ(10u, neurons.neurites_.capacity());
  EXPECT_EQ(10u, neurons.diameter_.capacity());
  EXPECT_EQ(10u, neurons.position_.capacity());
}

TEST(SimulationObjectUtilTest, Soa_AssignmentOperator) {
  std::vector<Neurite> neurites;
  neurites.push_back(Neurite(2));
  neurites.push_back(Neurite(3));

  Neuron neuron1(neurites, std::array<double, 3>{4, 5, 6});

  neurites.push_back(Neurite(4));
  Neuron new_neuron1(neurites, std::array<double, 3>{9, 8, 7});
  new_neuron1.SetDiameter(123);

  auto neurons = Neuron::NewEmptySoa();
  neurons.push_back(neuron1);

  EXPECT_EQ(1u, neurons.size());

  neurons[0] = new_neuron1;
  EXPECT_EQ(123u, neurons[0].GetDiameter());
  auto& position = neurons[0].GetPosition();
  EXPECT_EQ(9, position[0]);
  EXPECT_EQ(8, position[1]);
  EXPECT_EQ(7, position[2]);
  EXPECT_EQ(3u, neurons[0].GetNeurites().size());
}

template <typename TContainer>
void RunDivideTest(TContainer* neurons) {
  Neuron neuron;
  neurons->push_back(neuron);

  auto new_neuron_idx = Divide((*neurons)[0], neurons, 1.0, 2.0, 3.0);

  EXPECT_EQ(1u, new_neuron_idx);
  EXPECT_EQ(1u, neurons->size()); // not visible yet
  EXPECT_EQ(987u, (*neurons)[new_neuron_idx].GetNeurites()[0].id_);
  EXPECT_EQ(5, (*neurons)[new_neuron_idx].GetPosition()[0]);
  EXPECT_EQ(4, (*neurons)[new_neuron_idx].GetPosition()[1]);
  EXPECT_EQ(3, (*neurons)[new_neuron_idx].GetPosition()[2]);

  // commit invalidates new_neuron
  neurons->Commit();

  ASSERT_EQ(2u, neurons->size());
  // new_neuron got invalidated by `Commit()`, but is now accessible in neurons
  EXPECT_EQ(987u, (*neurons)[1].GetNeurites()[0].id_);
  EXPECT_EQ(5, (*neurons)[1].GetPosition()[0]);
  EXPECT_EQ(4, (*neurons)[1].GetPosition()[1]);
  EXPECT_EQ(3, (*neurons)[1].GetPosition()[2]);
  EXPECT_EQ(1.123, (*neurons)[0].GetDiameter());
}

TEST(SimulationObjectUtilTest, Aos_Divide) {
  TransactionalVector<Neuron> neurons;
  RunDivideTest(&neurons);
}

TEST(SimulationObjectUtilTest, Soa_Divide) {
  auto neurons = Neuron::NewEmptySoa();
  RunDivideTest(&neurons);
}

// Tests overloaded Divide function which adds new daughter cell to the
// container managed by the ResourceManager with default template parameters
TEST(SimulationObjectUtilTest, Soa_DivideWithResourceManager) {
  auto rm = ResourceManager<>::Get();
  rm->Clear();

  auto neurons = rm->Get<Neuron>();
  Neuron neuron;
  neurons->push_back(neuron);

  auto&& new_neuron = Divide((*neurons)[0], 1.0, 2.0, 3.0);

  EXPECT_EQ(987u, new_neuron.Get().GetNeurites()[0].id_);
  EXPECT_EQ(5, new_neuron.Get().GetPosition()[0]);
  EXPECT_EQ(4, new_neuron.Get().GetPosition()[1]);
  EXPECT_EQ(3, new_neuron.Get().GetPosition()[2]);

  // commit invalidates new_neuron
  neurons->Commit();

  ASSERT_EQ(2u, neurons->size());
  // new_neuron got invalidated by `Commit()`, but is now accessible in neurons
  EXPECT_EQ(987u, (*neurons)[1].GetNeurites()[0].id_);
  EXPECT_EQ(5, (*neurons)[1].GetPosition()[0]);
  EXPECT_EQ(4, (*neurons)[1].GetPosition()[1]);
  EXPECT_EQ(3, (*neurons)[1].GetPosition()[2]);
  EXPECT_EQ(1.123, (*neurons)[0].GetDiameter());
}

template <typename TContainer>
void RunDeleteTest(TContainer* neurons) {
  Neuron neuron;
  neurons->push_back(neuron);

  Delete(neurons, 0);

  neurons->Commit();

  EXPECT_EQ(0u, neurons->size());
}

TEST(SimulationObjectUtilTest, Aos_Delete) {
  TransactionalVector<Neuron> neurons;
  RunDeleteTest(&neurons);
}

TEST(SimulationObjectUtilTest, Soa_Delete) {
  auto neurons = Neuron::NewEmptySoa();
  RunDeleteTest(&neurons);
}

TEST(SimulationObjectUtilTest, Soa_IO) { RunSoaIOTest(); }

// ----------------------------------------------------------------------------
// SoPointer tests
/// This function is before the decleration of `SoPointerTestClass` to test
/// if the SoPointer implementation can deal with incomplete types
template <typename T, typename TBackend>
void RunSoPointerTest(T* sim_objects) {
  using SO = typename T::value_type;

  SoPointer<SO, TBackend> null_so_pointer;
  EXPECT_TRUE(null_so_pointer.IsNullPtr());

  SoPointer<SO, TBackend> so_ptr(sim_objects, 0);

  EXPECT_FALSE(so_ptr.IsNullPtr());
  EXPECT_EQ(123u, so_ptr.Get().GetId());
}

BDM_SIM_CLASS(SoPointerTestClass, SimulationObject) {
  BDM_CLASS_HEADER(SoPointerTestClassExt, 1, id_);

 public:
  SoPointerTestClassExt() {}
  SoPointerTestClassExt(uint64_t id) { id_[kIdx] = id; }

  uint64_t GetId() const { return id_[kIdx]; }

 private:
  vec<uint64_t> id_;
};

TEST(SimulationObjectUtilTest, Aos_SoPointer) {
  TransactionalVector<SoPointerTestClass> sim_objects;
  SoPointerTestClass so(123);
  sim_objects.push_back(so);
  RunSoPointerTest<decltype(sim_objects), Scalar>(&sim_objects);
}

TEST(SimulationObjectUtilTest, Soa_SoPointer) {
  auto sim_objects = SoPointerTestClass::NewEmptySoa();
  SoPointerTestClass so(123);
  sim_objects.push_back(so);
  RunSoPointerTest<decltype(sim_objects), Soa>(&sim_objects);
}

// }  // namespace simulation_object_util_test_internal
}  // namespace bdm

int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
