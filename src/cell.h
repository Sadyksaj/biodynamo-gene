#ifndef CELL_H_
#define CELL_H_

#include <array>
#include <cmath>
#include <type_traits>
#include <vector>

#include "backend.h"
#include "biology_module_util.h"
#include "cell.h"
#include "default_force.h"
#include "inline_vector.h"
#include "math_util.h"
#include "matrix.h"
#include "bdm_param.h"
#include "simulation_object.h"
#include "simulation_object_util.h"

namespace bdm {

using std::array;
using std::vector;

template <typename Base = SimulationObject<>,
          typename TBiologyModuleVariant = variant<NullBiologyModule>>
class CellExt : public Base {
 public:
  BDM_CLASS_HEADER_ADV(CellExt,
                       CellExt<typename Base::template Self<TTBackend> COMMA()
                                   TBiologyModuleVariant>,
                       template <typename COMMA() typename>, position_, gpu_position_,
                       mass_location_, tractor_force_, diameter_, volume_,
                       adherence_, density_, x_axis_, y_axis_, z_axis_,
                       neighbors_, biology_modules_);

 public:
  CellExt() {}
  explicit CellExt(float diameter) : diameter_(diameter) { UpdateVolume(); }
  explicit CellExt(const array<float, 3>& position)
      : position_(position), mass_location_(position) {}

  virtual ~CellExt() {}

  /// Add a biology module to this cell
  /// @tparam TBiologyModule type of the biology module. Must be in the set of
  ///         types specified in `TBiologyModuleVariant`
  template <typename TBiologyModule>
  void AddBiologyModule(TBiologyModule&& module);

  /// Execute all biology modules
  void RunBiologyModules();

  /// Divide the cell. Of the two daughter cells, one is this one (but smaller,
  /// with half GeneSubstances etc.),
  /// and the other one is given as parameter and initialzed accordingly. Both
  /// cells have more or less the same volume,
  /// the axis of division is random.
  /// @param daughter - second daughter cell = scalar instance which will be
  /// initialized in this method
  void Divide(Self<Scalar>* daughter);

  /// Divide the cell. Of the two daughter cells, one is this one (but smaller,
  /// with half GeneSubstances etc.),
  /// and the other one is given as parameter and initialzed accordingly. The
  /// axis of division is random.
  /// @param daughter second daughter cell = scalar instance which will be
  /// initialized in this method
  /// @param volume_ratio the ratio (Volume daughter 1)/(Volume daughter 2). 1.0
  /// gives equal cells.
  void Divide(Self<Scalar>* daughter, float volume_ratio);

  /// Divide the cell. Of the two daughter cells, one is this one (but smaller,
  /// with half GeneSubstances etc.),
  /// and the other one is is given as parameter and initialzed accordingly.
  /// @param daughter second daughter cell = scalar instance which will be
  /// initialized in this method
  /// @param axis specifies direction of division
  void Divide(Self<Scalar>* daughter, const array<float, 3>& axis);

  /// Divide the cell. Of the two daughter cells, one is this one (but smaller,
  /// with half GeneSubstances etc.),
  /// and the other one is instantiated de novo and is returned.
  /// @param daughter second daughter cell = scalar instance which will be
  /// initialized in this method
  /// @param volume_ratio the ratio (Volume daughter 1)/(Volume daughter 2). 1.0
  /// gives equal cells.
  /// @param axis specifies direction of division
  void Divide(Self<Scalar>* daughter, float volume_ratio,
              const array<float, 3>& axis);

  /// Forwards call to DivideImpl @see DivideImpl
  void Divide(Self<Scalar>* daughter, float volume_ratio, float phi,
              float theta);

  /// Divide mother cell in two daughter cells\n
  /// When mother cell divides, by definition:\n
  /// 1) the mother cell becomes the 1st daughter cell\n
  /// 2) the new cell becomes the 2nd daughter cell and inherits a equal or
  /// bigger volume than the 1st
  ///    daughter cell, which means that this cell will eventually inherit more
  ///    differentiating factors
  ///    and will be recorded in the left side of the lineage tree.
  ///
  /// @param daughter second daughter cell = scalar instance which will be
  /// initialized in this method
  /// @param volume_ratio the ratio (Volume daughter 1)/(Volume daughter 2). 1.0
  /// gives equal cells.
  /// @param phi azimuthal angle (polar coordinate)
  /// @param theta polar angle (polar coordinate)
  /// @see \link simulation_object_util.h Divide \endlink
  virtual void DivideImpl(Self<Scalar>* daughter, float volume_ratio,
                          float phi, float theta);

  float GetAdherence() const { return adherence_[kIdx]; }

  float GetDiameter() const { return diameter_[kIdx]; }

  float GetMass() const { return density_[kIdx] * volume_[kIdx]; }

  float GetDensity() const { return density_[kIdx]; }

  const array<float, 3>& GetMassLocation() const {
    return mass_location_[kIdx];
  }

  const array<float, 3>& GetPosition() const { return position_[kIdx]; }

  array<float, 3>& GetPosition() { return position_[kIdx]; }

  vec<array<float, 3>>& GetAllPositions() { return position_; }

  float* GetGpuPositions() { return &((gpu_position_.data())[0][0]); }

  const array<float, 3>& GetTractorForce() const {
    return tractor_force_[kIdx];
  }

  float GetVolume() const { return volume_[kIdx]; }

  const InlineVector<int, 8>& GetNeighbors() const { return neighbors_[kIdx]; }

  void SetAdherence(float adherence) { adherence_[kIdx] = adherence; }

  void SetDiameter(float diameter) {
    diameter_[kIdx] = diameter;
    UpdateVolume();
  }

  void SetMass(float mass) { density_[kIdx] = mass / volume_[kIdx]; }

  void SetDensity(float density) { density_[kIdx] = density; }

  void SetGpuPosition(const array<float, 4>& gpu_position) {
    gpu_position_[kIdx] = gpu_position;
  }

  void SetMassLocation(const array<float, 3>& mass_location) {
    mass_location_[kIdx] = mass_location;
  }

  void SetPosition(const array<float, 3>& position) {
    position_[kIdx] = position;
  }

  void SetTractorForce(const array<float, 3>& tractor_force) {
    tractor_force_[kIdx] = tractor_force;
  }

  void SetNeighbors(const InlineVector<int, 8>& neighbors) {
    neighbors_[kIdx] = neighbors;
  }

  void ChangeVolume(float speed) {
    // scaling for integration step
    float delta = speed * Param::kSimulationTimeStep;
    volume_[kIdx] += delta;
    if (volume_[kIdx] < 5.2359877E-7) {
      volume_[kIdx] = 5.2359877E-7;
    }
    UpdateDiameter();
  }

  void UpdateDiameter() {
    // V = (4/3)*pi*r^3 = (pi/6)*diameter^3
    diameter_[kIdx] = std::cbrt(volume_[kIdx] * 6 / Math::kPi);
  }

  void UpdateVolume() {
    // V = (4/3)*pi*r^3 = (pi/6)*diameter^3
    volume_[kIdx] = Math::kPi / 6 * std::pow(diameter_[kIdx], 3);
  }

  void UpdateMassLocation(const array<float, 3>& delta) {
    mass_location_[kIdx][0] += delta[0];
    mass_location_[kIdx][1] += delta[1];
    mass_location_[kIdx][2] += delta[2];
  }

  void GetForceOn(const array<float, 3>& ref_mass_location,
                  float ref_diameter, array<float, 3>* force) const {
    DefaultForce default_force;
    float iof_coefficient = Param::kSphereDefaultInterObjectCoefficient;

    default_force.ForceBetweenSpheres(ref_mass_location, ref_diameter,
                                      iof_coefficient, mass_location_[kIdx],
                                      diameter_[kIdx], iof_coefficient, force);
  }

 protected:
  /// Returns the position in the polar coordinate system (cylindrical or
  /// spherical) of a point expressed in global cartesian coordinates
  /// ([1,0,0],[0,1,0],[0,0,1]).
  /// @param coord: position in absolute coordinates - [x,y,z] cartesian values
  /// @return the position in local coordinates
  array<float, 3> TransformCoordinatesGlobalToPolar(
      const array<float, 3>& coord) const;

  vec<array<float, 3>> position_;
  vec<array<float, 4>> gpu_position_;
  vec<array<float, 3>> mass_location_;
  vec<array<float, 3>> tractor_force_;
  vec<float> diameter_;
  vec<float> volume_;
  vec<float> adherence_;
  vec<float> density_;

  /// First axis of the local coordinate system.
  vec<array<float, 3>> x_axis_ = {array<float, 3>{1.0, 0.0, 0.0}};
  /// Second axis of the local coordinate system.
  vec<array<float, 3>> y_axis_ = {array<float, 3>{0.0, 1.0, 0.0}};
  /// Third axis of the local coordinate system.
  vec<array<float, 3>> z_axis_ = {array<float, 3>{0.0, 0.0, 1.0}};

  /// stores a list of neighbor ids for each scalar cell
  vec<InlineVector<int, 8>> neighbors_;

  /// collection of biology modules which define the internal behavior
  vec<vector<TBiologyModuleVariant>> biology_modules_;
};

template <typename Backend = Scalar,
          typename TBiologyModuleVariant = variant<NullBiologyModule>>
using Cell = CellExt<SimulationObject<Backend>, TBiologyModuleVariant>;

// ----------------------------------------------------------------------------
// Implementation -------------------------------------------------------------

template <typename T, typename U>
template <typename TBiologyModule>
inline void CellExt<T, U>::AddBiologyModule(TBiologyModule&& module) {
  biology_modules_[kIdx].emplace_back(module);
}

template <typename T, typename U>
inline void CellExt<T, U>::RunBiologyModules() {
  RunVisitor<Self<Backend>> visitor(this);
  for (auto& module : biology_modules_[kIdx]) {
    visit(visitor, module);
  }
}

template <typename T, typename U>
inline void CellExt<T, U>::Divide(Self<Scalar>* daughter) {
  Divide(daughter, 0.9 + 0.2 * gRandom.NextDouble());
}

template <typename T, typename U>
inline void CellExt<T, U>::Divide(Self<Scalar>* daughter, float volume_ratio) {
  // find random point on sphere (based on :
  // http://mathworld.wolfram.com/SpherePointPicking.html)
  float theta = 2 * Math::kPi * gRandom.NextDouble();
  float phi = std::acos(2 * gRandom.NextDouble() - 1);
  DivideImpl(daughter, volume_ratio, phi, theta);
}

template <typename T, typename U>
inline void CellExt<T, U>::Divide(Self<Scalar>* daughter,
                                  const array<float, 3>& axis) {
  auto polarcoord = TransformCoordinatesGlobalToPolar(
      Matrix::Add(axis, mass_location_[kIdx]));
  DivideImpl(daughter, 0.9 + 0.2 * gRandom.NextDouble(), polarcoord[1],
             polarcoord[2]);
}

template <typename T, typename U>
inline void CellExt<T, U>::Divide(Self<Scalar>* daughter, float volume_ratio,
                                  const array<float, 3>& axis) {
  auto polarcoord = TransformCoordinatesGlobalToPolar(
      Matrix::Add(axis, mass_location_[kIdx]));
  DivideImpl(daughter, volume_ratio, polarcoord[1], polarcoord[2]);
}

template <typename T, typename U>
inline void CellExt<T, U>::Divide(Self<Scalar>* daughter, float volume_ratio,
                                  float phi, float theta) {
  DivideImpl(daughter, volume_ratio, phi, theta);
}

template <typename T, typename TBiologyModuleVariant>
inline void CellExt<T, TBiologyModuleVariant>::DivideImpl(
    Self<Scalar>* daughter, float volume_ratio, float phi, float theta) {
  // A) Defining some values
  // ..................................................................
  // defining the two radii s.t total volume is conserved
  // * radius^3 = r1^3 + r2^3 ;
  // * volume_ratio = r2^3 / r1^3
  float radius = diameter_[kIdx] * 0.5;
  float r1 = radius / std::pow(1.0 + volume_ratio, 1.0 / 3.0);
  float r2 = radius / std::pow(1.0 + 1 / volume_ratio, 1.0 / 3.0);

  // define an axis for division (along which the nuclei will move)
  float x_coord = std::cos(theta) * std::sin(phi);
  float y_coord = std::sin(theta) * std::sin(phi);
  float z_coord = std::cos(phi);
  float total_length_of_displacement = radius / 4.0;
  array<float, 3> axis_of_division{
      total_length_of_displacement *
          (x_coord * x_axis_[kIdx][0] + y_coord * y_axis_[kIdx][0] +
           z_coord * z_axis_[kIdx][0]),
      total_length_of_displacement *
          (x_coord * x_axis_[kIdx][1] + y_coord * y_axis_[kIdx][1] +
           z_coord * z_axis_[kIdx][1]),
      total_length_of_displacement *
          (x_coord * x_axis_[kIdx][2] + y_coord * y_axis_[kIdx][2] +
           z_coord * z_axis_[kIdx][2])};

  // two equations for the center displacement :
  //  1) d2/d1= v2/v1 = volume_ratio (each sphere is shifted inver.
  //  proportionally to its volume)
  //  2) d1 + d2 = TOTAL_LENGTH_OF_DISPLACEMENT
  float d_2 = total_length_of_displacement / (volume_ratio + 1);
  float d_1 = total_length_of_displacement - d_2;

  // B) Instantiating a new sphere = 2nd daughter
  daughter->x_axis_[0] = x_axis_[kIdx];
  daughter->y_axis_[0] = y_axis_[kIdx];
  daughter->z_axis_[0] = z_axis_[kIdx];
  daughter->adherence_[0] = adherence_[kIdx];
  daughter->density_[0] = density_[kIdx];

  daughter->diameter_[0] = r2 * 2;
  daughter->UpdateVolume();

  // Mass Location
  array<float, 3> new_mass_location{
      mass_location_[kIdx][0] + d_2 * axis_of_division[0],
      mass_location_[kIdx][1] + d_2 * axis_of_division[1],
      mass_location_[kIdx][2] + d_2 * axis_of_division[2]};
  daughter->mass_location_[0] = new_mass_location;
  daughter->position_[0][0] = daughter->mass_location_[0][0];
  daughter->position_[0][1] = daughter->mass_location_[0][1];
  daughter->position_[0][2] = daughter->mass_location_[0][2];

  CopyVisitor<vector<TBiologyModuleVariant>> visitor(
      Event::kCellDivision, &(daughter->biology_modules_[0]));
  for (auto& module : biology_modules_[kIdx]) {
    visit(visitor, module);
  }

  // E) This sphere becomes the 1st daughter
  // move these cells on opposite direction
  position_[kIdx][0] -= d_1 * axis_of_division[0];
  position_[kIdx][1] -= d_1 * axis_of_division[1];
  position_[kIdx][2] -= d_1 * axis_of_division[2];
  mass_location_[kIdx][0] = position_[kIdx][0];
  mass_location_[kIdx][1] = position_[kIdx][1];
  mass_location_[kIdx][2] = position_[kIdx][2];

  // F) change properties of this cell
  diameter_[kIdx] = r1 * 2;
  UpdateVolume();

  // G) TODO(lukas) Copy the intracellular and membrane bound Substances
}

template <typename T, typename U>
array<float, 3> CellExt<T, U>::TransformCoordinatesGlobalToPolar(
    const array<float, 3>& pos) const {
  auto vector_to_point = Matrix::Subtract(pos, mass_location_[kIdx]);
  array<float, 3> local_cartesian{Matrix::Dot(x_axis_[kIdx], vector_to_point),
                                   Matrix::Dot(y_axis_[kIdx], vector_to_point),
                                   Matrix::Dot(z_axis_[kIdx], vector_to_point)};
  return {std::sqrt(local_cartesian[0] * local_cartesian[0] +
                    local_cartesian[1] * local_cartesian[1] +
                    local_cartesian[2] * local_cartesian[2]),
          std::atan2(std::sqrt(local_cartesian[0] * local_cartesian[0] +
                               local_cartesian[1] * local_cartesian[1]),
                     local_cartesian[2]),
          std::atan2(local_cartesian[1], local_cartesian[0])};
}

}  // namespace bdm

#endif  // CELL_H_
