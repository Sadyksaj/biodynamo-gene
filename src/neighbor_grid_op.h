#ifndef NEIGHBOR_GRID_OP_H_
#define NEIGHBOR_GRID_OP_H_

#include <utility>
#include <vector>

#include "grid.h"
#include "inline_vector.h"

namespace bdm {

/// A class that sets up an uniform grid to perform operations that require
/// knowledge about neighboring simulation objects
class NeighborGridOp {
 public:
  explicit NeighborGridOp(Grid::Adjacency adjacency = Grid::kHigh,
                          bool set_local_neighbors = false,
                          float radius = 3000)
      : adjacency_(adjacency),
        set_local_neighbors_(set_local_neighbors),
        radius_(radius) {}
  virtual ~NeighborGridOp() {}

  template <typename TContainer>
  void Compute(TContainer* cells) const {
    // Construct a 3D grid with the current positions for the simulation objects
    auto& grid = Grid::GetInstance();
    grid.Initialize(cells, adjacency_);

    if (set_local_neighbors_) {
      // Initiate the operation
      grid.SetNeighborsWithinRadius(cells, radius_);
    }
  }

 private:
  /// Determines how many neighboring boxes to consider for neighbor operations
  Grid::Adjacency adjacency_;
  /// Boolean to cache the local neighbors for each simulation object or not
  bool set_local_neighbors_;
  /// The searching radius for which to set the local neighbors to
  float radius_;
};

}  // namespace bdm

#endif  // NEIGHBOR_GRID_OP_H_
