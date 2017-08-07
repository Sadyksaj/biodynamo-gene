#ifndef SCHEDULER_H_
#define SCHEDULER_H_

#include "biology_module_op.h"
#include "displacement_op.h"
#include "neighbor_grid_op.h"
#include "neighbor_gpu_grid_op.h"
#include "neighbor_nanoflann_op.h"
#include "op_timer.h"
#include "resource_manager.h"

namespace bdm {

class Scheduler {
 public:
  Scheduler() {}

  virtual ~Scheduler() {}

  template <typename TCellContainer>
  void Simulate(unsigned steps) {
    OpTimer<NeighborGpuGridOp> neighbor("neighbor", NeighborGpuGridOp());
    OpTimer<BiologyModuleOp> biology("biology");
    OpTimer<DisplacementOp> physics("physics");

    auto cells = ResourceManager<TCellContainer>::Get()->GetCells();

    while (steps-- > 0) {
      neighbor.Compute(cells);
      biology.Compute(cells);
      physics.Compute(cells);
      cells->Commit();
    }
  }
};

}  // namespace bdm

#endif  // SCHEDULER_H_
