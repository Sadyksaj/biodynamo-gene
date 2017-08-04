#ifndef DEMO_BACKUP_RESTORE_H_
#define DEMO_BACKUP_RESTORE_H_

#include <unistd.h>

#include "biology_module_util.h"
#include "cell.h"
#include "command_line_options.h"
#include "resource_manager.h"
#include "scheduler.h"
#include "simulation_object_util.h"

namespace bdm {

struct TestBehaviour {
  template <typename T>
  void Run(T* cell) {
    usleep(35000);  // 35 ms -> one iteration will take 350 ms
    cell->SetDiameter(cell->GetDiameter() + 1);
  }

  bool IsCopied(Event event) const { return true; }
  ClassDefNV(TestBehaviour, 1);
};

typedef Variant<TestBehaviour> BiologyModules;

template <typename Backend = Scalar>
using MyCell = Cell<Backend, BiologyModules>;

inline int Simulate(const CommandLineOptions& options) {
  Param::backup_every_x_seconds_ = 1;

  auto& cells = ResourceManager<MyCell<Soa>>::Get()->GetCells();
  for (size_t i = 0; i < 10; i++) {
    MyCell<> cell({100.0 * i, 100.0 * i, 100.0 * i});  // no colliding cells
    cell.SetDiameter(i);
    cell.AddBiologyModule(TestBehaviour());
    cells->push_back(cell);
  }

  Scheduler<MyCell<Soa>> scheduler(options.backup_file_, options.restore_file_);

  // will perform backup after iteration 3
  scheduler.Simulate(3);  // 1050 ms

  // application crash will happen inside this call
  scheduler.Simulate(11);  // 3850 ms

  // another call to Simulate after recovery
  scheduler.Simulate(2);

  // check result
  for (size_t i = 0; i < 10; i++) {
    if ((*cells)[i].GetDiameter() != 16 + i) {
      std::cerr << "Test failure: result incorrect" << std::endl;
      return 1;
    }
  }
  std::cout << "Test finished successfully" << std::endl;
  return 0;
}

}  // namespace bdm

#endif  // DEMO_BACKUP_RESTORE_H_