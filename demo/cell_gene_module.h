#ifndef DEMO_CELL_GENE_MODULE_H_
#define DEMO_CELL_GENE_MODULE_H_

#include "biodynamo.h"

namespace bdm {

  using std::array;
  using std::vector;
  using std::string;
  // using bdm::Param::step_global_;
  //This constant specifies the amount of different proteins in the simulation
  const int protein_amount = 3;
  //This constant specifies method using for solving differential equation {"Euler", "RK4"}.
  const string DE_solve_method = "RK4";

  BDM_SIM_OBJECT(GeneCell, Cell) {
    BDM_SIM_OBJECT_HEADER(GeneCellExt, 1, substances_);

   public:
    vec<array<double, protein_amount>> substances_;
    GeneCellExt() {}
    explicit GeneCellExt(const array<double, 3>& position, array<double, protein_amount> init_substances)
     : Base(position), substances_(init_substances){}

    virtual ~GeneCellExt() {}

    // In this method described functions. These functions determine changes in concentration of proteins.
    // You should describe your own function and push it to the update_value vector.
    // Value of protein concentration inside substances_ vector.
    // update_value.push_back(<"your function">);
    vector<double> evolve(double curr_time){
      vector<double> update_value;
      update_value.push_back(1/substances_[0][0]);
      update_value.push_back(1/substances_[0][1]);
      update_value.push_back(1/substances_[0][2]);
      assert(update_value.size() == protein_amount && "Amount of functions does not equal to amount of proteins\n");
      return update_value;
    }

    // void DivideImpl(typename Base::template Self<Scalar> * daughter, double volume_ratio, double phi, double theta) override {
    //   std::cout<<"I'm out\n";
    //   // daughter->substances_ = substances_ * 0.5;
    //   // substances_ = substances_ * 0.5;
    //   Base::DivideImpl(daughter, volume_ratio, phi, theta);
    // }
    void DivideImpl(void* daughter, double volume_ratio, double phi, double theta)
      override {
        auto daughter_cast = static_cast<Self<Scalar>*>(daughter);
        for (int j =  0; j < protein_amount; j++){
          daughter_cast->substances_[0][j] = substances_[0][j] * 0.4;
          substances_[0][j] = substances_[0][j] * 0.6;
        }
        // forward call to implementation in CellExt
        Base::DivideImpl(daughter_cast, volume_ratio, phi, theta);
      }
  };

  struct GeneCalculation : public BaseBiologyModule {
    double time_step = Param::simulation_time_step_;
    // size_t& step = Param::step_global_;

    GeneCalculation() : BaseBiologyModule(gAllBmEvents){}

    template <typename T>
    void Run(T* cell) {
      std::cout<<"evolved\n";
      if (DE_solve_method == "Euler"){
        // evolve method needs in current time of simulation
        vector<double> update_value = cell->evolve(Param::step_global_ * time_step);
        std::cout<<cell->substances_[0][1] << " - evolved to 1\n";
        for (int i = 0; i < protein_amount; i++){
          cell->substances_[0][i] += update_value[i] * time_step;
        }
      }
      else if (DE_solve_method == "RK4"){
        std::cout<<cell->substances_[0][1] << " - evolved to 1\n";
        vector<double> k1 = cell->evolve(Param::step_global_ * time_step);
        for (int i = 0; i < protein_amount; i++)
          cell->substances_[0][i] += time_step*k1[i]/2.0f;
        vector<double> k2 = cell->evolve(Param::step_global_ * time_step + time_step/2.0f);
        for (int i = 0; i < protein_amount; i++)
          cell->substances_[0][i] += time_step*k2[i]/2.0f - time_step*k1[i]/2.0f;
        vector<double> k3 = cell->evolve(Param::step_global_ * time_step + time_step/2.0f);
        for (int i = 0; i < protein_amount; i++)
          cell->substances_[0][i] += time_step*k3[i] - time_step*k2[i]/2.0f;
        vector<double> k4 = cell->evolve(Param::step_global_ * time_step + time_step);
        for (int i = 0; i < protein_amount; i++)
          cell->substances_[0][i] += time_step*(k1[i] + 2*k2[i] + 2*k3[i] + k4[i])/6.0f;
      }
      cell->ChangeVolume(300);
      Divide(*cell);
    }

    ClassDefNV(GeneCalculation, 1);
  };

// 2. Define compile time parameter
template <typename Backend>
struct CompileTimeParam : public DefaultCompileTimeParam<Backend> {
  using BiologyModules = Variant<GeneCalculation>;
  using AtomicTypes = VariadicTypedef<GeneCell>;
};

inline int Simulate(int argc, const char** argv) {
  // 2. Initialize BioDynaMo
  InitializeBioDynamo(argc, argv);
  size_t cells_per_dim = 1;
  auto construct = [](const std::array<double, 3>& position) {
    array<double, protein_amount> init_vals = {10.0, 10.0, 10.0};
    GeneCell cell(position, init_vals);
    cell.SetDiameter(30);
    cell.SetAdherence(0.4);
    cell.SetMass(1.0);
    cell.AddBiologyModule(GeneCalculation());
    return cell;
  };
  ModelInitializer::Grid3D(cells_per_dim, 20, construct);

    Scheduler<> scheduler;
    scheduler.Simulate(2);
  return 0;
}

}  // namespace bdm

#endif  // DEMO_CELL_GENE_MODULE_H_
