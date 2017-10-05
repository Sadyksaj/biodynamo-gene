#ifndef VISUALIZATION_CATALYST_ADAPTOR_H_
#define VISUALIZATION_CATALYST_ADAPTOR_H_

#include <TError.h>
#include <string>

// check for ROOTCLING was necessary, due to ambigous reference to namespace
// detail when using ROOT I/O
#if defined(USE_CATALYST) && !defined(__ROOTCLING__)
#include <vtkCPDataDescription.h>
#include <vtkCPInputDataDescription.h>
#include <vtkCPProcessor.h>
#include <vtkCPPythonScriptPipeline.h>
#include <vtkDoubleArray.h>
#include <vtkFieldData.h>
#include <vtkIdTypeArray.h>
#include <vtkNew.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkStringArray.h>
#include <vtkUnstructuredGrid.h>

#endif  // defined(USE_CATALYST) && !defined(__ROOTCLING__)

namespace bdm {

#if defined(USE_CATALYST) && !defined(__ROOTCLING__)

/// The class that bridges the simulation code with ParaView
class CatalystAdaptor {
 public:
  static CatalystAdaptor* GetInstance() {
    static CatalystAdaptor kInstance;
    return &kInstance;
  }

  ///
  /// @brief      Builds a VTK grid.
  ///
  /// @param      sim_objects  The simulation objects
  ///
  /// @tparam     TContainer   { The container that holds the simulation objects
  ///                          }
  ///
  template <typename TContainer>
  void BuildVTKGrid(TContainer* sim_objects) {
    // // create the points information
    // vtkNew<vtkDoubleArray> position_array;
    // vtkDoubleArray* diameter_array = vtkDoubleArray::New();
    // diameter_array->SetName("Diameter");
    //
    // position_array->SetNumberOfComponents(3);
    // position_array->SetArray(sim_objects->GetPositionPtr(),
    //                          static_cast<vtkIdType>(sim_objects->size() * 3),
    //                          1);
    // diameter_array->SetArray(sim_objects->GetDiameterPtr(),
    //                          static_cast<vtkIdType>(sim_objects->size()), 1);
    //
    // vtkNew<vtkPoints> points;
    //
    // points->SetData(position_array.GetPointer());
    //
    // g_vtk_grid_->SetPoints(points.GetPointer());
    // g_vtk_grid_->GetPointData()->AddArray(diameter_array);
  }

  ///
  /// @brief      Wrapper around @ref BuildVTKGRid to define
  ///
  /// @param      sim_objects  The simulation objects
  ///
  /// @tparam     TContainer   { The container that holds the simulation objects
  ///                          }
  ///
  template <typename TContainer>
  void BuildVTKDataStructures(TContainer* sim_objects) {
    if (g_vtk_grid_ == NULL) {
      // The grid structure isn't changing so we only build it
      // the first time it's needed. If we needed the memory
      // we could delete it and rebuild as necessary.
      g_vtk_grid_ = vtkUnstructuredGrid::New();
    }
    BuildVTKGrid(sim_objects);
  }

  ///
  /// @brief      Initializes Catalyst with the predefined pipeline
  ///
  /// @param[in]  script  The Python script that contains the pipeline
  ///
  inline void Initialize(const std::string& script) {
    if (g_processor_ == NULL) {
      g_processor_ = vtkCPProcessor::New();
      g_processor_->Initialize();
    } else {
      g_processor_->RemoveAllPipelines();
    }
    vtkNew<vtkCPPythonScriptPipeline> pipeline;
    pipeline->Initialize(script.c_str());
    g_processor_->AddPipeline(pipeline.GetPointer());
  }

  /// Cleans up allocated memory
  inline void Finalize() {
    if (g_processor_) {
      // Call made to MPI_FINALIZE
      g_processor_->Delete();
      g_processor_ = NULL;
    }
    if (g_vtk_grid_) {
      g_vtk_grid_->Delete();
      g_vtk_grid_ = NULL;
    }
  }

  ///
  /// @brief      Applies the pipeline to the data structures in the VTK Grid
  ///
  /// @param      sim_objects     The simulation objects
  /// @param[in]  step            The simulation step
  /// @param[in]  time_step       The time step duration
  /// @param[in]  last_time_step  Last time step or not
  ///
  /// @tparam     TContainer      { The container that holds the simulation
  ///                             objects }
  ///
  template <typename TContainer>
  inline void CoProcess(TContainer* sim_objects, double step, size_t time_step,
                        bool last_time_step) {
    vtkNew<vtkCPDataDescription> data_description;
    data_description->AddInput("input");
    data_description->SetTimeData(step, time_step);
    if (last_time_step == true) {
      // assume that we want to all the pipelines to execute if it
      // is the last time step.
      data_description->ForceOutputOn();
    }

    // If we segfault at here it probably means that the pipeline was not
    // initialized (with a python script)
    if (g_processor_->RequestDataDescription(data_description.GetPointer()) !=
        0) {
      BuildVTKDataStructures(sim_objects);
      data_description->GetInputDescriptionByName("input")->SetGrid(
          g_vtk_grid_);

      g_processor_->CoProcess(data_description.GetPointer());
    }

    // // ----------------- User changes propagation
    // ------------------------------

    // vtkFieldData* user_data = data_description->GetUserData();
    // if (!user_data) {
    //   // no user changes
    //   return;
    // }

    // // Which properties/attribute the user changed
    // vtkStringArray* prop_arrays =
    //     vtkStringArray::SafeDownCast(user_data->GetAbstractArray("PropArrays"));
    // if (!prop_arrays) {
    //   std::cout << "Warning: Cannot find propagated array names" << endl;
    //   return;
    // }

    // // Get every changed attribute
    // vtkIdTypeArray* idx_array;
    // vtkDoubleArray* val_array;
    // for (int j = 0; j < prop_arrays->GetSize(); j++) {
    //   auto attribute = prop_arrays->GetValue(j);
    //   idx_array = vtkIdTypeArray::SafeDownCast(user_data->GetAbstractArray(
    //       (std::string("PropIdx") + std::string(attribute)).c_str()));
    //   val_array = vtkDoubleArray::SafeDownCast(user_data->GetAbstractArray(
    //       (std::string("PropVals") + std::string(attribute)).c_str()));

    //   if (!idx_array || !val_array) {
    //     std::cerr << "Warning: null pointer returned while fetching '"
    //               << attribute << "' array " << endl;
    //   }

    //   // Update changed sim_objects
    //   for (int i = 0; i < idx_array->GetNumberOfTuples(); i++) {
    //     std::cout << "sim_objects[" << idx_array->GetValue(i)
    //               << "] = " << val_array->GetValue(i) << endl;

    //     // reflection here!
    //     (*sim_objects)[idx_array->GetValue(i)].SetDiameter(
    //         val_array->GetValue(i));
    //   }
    // }
  }

 private:
  vtkCPProcessor* g_processor_ = nullptr;
  vtkUnstructuredGrid* g_vtk_grid_;
};

#else

/// False front (to ignore Catalyst in gtests)
class CatalystAdaptor {
 public:
  static CatalystAdaptor* GetInstance() {
    static CatalystAdaptor kInstance;
    return &kInstance;
  }

  template <typename TContainer>
  void BuildVTKGrid(TContainer* sim_objects) {
    Fatal("",
          "Simulation was compiled without ParaView support, but you are "
          "trying to use it.");
  }

  template <typename TContainer>
  void BuildVTKDataStructures(TContainer* sim_objects) {
    Fatal("",
          "Simulation was compiled without ParaView support, but you are "
          "trying to use it.");
  }

  void Initialize(const std::string& script) {
    Fatal("",
          "Simulation was compiled without ParaView support, but you are "
          "trying to use it.");
  }

  void Finalize() {
    Fatal("",
          "Simulation was compiled without ParaView support, but you are "
          "trying to use it.");
  }

  template <typename TContainer>
  void CoProcess(TContainer* sim_objects, double time, size_t time_step,
                 bool last_time_step) {
    Fatal("",
          "Simulation was compiled without ParaView support, but you are "
          "trying to use it.");
  }
};

#endif  // defined(USE_CATALYST) && !defined(__ROOTCLING__)

}  // namespace bdm

#endif  // VISUALIZATION_CATALYST_ADAPTOR_H_
