#ifndef _SPHAEROSIM_SIMULATION_VTK_INPUT_H_
#define _SPHAEROSIM_SIMULATION_VTK_INPUT_H_

#ifndef _INCLUDED_STRING_H_

#include <string>

#define _INCLUDED_STRING_H_
#endif

#ifndef _INCLUDED_VECTOR_H_

#include <vector>

#define _INCLUDED_VECTOR_H_
#endif

#ifndef _SPHAEROSIM_EIGEN_LIBRARY_H_

#include "EigenLibrary.h"

#endif

#ifndef _SPHAEROSIM_SIMULATION_H_

#include "Simulation.h"

#endif

#ifndef _SPHAEROSIM_SIMULATION_INPUT_FILES_H_

#include "SimulationInputFiles.h"

#endif

namespace SphaeroSim {

    class SimulationVTKInput : public SimulationInputFiles {
    public:
        SimulationVTKInput(const double timestep,
                           const std::vector<std::string> &folder_list,
                           const std::vector<double> &folder_timesteps,
                           const std::vector<bool> &folder_start_time_zero,
                           const std::vector<std::string> &folder_temperature_indicator,
                           const std::vector<std::string> &folder_velocity_indicator,
                           const std::vector<uint64_t> &element_ids,
                           const std::string &name,
                           const ResultStep convergence_criterium,
                           const double cell_size_m,
                           const InterpolationType temporal_interpolation,
                           const InterpolationType spatial_interpolation,
                           const std::string &output_folder,
                           const bool create_new_subfolder,
                           const VTKResultOption &vtk_result_option,
                           const char csv_separator,
                           const std::vector<bool> &noslip_neighbourhood,
                           const std::string &calibration_file,
                           const std::vector<uint32_t> &rng_state_registers,
                           const double integral_error,
                           const bool save_memory,
                           const bool calculate_flow_energy,
                           const double cache_convert_temperature,
                           const double cache_convert_velocity,
                           const ResultStep &result_step,
                           const ResultStep &statistics_step,
                           const NumericalOptions &numerical_options,
                           GrowthModel *growth_model,
                           NucleationModel *nucleation_model,
                           Material *material);

        // interface
        void BuildElements(std::vector<Element *> *dst,
                           RandomNumberGenerator *rng);

    private:
        void FilenameArray(const std::string &folderpath,
                           std::vector<std::string> *dst) const;

        const std::vector<double> folder_timesteps_;
        const std::vector<bool> folder_start_time_zero_;
        const std::vector<uint64_t> element_ids_;
        const std::vector<std::string> folder_temperature_indicator_;
        const std::vector<std::string> folder_velocity_indicator_;

        std::vector<std::string> folder_list_;
    };

}  // namespace SphaeroSim

#endif  // _SPHAEROSIM_SIMULATION_VTK_INPUT_H_
