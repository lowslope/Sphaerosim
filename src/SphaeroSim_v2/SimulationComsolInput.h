#ifndef _SPHAEROSIM_SIMULATION_COMSOL_INPUT_H_
#define _SPHAEROSIM_SIMULATION_COMSOL_INPUT_H_

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

    class SimulationComsolInput : public SimulationInputFiles {
    public:
        SimulationComsolInput(const double timestep,
                              const std::vector<std::string> &file_list,
                              const std::string &coordinate_x_identifier,
                              const std::string &coordinate_y_identifier,
                              const std::string &coordinate_z_identifier,
                              const std::string &velocity_x_identifier,
                              const std::string &velocity_y_identifier,
                              const std::string &velocity_z_identifier,
                              const std::string &temperature_identifier,
                              const std::string &levelset_identifier,
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

        const std::string coordinate_x_identifier_;
        const std::string coordinate_y_identifier_;
        const std::string coordinate_z_identifier_;
        const std::string velocity_x_identifier_;
        const std::string velocity_y_identifier_;
        const std::string velocity_z_identifier_;
        const std::string temperature_identifier_;
        const std::string levelset_identifier_;

        std::vector<std::string> file_list_;
    };

}  // namespace SphaeroSim

#endif  // _SPHAEROSIM_SIMULATION_COMSOL_INPUT_H_
