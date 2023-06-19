#ifndef _SPHAEROSIM_SIMULATION_TESTSCENARIO_H_
#define _SPHAEROSIM_SIMULATION_TESTSCENARIO_H_

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

namespace SphaeroSim {

    class SimulationTestscenario : public Simulation {
    public:
        SimulationTestscenario(const std::vector<double> &temperatures,
                               const std::vector<double> &shearrates,
                               const std::vector<double> &strainrates,
                               const std::vector<Eigen::Vector2d> &temperature_profile,
                               const std::vector<Eigen::Vector2d> &velocity_profile,
                               const Eigen::Vector3s num_cells,
                               const std::string &name,
                               const ResultStep convergence_criterium,
                               const std::vector<double> &timesteps,
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

        void SetElementNeighbourhood(std::map<const uint64_t, Element *> *dst) const;

        const double Timestep();

    private:
        void EnsureConsistentDatapoints();

        std::vector<double> timesteps_;
        std::vector<double> result_timesteps_;
        std::vector<double> temperatures_;
        std::vector<double> shearrates_;
        std::vector<double> strainrates_;
        std::vector<Eigen::Vector2d> temperature_profile_;
        std::vector<Eigen::Vector2d> velocity_profile_;
        Eigen::Vector3s num_cells_;
    };

}  // namespace SphaeroSim

#endif  // _SPHAEROSIM_SIMULATION_TESTSCENARIO_H_
