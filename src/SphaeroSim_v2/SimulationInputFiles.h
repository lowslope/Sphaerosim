#ifndef _SPHAEROSIM_SIMULATION_INPUT_FILES_H_
#define _SPHAEROSIM_SIMULATION_INPUT_FILES_H_

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

    class SimulationInputFiles : public Simulation {
    public:
        SimulationInputFiles(const double timestep,
                             const std::string &name,
                             const ResultStep convergence_criterium,
                             const double cell_size_m,
                             const InterpolationType temporal_interpolation_type,
                             const InterpolationType spatial_interpolation_type,
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

        virtual ~SimulationInputFiles();

        // interface
        virtual void BuildElements(std::vector<Element *> *dst,
                                   RandomNumberGenerator *rng) = 0;

        virtual const double Timestep();

        virtual void
        SetElementNeighbourhood(std::map<const uint64_t, Element *> *dst) const;

    protected:
        const double timestep_;

        struct GeometryInformation {
            std::vector<Eigen::Vector3d> locations;
            std::vector<Eigen::Vector3d> sizes;
            std::vector<std::vector<std::size_t> > indices;
        };
        GeometryInformation geometry_information_;

    private:
        void ElementSharingPoints(const std::size_t neighbour,
                                  std::vector<std::size_t> *element_points,
                                  std::vector<std::size_t> *neighbour_points) const;

        bool
        CheckSameIndices(const uint64_t first_element_id,
                         const std::vector<std::size_t> &first_checkpoints,
                         const uint64_t second_element_id,
                         const std::vector<std::size_t> &second_checkpoints) const;

        bool
        CheckSameCornerPoints(const uint64_t first_element_id,
                              const std::vector<std::size_t> &first_checkpoints,
                              const uint64_t second_element_id,
                              const std::vector<std::size_t> &second_checkpoints) const;
    };

}  // namespace SphaeroSim

#endif  // _SPHAEROSIM_SIMULATION_TESTSCENARIO_H_
