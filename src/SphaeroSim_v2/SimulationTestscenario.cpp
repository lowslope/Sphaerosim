#include "SimulationTestscenario.h"

#ifndef _INCLUDED_STRING_H_
#include <string>
#define _INCLUDED_STRING_H_
#endif

#ifndef _INCLUDED_VECTOR_H_
#include <vector>
#define _INCLUDED_VECTOR_H_
#endif

#ifndef _SPHAEROSIM_ELEMENT_H_
#include "Element.h"
#endif

#ifndef _SPHAEROSIM_INTERPOLATION_H_
#include "Interpolation.h"
#endif

#ifndef _SPHAEROSIM_BOUNDARY_CONDITION_H_
#include "BoundaryCondition.h"
#endif

#ifndef _SPHAEROSIM_GROWTH_MODEL_H_

#include "GrowthModel.h"

#endif

namespace SphaeroSim {

    SimulationTestscenario::
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
                           Material *material) :
            Simulation(name,
                       convergence_criterium,
                       cell_size_m,
                       temporal_interpolation,
                       spatial_interpolation,
                       output_folder,
                       create_new_subfolder,
                       vtk_result_option,
                       csv_separator,
                       noslip_neighbourhood,
                       calibration_file,
                       rng_state_registers,
                       integral_error,
                       save_memory,
                       calculate_flow_energy,
                       cache_convert_temperature,
                       cache_convert_velocity,
                       result_step,
                       statistics_step,
                       numerical_options,
                       growth_model,
                       nucleation_model,
                       material) {

        // boundary conditions for the testscenario
        timesteps_ = timesteps;
        temperatures_ = temperatures;
        shearrates_ = shearrates;
        strainrates_ = strainrates;
        velocity_profile_ = velocity_profile;
        temperature_profile_ = temperature_profile;
        num_cells_ = num_cells;

        EnsureConsistentDatapoints();

        SetNumRuns(timesteps_.size());
    }

    void SimulationTestscenario::EnsureConsistentDatapoints() {
        // ensure the same arraysize
        std::vector<double> *value_lists[] = {&timesteps_,
                                              &result_timesteps_,
                                              &shearrates_,
                                              &strainrates_,
                                              &temperatures_};
        const std::size_t num_entries =
                sizeof(value_lists) / sizeof(std::vector<double> *);

        // find the maximum size
        std::size_t maximum_index = 0;
        for (std::size_t cur = 1; cur < num_entries; ++cur) {
            if (value_lists[maximum_index]->size() < value_lists[cur]->size()) {
                maximum_index = cur;
            }
        }

        // add values to the lists that do not have enough points
        for (std::size_t cur = 0; cur < num_entries; ++cur) {
            std::vector<double> *dataset = value_lists[cur];

            // if no value in the list, assume the list takes the values '0.0'
            if (dataset->size() == 0) {
                dataset->push_back(0.0);
            }

            std::size_t cur_array_index = 0;
            std::size_t original_size = dataset->size();
            while (dataset->size() < value_lists[maximum_index]->size()) {
                dataset->push_back(dataset->at(cur_array_index));
                ++cur_array_index;
                if (cur_array_index >= original_size) {
                    cur_array_index = 0;
                }
            }
        }

        for (std::size_t cur = 0; cur < timesteps_.size(); ++cur) {
            if (result_timesteps_[cur] < timesteps_[cur]) {
                result_timesteps_[cur] = timesteps_[cur];
            }
        }
    }

    void SimulationTestscenario::
    SetElementNeighbourhood(std::map<const uint64_t, Element *> *dst) const {
        std::vector<Element::neighbour_info> infos;
        for (std::size_t cur = 0; cur < 6; ++cur) {
            if (GetNoslipNeighbourhood(cur) == true) {
                infos.push_back(Element::neighbour_info(std::numeric_limits<uint64_t>::max(),
                                                        false));
            } else {
                infos.push_back(Element::neighbour_info(cur + 1, false));
            }
        }

        std::map<const uint64_t, Element *>::iterator it = dst->begin();
        while (it != dst->end()) {
            (*it).second->SetNeighbourElementInfos(infos, *dst);
            ++it;
        }
    }

    void SimulationTestscenario::BuildElements(std::vector<Element *> *dst,
                                               RandomNumberGenerator *rng) {
        Element *new_element = new Element(0,
                                           Eigen::Vector3d(0.0, 0.0, 0.0),
                                           GetCellSizeM(),
                                           num_cells_,
                                           GetTemporalInterpolationType(),
                                           GetSpatialInterpolationType(),
                                           GetGrowthModel()->GetGrowthType(),
                                           rng,
                                           this);

        std::vector<Element::temperature_point> temperatures;
        std::vector<Element::velocity_point> velocities;

        // Boundary conditions
        // starting temperatures
        std::vector<double> initial_temperatures;
        initial_temperatures.resize(8);
        for (std::size_t cur = 0; cur < 8; ++cur)
            initial_temperatures[cur] = temperatures_[GetCurrentRunIndex()];

        // initial temperature
        temperatures.push_back(Element::temperature_point(0.0,
                                                          initial_temperatures));

        // start velocities
        std::vector<Eigen::Vector3d> initial_velocity_field;
        initial_velocity_field.resize(8);
        for (std::size_t cur = 0; cur < 8; ++cur) {
            double shearrate = shearrates_[GetCurrentRunIndex()];
            double strainrate = strainrates_[GetCurrentRunIndex()];

            // shear flow: vx = 0, vy = xpos*shearrate, vz = 0
            Eigen::Vector3d cur_vel_shear(0.0, 0.0, 0.0);
            cur_vel_shear[1] = shearrate * new_element->GetCornerPosition(cur)[0];

            // elongation flow
            // vx = strainrate*x, vy = -0.5*strainrate*y, vz = -0.5*strainrate*z
            Eigen::Vector3d cur_vel_elongational(0.0, 0.0, 0.0);
            cur_vel_elongational[0] =
                    strainrate * new_element->GetCornerPosition(cur)[0];
            cur_vel_elongational[1] =
                    -0.5 * strainrate * new_element->GetCornerPosition(cur)[1];
            cur_vel_elongational[2] =
                    -0.5 * strainrate * new_element->GetCornerPosition(cur)[2];

            initial_velocity_field[cur] = cur_vel_shear + cur_vel_elongational;
        }
        velocities.push_back(Element::velocity_point(0.0, initial_velocity_field));

        // temperature profile
        std::vector<double> temperature_field;
        temperature_field.resize(8);
        std::size_t num_points = temperature_profile_.size();
        for (std::size_t cur_point = 0; cur_point < num_points; ++cur_point) {
            temperature_field = initial_temperatures;
            for (std::size_t cur = 0; cur < temperature_field.size(); ++cur)
                temperature_field[cur] *= temperature_profile_[cur_point][1];
            temperatures.push_back(
                    Element::temperature_point(temperature_profile_[cur_point][0],  // point in time
                                               temperature_field));
        }
        // velocity profile
        std::vector<Eigen::Vector3d> velocity_field;
        velocity_field.resize(8);
        num_points = velocity_profile_.size();
        for (std::size_t cur_point = 0; cur_point < num_points; ++cur_point) {
            velocity_field = initial_velocity_field;
            for (std::size_t cur = 0; cur < velocity_field.size(); ++cur)
                velocity_field[cur] *= velocity_profile_[cur_point][1];
            velocities.push_back(
                    Element::velocity_point(velocity_profile_[cur_point][0],  // point in time
                                            velocity_field));  // velocitiy field
        }

        // only one element in a testscenario
        new_element->SetTemperatureField(temperatures);
        new_element->SetVelocityField(velocities);

        dst->push_back(new_element);
    }

    const double SimulationTestscenario::Timestep() {
        return timesteps_[GetCurrentRunIndex()];
    }
}  // namespace SphaeroSim