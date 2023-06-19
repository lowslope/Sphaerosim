#include "SimulationInputFiles.h"

#ifndef _INCLUDED_STDINT_H_
#include <stdint.h>
#define _INCLUDED_STDINT_H_
#endif

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

#ifndef _SPHAEROSIM_VTKFILE_H_

#include "VTKFile.h"

#endif

#ifndef __DIRECTORY_H

#include "directory.h"

#endif

#ifndef _SPHAEROSIM_ELEMENT_H_
#include "Element.h"
#endif

#ifndef _SPHAEROSIM_GROWTH_MODEL_H_

#include "GrowthModel.h"

#endif

namespace SphaeroSim {
    SimulationInputFiles::
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
                         Material *material) :
            Simulation(name,
                       convergence_criterium,
                       cell_size_m,
                       temporal_interpolation_type,
                       spatial_interpolation_type,
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
                       material),
            timestep_(timestep) {
    }

    SimulationInputFiles::~SimulationInputFiles() {

    }

//                       6--------7
//                      /|       /|                  
//                     2--------3 |                  y
//                     | |      | |                  |  z
//                     | 4------|-5                  | /
//                     |/       |/                   |/
//                     0--------1                    ------x
    void SimulationInputFiles::
    ElementSharingPoints(const std::size_t neighbour,
                         std::vector<std::size_t> *element_points,
                         std::vector<std::size_t> *neighbour_points) const {
        // the neighbour elements (+/- X_AXIS; +/- Y_AXIS; +/- Z_AXIS)
        if (neighbour == 0 || neighbour == 1) {
            std::size_t indices_element[] = {1, 3, 5, 7};
            std::size_t indices_neighbour[] = {0, 2, 4, 6};
            for (std::size_t cur = 0; cur < 4; ++cur) {
                element_points->push_back(neighbour == 0 ?
                                          indices_element[cur] :
                                          indices_neighbour[cur]);
                neighbour_points->push_back(neighbour == 0 ?
                                            indices_neighbour[cur] :
                                            indices_element[cur]);
            }
        } else if (neighbour == 2 || neighbour == 3) {
            std::size_t indices_element[] = {2, 3, 6, 7};
            std::size_t indices_neighbour[] = {0, 1, 4, 5};
            for (std::size_t cur = 0; cur < 4; ++cur) {
                element_points->push_back(neighbour == 2 ?
                                          indices_element[cur] :
                                          indices_neighbour[cur]);
                neighbour_points->push_back(neighbour == 2 ?
                                            indices_neighbour[cur] :
                                            indices_element[cur]);
            }
        } else if (neighbour == 4 || neighbour == 5) {
            std::size_t indices_element[] = {4, 5, 6, 7};
            std::size_t indices_neighbour[] = {0, 1, 2, 3};
            for (std::size_t cur = 0; cur < 4; ++cur) {
                element_points->push_back(neighbour == 4 ?
                                          indices_element[cur] :
                                          indices_neighbour[cur]);
                neighbour_points->push_back(neighbour == 4 ?
                                            indices_neighbour[cur] :
                                            indices_element[cur]);
            }
        }
    }

    bool SimulationInputFiles::
    CheckSameIndices(const uint64_t first_element_id,
                     const std::vector<std::size_t> &first_checkpoints,
                     const uint64_t second_element_id,
                     const std::vector<std::size_t> &second_checkpoints) const {
        for (std::size_t cur = 0; cur < first_checkpoints.size(); ++cur) {
            const std::size_t cur_index_1 =
                    geometry_information_.indices[first_element_id][first_checkpoints[cur]];
            const std::size_t cur_index_2 =
                    geometry_information_.indices[second_element_id][second_checkpoints[cur]];
            if (cur_index_1 != cur_index_2) {
                return false;
            }
        }
        return true;
    }

    bool SimulationInputFiles::
    CheckSameCornerPoints(const uint64_t first_element_id,
                          const std::vector<std::size_t> &first_checkpoints,
                          const uint64_t second_element_id,
                          const std::vector<std::size_t> &second_checkpoints) const {
        const double delta = 1e-8;
        for (std::size_t cur = 0; cur < first_checkpoints.size(); ++cur) {
            Eigen::Vector3d point_1 =
                    GetElement(first_element_id)->GetCornerPosition(first_checkpoints[cur]);
            Eigen::Vector3d point_2 =
                    GetElement(second_element_id)->GetCornerPosition(second_checkpoints[cur]);

            for (std::size_t dim = 0; dim < 3; ++dim) {
                if (fabs(point_1[dim] - point_2[dim]) > delta) {
                    return false;
                }
            }
        }
        return true;
    }

    void SimulationInputFiles::
    SetElementNeighbourhood(std::map<const uint64_t, Element *> *dst) const {
        for (auto it:*dst) {
            std::vector<Element::neighbour_info> infos;
            // the neighbour elements (+/- X_AXIS; +/- Y_AXIS; +/- Z_AXIS)
            for (std::size_t neighbour = 0; neighbour < 6; ++neighbour) {
                std::vector<std::size_t> checkpoints_element;
                std::vector<std::size_t> checkpoints_neighbour;

                ElementSharingPoints(neighbour,
                                     &checkpoints_element,
                                     &checkpoints_neighbour);

                bool found_neighbour = false;
                // first check the neighbours for which we have data
                for (auto check_it:*dst) {
                    if (it.first != check_it.first) {
                        if (CheckSameCornerPoints(it.first,
                                                  checkpoints_element,
                                                  check_it.first,
                                                  checkpoints_neighbour)) {
                            infos.push_back(Element::neighbour_info(check_it.first,  // id of the neighbour
                                                                    true));  // data exists
                            found_neighbour = true;
                            break;
                        }
                        if (CheckSameIndices(it.first,
                                             checkpoints_element,
                                             check_it.first,
                                             checkpoints_neighbour)) {
                            infos.push_back(Element::neighbour_info(check_it.first,  // id of the neighbour
                                                                    false));  // no direct data at the corner points exists
                            found_neighbour = true;
                            break;
                        }
                    }
                }

                if (!found_neighbour) {
                    // check all the elements from the vtk file if one is the neighbouring
                    // element
                    const std::size_t num_elements = geometry_information_.indices.size();
                    for (std::size_t ele_id = 0; ele_id < num_elements; ++ele_id) {
                        if (CheckSameIndices(it.first,
                                             checkpoints_element,
                                             ele_id,
                                             checkpoints_neighbour)) {
                            infos.push_back(Element::neighbour_info(ele_id,  // id of the neighbour
                                                                    false));  // no direct data at the corner points exists
                            found_neighbour = true;
                            break;
                        }
                    }
                }
                const uint64_t no_element = UINT64_MAX;
                if (!found_neighbour) { // no neighbour exists
                    if (GetNoslipNeighbourhood(neighbour)) {
                        infos.push_back(Element::neighbour_info(
                                no_element,  // no slip conditions
                                false));  // no direct data at the corner points exists
                    } else {
                        infos.push_back(Element::neighbour_info(
                                no_element - 1,  // dummy id
                                false));
                    }  // no direct data at the corner points exists
                }
            }
            it.second->SetNeighbourElementInfos(infos, *dst);

        }
    }

    const double SimulationInputFiles::Timestep() {
        return timestep_;
    }
}  // namespace SphaeroSim