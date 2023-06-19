#include "SimulationComsolInput.h"

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

#ifndef _COMSOL_FILE_H_

#include "ComsolFile.h"

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
    SimulationComsolInput::
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
                          Material *material) :
            SimulationInputFiles(timestep,
                                 name,
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
                                 material),
            coordinate_x_identifier_(coordinate_x_identifier),
            coordinate_y_identifier_(coordinate_y_identifier),
            coordinate_z_identifier_(coordinate_z_identifier),
            velocity_x_identifier_(velocity_x_identifier),
            velocity_y_identifier_(velocity_y_identifier),
            velocity_z_identifier_(velocity_z_identifier),
            temperature_identifier_(temperature_identifier),
            levelset_identifier_(levelset_identifier) {

        file_list_ = file_list;
        if (file_list_.size() == 0) {
            throw Exception("Wrong input parameters",
                            "The number of files has to be larger than zero");
        }

        // check for the correct foldersyntax
        for (std::size_t cur = 0; cur < file_list_.size(); ++cur) {
            std::replace(file_list_[cur].begin(), file_list_[cur].end(), '\\', '/');
        }
        SetNumRuns(1);
    }

    void SimulationComsolInput::BuildElements(std::vector<Element *> *dst,
                                              RandomNumberGenerator *rng) {
        // create all the vtkfiles
        std::vector<ComsolFile *> comsolfiles;
        for (std::size_t cur = 0; cur < file_list_.size(); ++cur) {
            ComsolFile *new_file = new ComsolFile(file_list_[cur],
                                                  coordinate_x_identifier_,
                                                  coordinate_y_identifier_,
                                                  coordinate_z_identifier_,
                                                  velocity_x_identifier_,
                                                  velocity_y_identifier_,
                                                  velocity_z_identifier_,
                                                  temperature_identifier_,
                                                  levelset_identifier_,
                                                  GetCellSizeM());
            comsolfiles.push_back(new_file);
        }

        if (comsolfiles.size() == 0) {
            throw Exception("ComsolFile",
                            "No comsol files as input data were found.");
        }

        // read the files in
        //TODO: analyse
//#pragma omp parallel for
        for (auto &comsolfile : comsolfiles)
            comsolfile->ReadFromFile();

        // ensure all files have the same number of elements
        const int64_t num_elements = comsolfiles[0]->GetNumberElements();
        for (auto &comsolfile : comsolfiles) {
            if (comsolfile->GetNumberElements() != num_elements) {
                throw Exception("ComsolFile",
                                "Not the same geometry information between the different comsol files");
            }
        }

        // create the geometry information
        std::vector<std::vector<std::size_t> > indices;
        std::vector<Eigen::Vector3d> locations;
        std::vector<Eigen::Vector3d> sizes;
        for (int64_t cur = 0; cur < num_elements; ++cur) {
            locations.push_back(comsolfiles[0]->GetElementPosition(cur));
            sizes.push_back(comsolfiles[0]->GetElementSize(cur));

            std::vector<std::size_t> cur_indices;
            comsolfiles[0]->ElementIndices(cur, &cur_indices);
            indices.push_back(cur_indices);
        }
        geometry_information_.locations = locations;
        geometry_information_.sizes = sizes;
        geometry_information_.indices = indices;

        //TODO: analyse
//#pragma omp parallel for
        for (int64_t cur_el = 0; cur_el < static_cast<int64_t>(num_elements); ++cur_el) {
            const std::size_t cur_element_id = cur_el;
            Element *new_element =
                    new Element(cur_element_id,
                                geometry_information_.locations[cur_element_id],
                                GetCellSizeM(),
                                geometry_information_.sizes[cur_element_id],
                                GetTemporalInterpolationType(),
                                GetSpatialInterpolationType(),
                                GetGrowthModel()->GetGrowthType(),
                                rng,
                                this);

            std::vector<Element::temperature_point> temperature_field;
            std::vector<Element::velocity_point> velocity_field;

            // get the temperature and velocity fields
            for (auto &curfile : comsolfiles) {
                std::vector<bool> valid_fields;
                std::vector<double> pt_in_time;

                if (curfile->GetHasField(ComsolFile::kLevelSetVariable)) {
                    std::vector<std::vector<double> > levelset_values;
                    // check the levelset value > 0.5 for all values
                    curfile->ElementFields(cur_element_id,
                                           ComsolFile::kLevelSetVariable,
                                           &pt_in_time,
                                           &levelset_values);
                    valid_fields.resize(pt_in_time.size(), true);

                    for (std::size_t cur_time = 0; cur_time < pt_in_time.size(); ++cur_time) {
                        for (int32_t cur_val = 0; cur_val < 8; ++cur_val) {
                            if (levelset_values[cur_time][cur_val] < 0.5) {
                                valid_fields[cur_time] = false;
                                break;
                            }
                        }
                    }
                }

                if (curfile->GetHasField(ComsolFile::kTemperature)) {
                    // temperatures
                    std::vector<std::vector<double> > temperature_values;
                    curfile->ElementFields(cur_element_id,
                                           ComsolFile::kTemperature,
                                           &pt_in_time,
                                           &temperature_values);
                    // velocities
                    std::vector<std::vector<double> > vel[3];
                    if (curfile->GetHasField(ComsolFile::kVelocityX)) {
                        curfile->ElementFields(cur_element_id,
                                               ComsolFile::kVelocityX,
                                               &pt_in_time,
                                               &vel[0]);
                    }
                    if (curfile->GetHasField(ComsolFile::kVelocityY)) {
                        curfile->ElementFields(cur_element_id,
                                               ComsolFile::kVelocityY,
                                               &pt_in_time,
                                               &vel[1]);
                    }
                    if (curfile->GetHasField(ComsolFile::kVelocityZ)) {
                        curfile->ElementFields(cur_element_id,
                                               ComsolFile::kVelocityZ,
                                               &pt_in_time,
                                               &vel[2]);
                    }

                    for (std::size_t cur_time = 0; cur_time < pt_in_time.size(); ++cur_time) {
                        if (!valid_fields.empty()) {
                            if (valid_fields[cur_time] == false) {
                                continue;
                            }
                        }
                        temperature_field.push_back(
                                Element::temperature_point(pt_in_time[cur_time],
                                                           temperature_values[cur_time]));

                        std::vector<Eigen::Vector3d> velocities;
                        for (int32_t cur_val = 0; cur_val < 8; ++cur_val) {
                            Eigen::Vector3d cur_vel(0.0, 0.0, 0.0);
                            if (vel[0].size() > 0) {
                                cur_vel[0] = vel[0][cur_time][cur_val];
                            }
                            if (vel[1].size() > 0) {
                                cur_vel[1] = vel[1][cur_time][cur_val];
                            }
                            if (vel[2].size() > 0) {
                                cur_vel[2] = vel[2][cur_time][cur_val];
                            }
                            velocities.push_back(cur_vel);
                        }
                        velocity_field.push_back(
                                Element::velocity_point(pt_in_time[cur_time],
                                                        velocities));
                    }
                } else {
                    throw Exception("Simulation - No temperatures specified", curfile->GetFilename());
                }
            }
            new_element->SetTemperatureField(temperature_field);
            new_element->SetVelocityField(velocity_field);
#pragma omp critical
            {
                dst->push_back(new_element);
            }
        }

        for (auto &comsolfile : comsolfiles) {
            delete comsolfile;
        }
    }

}  // namespace SphaeroSim
