#include "SimulationVTKInput.h"

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
    SimulationVTKInput::SimulationVTKInput(const double timestep,
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
            folder_timesteps_(folder_timesteps),
            folder_start_time_zero_(folder_start_time_zero),
            folder_temperature_indicator_(folder_temperature_indicator),
            folder_velocity_indicator_(folder_velocity_indicator),
            element_ids_(element_ids) {

        folder_list_ = folder_list;
        if (folder_list_.size() != folder_timesteps_.size()) {
            throw Exception("Wrong input parameters",
                            "The number of folder has to equal the number " \
                    "of folder timesteps.");
        }
        if (folder_list_.size() != folder_start_time_zero_.size()) {
            throw Exception("Wrong input parameters",
                            "The number of folder has to equal the number " \
                    "of folder start-time-zero flags.");
        }
        if (folder_list_.size() != folder_temperature_indicator_.size()) {
            throw Exception("Wrong input parameters",
                            "The number of folder has to equal the number " \
                    "of folder temperature indicators.");
        }
        if (folder_list_.size() != folder_velocity_indicator_.size()) {
            throw Exception("Wrong input parameters",
                            "The number of folder has to equal the number " \
                    "of folder velocity indicators.");
        }

        // check for the correct foldersyntax
        for (std::size_t cur = 0; cur < folder_list_.size(); ++cur) {
            std::replace(folder_list_[cur].begin(), folder_list_[cur].end(), '\\', '/');
            if (folder_list_[cur].back() != '/') {
                folder_list_[cur].push_back('/');
            }
        }
        SetNumRuns(1);
    }

    void SimulationVTKInput::FilenameArray(const std::string &folderpath,
                                           std::vector<std::string> *dst) const {
        Directory dir(folderpath.c_str());
        for (int cur_file = 0; cur_file < dir.elements.size(); ++cur_file) {
            const std::string filename = dir.getFileNameFromIndex(cur_file);
            if (filename.size() <= 3) {
                continue;
            }
            std::size_t last_pos = filename.size() - 1;
            if (filename[last_pos - 2] == 'v' &&
                filename[last_pos - 1] == 't' &&
                filename[last_pos] == 'k') {
                dst->push_back(filename);
            }
        }

        class Comperator {
        public:
            bool operator()(const std::string &left,
                            const std::string &right) {
                std::string::const_iterator lit, rit;
                for (lit = left.begin(), rit = right.begin(); lit != left.end() && rit != right.end(); ++lit, ++rit)
                    if (tolower(*lit) < tolower(*rit)) {
                        return true;
                    } else if (tolower(*lit) > tolower(*rit)) {
                        return false;
                    }
                if (left.size() < right.size()) {
                    return true;
                }
                return false;
            }
        } sort_obj;
        std::sort(dst->begin(), dst->end(), sort_obj);
    }

    void SimulationVTKInput::BuildElements(std::vector<Element *> *dst,
                                           RandomNumberGenerator *rng) {
        // create all the vtkfiles
        std::vector<std::pair<std::size_t, VTKFile *> > vtkfiles;
        for (std::size_t cur = 0; cur < folder_list_.size(); ++cur) {
            std::vector<std::string> filenames;
            FilenameArray(folder_list_[cur], &filenames);

            std::size_t start_file = 0;
            if (cur != 0 && folder_start_time_zero_[cur]) {
                start_file = 1;
            }
            for (std::size_t cur_file = start_file; cur_file < filenames.size(); ++cur_file) {
                std::stringstream ss;
                ss << folder_list_[cur] << filenames[cur_file];
                vtkfiles.push_back(std::pair<std::size_t, VTKFile *>(cur,
                                                                     new VTKFile(ss.str())));
            }
        }

        if (vtkfiles.empty()) {
            throw Exception("VTKFiles",
                            "No vtk files as input data were found.");
        }

        // read the files in
#pragma omp parallel for
        for (int64_t cur = 0; cur < static_cast<int64_t>(vtkfiles.size()); ++cur)
            vtkfiles[cur].second->ReadFromFile();

        // create the geometry information
        std::vector<std::vector<std::size_t> > indices;
        std::vector<Eigen::Vector3d> locations;
        std::vector<Eigen::Vector3d> sizes;
        vtkfiles[0].second->CellGeometry(&locations,
                                         &sizes,
                                         &indices);
        geometry_information_.locations = locations;
        geometry_information_.sizes = sizes;
        geometry_information_.indices = indices;

        const std::size_t num_elements = (element_ids_.empty()) ?
                                         indices.size() :
                                         element_ids_.size();
//#pragma omp parallel for
        for (int64_t cur_el = 0; cur_el < static_cast<int64_t>(num_elements); ++cur_el) {
            const std::size_t cur_element_id = (element_ids_.empty()) ?
                                               cur_el :
                                               element_ids_[cur_el];
            auto *new_element =
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
            double current_time = 0;
            // get the temperature and velocity fields
            for (auto &cur_file : vtkfiles) {
                VTKFile *vtkfile = cur_file.second;
                std::size_t index = cur_file.first;

                std::vector<double> temperatures;
                if (folder_temperature_indicator_[index] == "none") {
                    for (int iX = 0; iX < 8; ++iX)
                        temperatures.push_back(0.0);
                } else {
                    vtkfile->GetPointFieldValues(indices[cur_element_id],
                                                 folder_temperature_indicator_[index],
                                                 &temperatures);
                }

                std::vector<Eigen::Vector3d> velocities;
                if (folder_velocity_indicator_[index] == "none") {
                    for (int iX = 0; iX < 8; ++iX)
                        velocities.push_back(Eigen::Vector3d(0.0, 0.0, 0.0));
                } else {
                    vtkfile->GetPointFieldValues(indices[cur_element_id],
                                                 folder_velocity_indicator_[index],
                                                 &velocities);
                }

                bool should_store_values = true;
                for (double temperature : temperatures) {
                    if (temperature < 1e-6) {
                        should_store_values = false;
                        break;
                    }
                }
                if (should_store_values) {
                    // store the fields
                    temperature_field.push_back(Element::temperature_point(current_time,
                                                                           temperatures));
                    velocity_field.push_back(Element::velocity_point(current_time,
                                                                     velocities));
                }
                current_time += folder_timesteps_[index];
            }
            new_element->SetTemperatureField(temperature_field);
            new_element->SetVelocityField(velocity_field);
#pragma omp critical
            {
                dst->push_back(new_element);
            }
        }
    }

}  // namespace SphaeroSim
