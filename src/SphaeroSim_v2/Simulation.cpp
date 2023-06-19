#include "Simulation.h"

#ifndef _INCLUDED_STRING_H_
#include <string>
#define _INCLUDED_STRING_H_
#endif

#ifndef _INCLUDED_VECTOR_H_
#include <vector>
#define _INCLUDED_VECTOR_H_
#endif

#ifndef _INCLUDED_FSTREAM_H_

#include <fstream>

#define _INCLUDED_FSTREAM_H_
#endif

#ifndef _INCLUDED_SSTREAM_H_

#include <sstream>

#define _INCLUDED_SSTREAM_H_
#endif

#ifndef _INCLUDED_ALGORITHM_H_

#include <algorithm>

#define _INCLUDED_ALGORITHM_H_
#endif

#ifndef _INCLUDED_MUTEX_H_
#include <mutex>
#define _INCLUDED_MUTEX_H_
#endif

#ifndef _INCLUDED_CTIME_H_
#include <ctime>
#define _INCLUDED_CTIME_H_
#endif

#ifndef _SPHAEROSIM_STATE_EVENT_H_
#include "StateEvent.h"
#endif

#ifndef _SPHAEROSIM_ELEMENT_H_
#include "Element.h"
#endif

#ifndef __DIRECTORY_H

#include "directory.h"

#endif

#ifndef _SPHAEROSIM_PROBE_H_

#include "Probe.h"

#endif

#ifndef _SPHAEROSIM_RANDOM_NUMBER_GENERATOR_H_

#include "RandomNumberGenerator.h"

#endif

#ifndef _SPHAEROSIM_SPHERULITE_INFO_H_

#include "SpheruliteInfo.h"

#endif

#ifndef _SPHAEROSIM_GROWTH_MODEL_H_

#include "GrowthModel.h"

#endif

#ifndef _SPHAEROSIM_NUCLEATION_MODEL_H_

#include "NucleationModel.h"

#endif

#ifndef _SPHAEROSIM_VERSION_H_

#include "Version.h"

#endif

#undef max
#undef min

namespace SphaeroSim {

    Simulation::
    Simulation(const std::string &name,
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
            name_(name),
            convergence_criterium_(convergence_criterium),
            cell_size_m_(cell_size_m),
            temporal_interpolation_type_(temporal_interpolation_type),
            spatial_interpolation_type_(spatial_interpolation_type),
            create_subfolder_(create_new_subfolder),
            root_folder_(output_folder),
            vtk_result_option_(vtk_result_option),
            csv_separator_(csv_separator),
            calibration_file_(calibration_file),
            integral_error_(integral_error),
            statistics_step_(statistics_step),
            result_step_(result_step),
            save_memory_(save_memory),
            calculate_flow_energy_(calculate_flow_energy),
            cache_convert_temperature_(cache_convert_temperature),
            cache_convert_velocity_(cache_convert_velocity),
            numerical_options_(numerical_options) {

        material_ = material;

        current_run_ = 0;
        num_runs_ = 0;
        growth_model_ = growth_model;
        nucleation_model_ = nucleation_model;

        output_folder_ = "";
        noslip_neighbourhood_ = noslip_neighbourhood;

        // ensure at least entries for the the neighbourhood
        while (noslip_neighbourhood_.size() < 6) {
            noslip_neighbourhood_.push_back(true);
        }

        // random number generator
        rng_state_registers_ = rng_state_registers;
        if (!rng_state_registers_.empty()) {
            std::size_t cur_index = 0;
            while (rng_state_registers_.size() != R) {
                rng_state_registers_.push_back(rng_state_registers[cur_index]);
                ++cur_index;
                if (cur_index > rng_state_registers.size()) {
                    cur_index = 0;
                }
            }
        }
    }

    Simulation::~Simulation() {
        ClearElements();


        delete growth_model_;


        delete nucleation_model_;
        //TODO: Analyse

    }

    void Simulation::ClearElements() {
        for (auto &it:elements_) {
            delete it.second;
        }
        elements_.clear();

        for (auto &sph_it : spherulite_infos_) {
            delete sph_it.second;
        }
        spherulite_infos_.clear();

        for (auto &probe : probes_)
            delete probe;
        probes_.clear();
    }

// create the outputfolder
    const std::string Simulation::
    CreateOutputFolder(const bool create_subfolder,
                       const std::size_t current_run_index) const {
        std::stringstream ss;
        printf("%s\n", root_folder_.c_str());
        if (std::getenv("OUTPUT_DIR") != nullptr) {
            ss << std::getenv("OUTPUT_DIR");
        } else {
            ss << root_folder_;
        }
        if (root_folder_.back() != '\\' && root_folder_.back() != '/') {
            ss << "/";
        }

        std::string folder_path = ss.str();
        if (create_subfolder != 0) {
            ss << GetName() << "/";
            // current day
            if (std::getenv("SLURM_JOB_ID") == nullptr) {
                time_t now = time(0);
                struct tm time_struct{};
                char buf[80];
                time_struct = *localtime(&now);
                memset(buf, 0, 80);
                strftime(buf, sizeof(buf), "%Y%m%d", &time_struct);
                ss << buf;
            } else {
                ss << std::getenv("SLURM_JOB_ID");
            }
            ss << "_" << current_run_index + 1;

            if (Directory::folderExists(ss.str())) {
                std::string cur_folder;
                std::size_t folder_index = 1;
                do {
                    std::stringstream folder_stream("");
                    folder_stream << ss.str() << "_" << folder_index;
                    cur_folder = folder_stream.str();
                    folder_index++;
                }
                while (Directory::folderExists(cur_folder));
                folder_path = cur_folder;
            } else {
                folder_path = ss.str();
            }
        }
        std::replace(folder_path.begin(), folder_path.end(), '\\', '/');

        if (!Directory::folderExists(folder_path)) {
            if (!Directory::create_directory(folder_path)) {
                throw Exception("Internal Error", "Couldn't create output directory");
            }
        }

        if (folder_path.back() != '/') {
            folder_path.push_back('/');
        }
        return folder_path;
    }

// calculate the relative crystallinity of the entire simulation area
    const double Simulation::RelativeCrystallinity() const {
        // calculate the total relative crystallization degree
        double rel_cryst_deg = 0.0;
        for (const auto &it: elements_) {
            rel_cryst_deg += it.second->RelativeCrystallizationDegree();
        }
        rel_cryst_deg /= static_cast<double>(elements_.size());
        return rel_cryst_deg;
    }

// return the total number of spherulites
    const uint64_t Simulation::TotalNumberSpherulites() const {
        return spherulite_infos_.size();
    }

    void Simulation::WriteSpheruliteCenterFile() {
        std::stringstream ss;

        const uint64_t number_spherulites = TotalNumberSpherulites();

        // store coordinates
        ss << "#" << std::endl;
        ss << "# ID	X	Y	Z Koordinaten" << std::endl;
        ss << "#" << std::endl;
        ss << "/SPHERU CENTER, NUMBER=" << number_spherulites << std::endl;
        for (uint64_t cur_sph = 0; cur_sph < number_spherulites; ++cur_sph) {
            uint64_t sph_id = cur_sph + 2;
            auto info = SpheruliteCenter(sph_id);
            auto center = info.first->CellPosition(info.second, true);

            ss << sph_id << ",\t" << center[0] << ",\t";
            ss << center[1] << ",\t" << center[2] << std::endl;
        }
        ss << std::endl;

        // filename of the spherulite file
        std::stringstream ss_file;
        ss_file << GetOutputFolder() << "SpheruliteCenter.sph";
        std::string filename = ss_file.str();

        // save to file
        FILE *file;
        file = fopen(filename.c_str(), "w");
        fprintf(file, "%s", ss.str().c_str());
        fclose(file);
    }

// update the Spherulites.csv file
    void Simulation::UpdateSpheruliteFile() const {
        std::stringstream filename;
        if (GetVTKOptions()->GetSubfolder().empty()) {
            filename << GetOutputFolder();
        } else {
            filename << GetOutputFolder() << GetVTKOptions()->GetSubfolder() << ".";
        }
        filename << "Spherulites.csv";

        FILE *file;
        file = fopen(filename.str().c_str(), "w");

        if (file == nullptr) {
            printf("  -> Error: Could not create '%s'\n", filename.str().c_str());
            return;
        }

        // headerline
        fprintf(file,
                "Spherulite_ID%cCreation Time [s]%cDiameter [um]%cEq. Diameter [um]%cBox Size X [um]%cBox Size Y [um]%cBox Size Z [um]%c",
                GetCSVSeparator(),
                GetCSVSeparator(),
                GetCSVSeparator(),
                GetCSVSeparator(),
                GetCSVSeparator(),
                GetCSVSeparator(),
                GetCSVSeparator());
        fprintf(file,
                "Box Center X [mm]%cBox Center Y [mm]%cBox Center Z [mm]%cStart Location X [mm]%cStart Location Y [mm]%cStart Location Z [mm]\n",
                GetCSVSeparator(),
                GetCSVSeparator(),
                GetCSVSeparator(),
                GetCSVSeparator(),
                GetCSVSeparator());

        // go over all spherulites and execute the export
        for (const auto &it:spherulite_infos_) {
            auto info = it.second;

            Eigen::Vector3d box_size = info->GetBoxSize();
            if (info->GetNumberCellsInSpherulite() == 1) {
                box_size[0] = box_size[1] = box_size[2] =
                        info->GetOwnerElement()->GetCellSizeM();
            }

            fprintf(file, "%d%c%f%c%f%c%f%c%f%c%f%c%f%c%f%c%f%c%f%c%f%c%f%c%f\n",
                    static_cast<int32_t>(info->GetID()),
                    GetCSVSeparator(),
                    info->GetCreationTime(),
                    GetCSVSeparator(),
                    info->GetDiameter() * 1e6,
                    GetCSVSeparator(),
                    info->EquivalentDiameter() * 1e6,
                    GetCSVSeparator(),
                    box_size[0] * 1e6,
                    GetCSVSeparator(),
                    box_size[1] * 1e6,
                    GetCSVSeparator(),
                    box_size[2] * 1e6,
                    GetCSVSeparator(),
                    info->GetBoxCenter()[0] * 1e3,
                    GetCSVSeparator(),
                    info->GetBoxCenter()[1] * 1e3,
                    GetCSVSeparator(),
                    info->GetBoxCenter()[2] * 1e3,
                    GetCSVSeparator(),
                    info->GetStartLocationAbsolute()[0] * 1e3,
                    GetCSVSeparator(),
                    info->GetStartLocationAbsolute()[1] * 1e3,
                    GetCSVSeparator(),
                    info->GetStartLocationAbsolute()[2] * 1e3);
        }

        fclose(file);
    }

// update the Statistics.csv file
    void Simulation::UpdateStatisticsFile(const int32_t step_index,
                                          const double current_time) const {
        std::stringstream filename;
        if (GetVTKOptions()->GetSubfolder().empty()) {
            filename << GetOutputFolder();
        } else {
            filename << GetOutputFolder() << GetVTKOptions()->GetSubfolder() << ".";
        }
        filename << "Statistics.csv";

        FILE *file;
        if (Directory::file_exists(filename.str()) && step_index != 0) {
            file = fopen(filename.str().c_str(), "a");
            if (file == nullptr) {
                printf("  -> Error: Could not open '%s'\n", filename.str().c_str());
                return;
            }
        } else {
            file = fopen(filename.str().c_str(), "w");

            if (file == nullptr) {
                printf("  -> Error: Could not create '%s'\n", filename.str().c_str());
                return;
            }

            // print the headerline
            fprintf(file, "Step%c%s%cCalculation Time [s]%cSimulation Time [s]%cRel. Cryst. Deg [-]%cSpherulites [-]\n",
                    GetCSVSeparator(),
                    GetName().c_str(),
                    GetCSVSeparator(),
                    GetCSVSeparator(),
                    GetCSVSeparator(),
                    GetCSVSeparator());
        }

        // calculate the total relative crystallization degree
        const double rel_cryst_deg = RelativeCrystallinity();

        const double delta_time =
                static_cast<double>(clock() - last_clock_value_) / CLOCKS_PER_SEC;
        // append the data
        fprintf(file, "%d%c%lu%c%f%c%f%c%f%c%lu\n",
                step_index,
                GetCSVSeparator(),
                current_run_ + 1,
                GetCSVSeparator(),
                delta_time,
                GetCSVSeparator(),
                current_time,
                GetCSVSeparator(),
                rel_cryst_deg,
                GetCSVSeparator(),
                TotalNumberSpherulites());
        fclose(file);
    }

    void Simulation::CopyCalibrationFile() const {
        if (calibration_file_.empty()) {
            throw Exception("CopyCalibrationFile",
                            "No calibration file specified.");
        };

        std::ifstream stream_read;
        stream_read.open(calibration_file_.c_str());

        std::string filename = calibration_file_;
        std::replace(filename.begin(), filename.end(), '\\', '/');
        const std::size_t index = filename.find_last_of('/');
        if (index != std::string::npos && index != filename.size() - 1) {
            filename = filename.substr(index + 1);
        }

        std::stringstream ss_dst;
        if (GetVTKOptions()->GetSubfolder().empty()) {
            ss_dst << GetOutputFolder() << filename;
        } else {
            ss_dst << GetOutputFolder() << GetVTKOptions()->GetSubfolder() << "." << filename;
        }

        FILE *file_out = fopen(ss_dst.str().c_str(), "w");

        while (!stream_read.eof()) {
            std::string cur_line;
            std::getline(stream_read, cur_line);
            fprintf(file_out, _T("%s\n"), cur_line.c_str());
        }

        stream_read.close();
        fclose(file_out);
    }

    const double Simulation::UpdateElement(const double current_time,
                                           Element *element) {
        double last_executed_event = std::numeric_limits<double>::quiet_NaN();
        // execute the state events in the current interval
        for (auto it = element->GetStateEventIterator(); it != element->GetStateEventEnd();) {
            if ((*it)->GetPointInTime() > current_time) {
                ++it;
                break;
            }
            //TODO:Tasking
            (*it)->Execute();
            last_executed_event = (*it)->GetPointInTime();

            // remove the executed event
            it=element->RemoveStateEvent((it));
        }
        return last_executed_event;
    }

// create the probes for result report
    void Simulation::CreateProbes() {
        for (auto &probe : probes_)
            delete probe;
        probes_.clear();

        for (const auto &cur : probe_description_) {
            if (GetElement(cur.element_id_) == nullptr) {
                throw Exception("Probe:element_id", "Invalid element id.");
            }

            const Element *element = GetElement(cur.element_id_);

            // get the location and size of the probe-area
            Eigen::Vector3s start_cell, size_in_cells;

            start_cell =
                    IndexCoordinate(cur.location_type_,
                                    GetElement(cur.element_id_),
                                    cur.location_);
            if (cur.location_type_ == "Relative") {
                if (cur.type_ == "Element" ||
                    cur.type_ == "Cellblock") {
                    size_in_cells = element->GetNumCells();
                    for (std::size_t cur_dim = 0; cur_dim < 3; ++cur_dim)
                        size_in_cells[cur_dim] =
                                static_cast<std::size_t>(element->GetNumCells()[cur_dim] *
                                                         cur.size_[cur_dim]);
                }
                if (!element->CoordinateInElement(start_cell)) {
                    throw Exception("Probe:location", "Out of the element boundaries");
                }
            } else if (cur.location_type_ == "Absolute") {
                if (!element->CoordinateInElement(cur.location_)) {
                    throw Exception("Probe:location", "Out of the element boundaries");
                }
                if (cur.type_ == "Element" ||
                    cur.type_ == "Cellblock") {
                    Eigen::Vector3s end_cell =
                            element->CellIndicesFromPosition(cur.location_ +
                                                             cur.size_,
                                                             false);
                    if (end_cell[0] < start_cell[0] ||
                        end_cell[1] < start_cell[1] ||
                        end_cell[2] < start_cell[2]) {
                        throw Exception("Probe:location", "Invalid location size.");
                    }
                    size_in_cells = end_cell - start_cell;
                }
            } else if (cur.location_type_ == "Index") {
                if (cur.type_ == "Element" ||
                    cur.type_ == "Cellblock") {
                    size_in_cells = Eigen::Vector3s(
                            static_cast<std::size_t>(cur.size_[0]),
                            static_cast<std::size_t>(cur.size_[1]),
                            static_cast<std::size_t>(cur.size_[2]));
                }
                if (!element->CoordinateInElement(start_cell)) {
                    throw Exception("Probe:location", "Out of the element boundaries");
                }
            } else {
                throw Exception("Probe:location_type",
                                "Invalid probe-description. Allowed values are " \
                      "'Relative', 'Absolute' and 'Index'");
            }

            if (cur.type_ == "Element") {
                probes_.push_back(new ProbeElement(element->GetElementID(),
                                                   start_cell,
                                                   size_in_cells,
                                                   cur.name_));
            } else {
                if (cur.type_ == "Cell") {
                    probes_.push_back(new ProbeCell(element->GetElementID(),
                                                    element->AbsoluteCellIndex(start_cell),
                                                    cur.name_));
                } else {
                    for (std::size_t cur_x = 0; cur_x < size_in_cells[0]; ++cur_x) {
                        if (cur_x + start_cell[0] > element->GetNumCells()[0]) {
                            continue;
                        }
                        for (std::size_t cur_y = 0; cur_y < size_in_cells[1]; ++cur_y) {
                            if (cur_y + start_cell[1] > element->GetNumCells()[1]) {
                                continue;
                            }
                            for (std::size_t cur_z = 0; cur_z < size_in_cells[2]; ++cur_z) {
                                if (cur_z + start_cell[2] > element->GetNumCells()[2]) {
                                    continue;
                                }
                                const Eigen::Vector3s cur_cell(start_cell[0] + cur_x,
                                                               start_cell[1] + cur_y,
                                                               start_cell[2] + cur_z);
                                std::stringstream ss;
                                ss << cur.name_ << "(" << cur_cell[0];
                                ss << "," << cur_cell[1] << "," << cur_cell[2] << ")";
                                probes_.push_back(new ProbeCell(element->GetElementID(),
                                                                element->AbsoluteCellIndex(cur_cell),
                                                                ss.str()));
                            }
                        }
                    }
                }
            }
        }
    }

    void Simulation::ListOfElementIDs(std::vector<uint64_t> *dst) const {
        for (auto &it:elements_) {
            dst->push_back(it.first);
        }
    }

    void Simulation::SaveRNGRegisters(const std::string &filename,
                                      const RandomNumberGenerator *rng) const {
        std::vector<unsigned int> registers;
        rng->GetStateRegisters(&registers);

        std::stringstream filename_stream;
        filename_stream << GetOutputFolder() << GetVTKOptions()->GetSubfolder() << ".";
        filename_stream << filename;

        FILE *file = fopen(filename_stream.str().c_str(), "w");

        fprintf(file, "SphaeroSim v%d.%d.%d\n",
                MAJOR_VERSION, MINOR_VERSION, UPDATE_VERSION);
        for (std::size_t cur = 0; cur < registers.size(); ++cur) {
            fprintf(file, "%u", registers[cur]);
            if (cur != registers.size() - 1) {
                fprintf(file, "%c", GetCSVSeparator());
            }
        }

        fclose(file);
    }

    bool Simulation::ResultStepFullfilled(const ResultStep &result_step,
                                          const double current_time,
                                          double *last_evaluation_point) const {
        if (result_step.GetStepsize() <= 0.0) {
            return false;
        }

        double value_to_check;
        if (result_step.GetType() == ResultStep::kTime) {
            value_to_check = current_time;
        } else {  // crystallization degree
            value_to_check = RelativeCrystallinity();
        }
        if (value_to_check + 1e-9 >=
            *last_evaluation_point + result_step.GetStepsize()) {
            *last_evaluation_point = value_to_check;
            return true;
        }
        return false;
    }

    void Simulation::Run() {
        printf("SphaeroSim v%d.%d.%d\n",
               MAJOR_VERSION, MINOR_VERSION, UPDATE_VERSION);
        printf("Running simulation '%s'\n", GetName().c_str());

        // execute the runs
        for (current_run_ = 0; current_run_ < GetNumRuns(); ++current_run_) {
            output_folder_ = CreateOutputFolder(create_subfolder_,
                                                GetCurrentRunIndex());
            CopyCalibrationFile();
            printf("  Executing run %zd of %zd\n", GetCurrentRunIndex() + 1, GetNumRuns());

            // create the random number generator
            auto *rng = new RandomNumberGenerator();
            if (!rng_state_registers_.empty()) {
                rng->SetStateRegisters(rng_state_registers_);
            }

            // save the registers
            SaveRNGRegisters("RNGRegister.txt",
                             rng);

            // create the elements
            std::vector<Element *> elements;
            printf("  -> Loading element boundary conditions...");
            fflush(stdout);
            BuildElements(&elements, rng);
            for (auto &element : elements)
                elements_[element->GetElementID()] = element;
            printf("completed\n");

            // set the neighbourhood-information
            printf("  -> Loading element neighbourhood...");
            fflush(stdout);
            SetElementNeighbourhood(&elements_);
            printf("completed\n");

            // precalculate the boundary fields
            printf("  -> Precalculating spatial temperature and velocities...");
            fflush(stdout);
            for (auto &element : elements)
                elements_[element->GetElementID()]->PrecalculateBoundaryFields();
            printf("completed\n");

            // create the initial events
            printf("  -> Creating initial events...");
            fflush(stdout);
            CreateInitialEvents();
            printf("completed\n");

            // create the probes
            printf("  -> Creating probes...");
            fflush(stdout);
            CreateProbes();
            printf("completed\n");

            // spread the initial nuclei
            printf("  -> Creating predefined nuclei distributions...");
            fflush(stdout);
            for (auto &cur : predefined_nuclei_)
                cur->SpreadInElement(*this);
            printf("completed\n");

            double current_time = 0.0;
            double last_printout = 0.0;
            double last_statistic_printout = 0.0;
            double last_result_printout = 0.0;

            last_clock_value_ = clock();
            double last_convergence_value = 0.0;

            int32_t statistic_export_counter = 0;
            double start_glob = omp_get_wtime();
            while (!ResultStepFullfilled(convergence_criterium_,
                                         current_time,
                                         &last_convergence_value)) {
                // update the elements
                //TODO:parallelise
                for (auto &cur_element:elements_) {
                    if (cur_element.second->GetFullyCrystalline()) {
                        continue;
                    }
                    UpdateElement(current_time, (cur_element).second);
                }
                // result step file
                double new_last_result_printout = last_result_printout;
                if (ResultStepFullfilled(result_step_,
                                         current_time,
                                         &new_last_result_printout)) {
                    for (auto &cur_element:elements_) {
                        WriteResultEvent result_event(current_time,
                                                      -1.0,  // no automatic generation of new events
                                                      *this,
                                                      (cur_element).second);
                        result_event.Execute();
                    }
                    WriteSpheruliteCenterFile();
                    last_result_printout = new_last_result_printout;
                }

                current_time += Timestep();

                if (current_time > last_printout + Timestep()) {
                    printf("    %d ms of the calculation done (%.2f%% solid, %lu spherulites) (WTime: %fs)\n",
                           static_cast<int32_t>(current_time * 1000.0),
                           RelativeCrystallinity() * 100.0,
                           TotalNumberSpherulites(),
                           omp_get_wtime() - start_glob);
                    last_printout = current_time;
                }

                // statistics file
                if (ResultStepFullfilled(statistics_step_,
                                         current_time,
                                         &last_statistic_printout)) {
                    UpdateStatisticsFile(statistic_export_counter, current_time);
                    UpdateSpheruliteFile();
                    statistic_export_counter += 1;
                }
            }
            // write the final state of the elements
            for (auto &cur_element:elements_) {
                WriteResultEvent result_event(current_time,
                                              -1.0,  // no automatic generation of new events
                                              *this,
                                              (cur_element).second);
                result_event.Execute();
            }
            WriteSpheruliteCenterFile();

            // last update of the statistics
            UpdateStatisticsFile(statistic_export_counter, current_time);
            UpdateSpheruliteFile();

            delete rng;
            ClearElements();
        }
    }

    void Simulation::CreateInitialEvents() {
        for (auto &cur_element:elements_) {
            if (cur_element.second->GetCrystalGrowthType() == GrowthEvent::kCCG) {
                const std::size_t undefined = std::numeric_limits<std::size_t>::max();
                // initial growth event
                cur_element.second->AddStateEvent(
                        new GrowthEvent(undefined,
                                        complex_index(cur_element.second, undefined),
                                        Timestep(),
                                        Timestep(),
                                        GrowthEvent::kCCG,
                                        Eigen::Vector3i(0, 0, 0),
                                        cur_element.second));


            }

            // initial nucleation event
            cur_element.second->AddStateEvent(
                    new NucleationEvent(Timestep(),
                                        Timestep(),
                                        0.0,
                                        cur_element.second));
        }

    }

// add a new spherulite to the pool
    const uint64_t Simulation::NewSpherulite(const std::size_t start_location,
                                             const double creation_time,
                                             const Element *owner) {
        // the first spherulite ID is two. This allows to store the,
        // spherulite id of each cell in the phasestate_ array of the
        // element. Phasestates in the range from -1.0 to 1.0 are
        // reserved for the calculation of the nucleation process
        uint64_t spherulite_id;

        spherulite_mutex_.lock();

        union RngColor {
            int32_t random_color32;
            int16_t random_colors16[2];
        } rng_color{};
        rng_color.random_colors16[0] = static_cast<int16_t>(rand());
        rng_color.random_colors16[1] = static_cast<int16_t>(rand());

        spherulite_id = spherulite_infos_.size() + 2;
        auto new_spherulite = new SpheruliteInfo(spherulite_id,
                                                 start_location,
                                                 creation_time,
                                                 rng_color.random_color32,
                                                 owner);
        spherulite_infos_.insert(std::pair<const uint64_t, SpheruliteInfo *>(
                spherulite_id,
                new_spherulite));
        spherulite_mutex_.unlock();

        return spherulite_id;
    }

// retreive the spherulite center
    const complex_index
    Simulation::SpheruliteCenter(const uint64_t spherulite_id) {
        SpheruliteInfo *cur_info;

        spherulite_mutex_.lock();
        cur_info = spherulite_infos_.at(spherulite_id);
        spherulite_mutex_.unlock();

        return complex_index(cur_info->GetOwnerElement(),
                             cur_info->GetStartLocation());
    }

// add a cell to the spherulite
    void Simulation::AddCellToSpherulite(const uint64_t spherulite_id,
                                         const complex_index &cell_to_add,
                                         const Eigen::Vector3i &periodic_continuation,
                                         const bool spherulitesAdded) {
        if (spherulitesAdded) {
            spherulite_mutex_.lock();
            spherulite_infos_.at(spherulite_id)->AddCellToSpherulite(cell_to_add,
                                                                     periodic_continuation);
            spherulite_mutex_.unlock();
        } else {
            spherulite_infos_.at(spherulite_id)->AddCellToSpherulite(cell_to_add,
                                                                     periodic_continuation);
        }
    }

// retrieve the random color for a spherulite
    const int32_t Simulation::SpheruliteRandomColor(const uint64_t &spherulite_id) {
        spherulite_mutex_.lock();
        const int32_t rng_color = spherulite_infos_.at(spherulite_id)->GetRandomColor();
        spherulite_mutex_.unlock();
        return rng_color;
    }

    const int32_t Simulation::SpheruliteRandomColor(const SpheruliteInfo *const spherulite) {
        return spherulite->GetRandomColor();
    }

// get the diameter of a spherulite
    const double Simulation::SpheruliteDiameter(const uint64_t &spherulite_id) {
        spherulite_mutex_.lock();
        const double diameter = spherulite_infos_.at(spherulite_id)->GetDiameter();
        spherulite_mutex_.unlock();
        return diameter;
    }

    const double Simulation::SpheruliteDiameter(const SpheruliteInfo *const spherulite) {
        spherulite_mutex_.lock();
        const double diameter = spherulite->GetDiameter();
        spherulite_mutex_.unlock();
        return diameter;
    }

    //if no cells are added to Spherulite while this
    const double Simulation::SpheruliteDiameterNoDataracePossible(const SpheruliteInfo *const spherulite) {
        return spherulite->GetDiameter();
    }


// retrieve the aspect ratio of a spherulite
    const double Simulation::SpheruliteAspectRatio(const uint64_t &spherulite_id) {
        spherulite_mutex_.lock();
        const double aspect_ratio =
                spherulite_infos_.at(spherulite_id)->AspectRatio(GetCellSizeM());
        spherulite_mutex_.unlock();
        return aspect_ratio;
    }

    const double Simulation::SpheruliteAspectRatio(const SpheruliteInfo *const spherulite) {
        spherulite_mutex_.lock();
        const double aspect_ratio = spherulite->AspectRatio(GetCellSizeM());
        spherulite_mutex_.unlock();
        return aspect_ratio;
    }

    const double Simulation::SpheruliteAspectRatioNoDatraracePossible(const SpheruliteInfo *const spherulite) {
        return spherulite->AspectRatio(GetCellSizeM());
    }

// get the average diameter and standard deviation of the spherulite diameter
// for all spherulites which have their origin in a specific element
// return: <average diameter, standard deviation>
    const std::pair<const double, const double> Simulation::
    AverageSpheruliteDiameterForElement(const Element *element,
                                        const Eigen::Vector3s start_location,
                                        const Eigen::Vector3s num_cells,
                                        std::size_t *dst_num_spherulites) {
        spherulite_mutex_.lock();


        // average value
        double summed_diameter = 0.0;
        std::size_t num_spherulites = 0;
        for (auto &it:spherulite_infos_) {
            const SpheruliteInfo *info = it.second;

            const Eigen::Vector3s spherulite_start =
                    info->GetOwnerElement()->CellIndices(info->GetStartLocation());
            if (spherulite_start[0] >= start_location[0] &&
                spherulite_start[1] >= start_location[1] &&
                spherulite_start[2] >= start_location[2] &&
                spherulite_start[0] < start_location[0] + num_cells[0] &&
                spherulite_start[1] < start_location[1] + num_cells[1] &&
                spherulite_start[2] < start_location[2] + num_cells[2]) {
                if (info->GetOwnerElement() == element) {
                    summed_diameter += info->GetDiameter();
                    ++num_spherulites;
                }
            }
        }
        if (num_spherulites == 0) {
            spherulite_mutex_.unlock();
            return std::pair<const double, const double>(0.0, 0.0);
        } else if (num_spherulites == 1) {
            spherulite_mutex_.unlock();
            return std::pair<const double, const double>(summed_diameter, 0.0);
        }
        const double average_diameter =
                summed_diameter / static_cast<double>(num_spherulites);

        // standard deviation
        double error_sum = 0.0;
        for (auto &it: spherulite_infos_) {
            const SpheruliteInfo *info = it.second;
            const Eigen::Vector3s spherulite_start =
                    info->GetOwnerElement()->CellIndices(info->GetStartLocation());
            if (spherulite_start[0] >= start_location[0] &&
                spherulite_start[1] >= start_location[1] &&
                spherulite_start[2] >= start_location[2] &&
                spherulite_start[0] < start_location[0] + num_cells[0] &&
                spherulite_start[1] < start_location[1] + num_cells[1] &&
                spherulite_start[2] < start_location[2] + num_cells[2]) {
                if (info->GetOwnerElement() == element) {
                    const double delta = average_diameter - info->GetDiameter();
                    error_sum += delta * delta;
                }
            }
        }
        const double standard_deviation =
                sqrt(error_sum / static_cast<double>(num_spherulites - 1));
        spherulite_mutex_.unlock();

        if (dst_num_spherulites != nullptr) {
            *dst_num_spherulites = num_spherulites;
        }

        return std::pair<const double, const double>(average_diameter,
                                                     standard_deviation);
    }
}  // namespace SphaeroSim