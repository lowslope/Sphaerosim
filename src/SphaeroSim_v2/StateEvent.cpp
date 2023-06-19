#include "StateEvent.h"

#ifndef _INCLUDED_STDINT_H_
#include <stdint.h>
#define _INCLUDED_STDINT_H_
#endif

#ifndef _INCLUDED_SSTREAM_H_

#include <sstream>

#define _INCLUDED_SSTREAM_H_
#endif

#ifndef _INCLUDED_ALGORITHM_H_

#include <algorithm>

#define _INCLUDED_ALGORITHM_H_
#endif

#ifndef _INCLUDED_MAP_H_
#include <map>
#define _INCLUDED_MAP_H_
#endif

#ifndef _INCLUDED_IOMANIP_H_

#include <iomanip>

#define _INCLUDED_IOMANIP_H_
#endif

#ifndef _SPHAEROSIM_VTKFILE_H_

#include "VTKFile.h"

#endif

#ifndef _SPHAEROSIM_ELEMENT_H_

#include "Element.h"

#endif

#ifndef _SPHAEROSIM_EIGEN_LIBRARY_H_
#include "EigenLibrary.h"
#endif

#ifndef _SPHAEROSIM_SIMULATION_H_

#include "Simulation.h"

#endif

#ifndef __DIRECTORY_H

#include "directory.h"

#endif

#ifndef _SPHAEROSIM_PROBE_H_

#include "Probe.h"

#endif

#ifndef _UTILITY_HASH_VALUE_H_

#include "HashValue.h"

#endif

#ifndef _SPHAEROSIM_VERSION_H_

#include "Version.h"

#endif

namespace SphaeroSim {

    StateEvent::StateEvent(const Type type,
                           const double point_in_time,
                           const double offset_new_event,
                           Element *element) :
            type_(type),
            point_in_time_(point_in_time),
            offset_new_event_(offset_new_event) {
        element_ = element;
    }

    NullEvent::NullEvent(const double point_in_time,
                         Element *element) :
            StateEvent(kNullEvent,
                       point_in_time,
                       -1.0,
                       element) {
    }

    void NullEvent::Execute() {
    }

    bool NullEvent::SameEvent(const StateEvent *compare_to) const {
        if (compare_to->GetType() != kNullEvent) {
            return false;
        }
        if (compare_to->GetElement() != GetElement()) {
            return false;
        }
        if (fabs(compare_to->GetPointInTime() - GetPointInTime()) > 1e-10) {
            return false;
        }
        return true;
    }

    NucleationEvent::NucleationEvent(const double point_in_time,
                                     const double offset_new_event,
                                     const double last_nucleation,
                                     Element *element) :
            StateEvent(kNucleation,
                       point_in_time,
                       offset_new_event,
                       element),
            last_nucleation_(last_nucleation) {
    }

    void NucleationEvent::Execute() {
        const size_t number_cells = GetElement()->GetTotalNumCells();
        const double cell_volumne = GetElement()->GetCellSizeM() *
                                    GetElement()->GetCellSizeM() *
                                    GetElement()->GetCellSizeM();

        std::map<uint64_t, std::pair<std::size_t, double> > nucleation_rate_cache;
        std::vector<uint64_t> nucleation_keys;
        std::vector<double> temperature_times;
        std::vector<double> velocity_times;

        // get temperature_times and velocity_times
        GetElement()->CellHashkey(0,
                                  last_nucleation_,
                                  &temperature_times,
                                  &velocity_times);

        if (temperature_times.empty() && velocity_times.empty()) {
            // add the next nucleation event
            if (GetOffsetNewEvent() > 0.0) {
                GetElement()->AddStateEvent(
                        new NucleationEvent(GetPointInTime() + GetOffsetNewEvent(),
                                            GetOffsetNewEvent(),
                                            GetPointInTime(),
                                            GetElement()));
            }
            return;
        }

        // calculate all nucleation hashkeys
        nucleation_keys.resize(number_cells);

        std::vector<std::map<uint64_t, std::pair<std::size_t, double>>> thread_caches;
        thread_caches.resize(static_cast<unsigned long>(omp_get_max_threads()),
                             std::map<uint64_t, std::pair<std::size_t, double>>());

        #pragma omp parallel for
        for (size_t cur_cell = 0; cur_cell < number_cells; ++cur_cell) {
            nucleation_keys[cur_cell] = 0;
            if (!GetElement()->CellCrystalline(cur_cell)) {
                // calculate the cell-key
                nucleation_keys[cur_cell] =
                        GetElement()->CellHashkey(static_cast<std::size_t>(cur_cell),
                                                  GetPointInTime(),
                                                  &temperature_times,
                                                  &velocity_times);
                thread_caches[omp_get_thread_num()][nucleation_keys[cur_cell]] =
                        std::pair<std::size_t, double>(cur_cell, -1.0);
//#pragma omp critical
//      {
//      // create an entry in the cache
//      nucleation_rate_cache[nucleation_keys[cur_cell]] =
//        std::pair<std::size_t, double>(cur_cell, -1.0);
//      }
            }
        }
        //for (int64_t cur_cell = 0; cur_cell < number_cells; ++cur_cell) {
        //  if (GetElement()->CellCrystalline(cur_cell) == false)
        //    // create an entry in the cache
        //    nucleation_rate_cache[nucleation_keys[cur_cell]] =
        //      std::pair<std::size_t, double>(cur_cell, -1.0);
        //}
        nucleation_rate_cache = thread_caches[0];
        for (int cur = 1; cur < omp_get_max_threads(); ++cur)
            nucleation_rate_cache.insert(thread_caches[cur].begin(),
                                         thread_caches[cur].end());

        std::vector<uint64_t> key_list;
        key_list.resize(nucleation_rate_cache.size());
        // retrieve all keys from the map
        std::size_t counter = 0;
        for (auto it:nucleation_rate_cache) {
            key_list[counter++] = it.first;
        }

        // go over all cells and integrate the nucleation rates
        const auto number_keys = static_cast<int64_t>(key_list.size());
        #pragma omp parallel for
        for (int64_t cur_key = 0; cur_key < number_keys; ++cur_key) {
            const std::size_t cur_cell =
                    nucleation_rate_cache[key_list[cur_key]].first;
            const double integrate_cell_over_time =
                    GetElement()->IntegrateCellOverTime(last_nucleation_,
                                                        GetPointInTime(),
                                                        cur_cell,
                                                        true);

            nucleation_rate_cache[key_list[cur_key]].second = integrate_cell_over_time;
        }

        // assign the phasestates to all cells
        #pragma omp parallel for
        for (size_t cur_cell = 0; cur_cell < number_cells; ++cur_cell) {
            if (GetElement()->CellCrystalline(cur_cell)) {
                continue;
            }

            double delta_cell =
                    cell_volumne * nucleation_rate_cache[nucleation_keys[cur_cell]].second;

            if (delta_cell < 0.0) {
                throw Exception("Invalid nucleation model",
                                "The integration of the nucleation rate " \
                      "resulted in a negative value");
            }

            const auto cell_index = static_cast<std::size_t>(cur_cell);
            if (GetElement()->GetPhasestate(cell_index) + delta_cell > 1.0) {
                delta_cell = 1.0 - GetElement()->GetPhasestate(cell_index);
            }

            GetElement()->AdjustCellPhasestate(cell_index, delta_cell);
        }

        // execute a monte-carlo step for each cell
        bool has_new_nucleus = false;
        for (size_t cur_cell = 0; cur_cell < number_cells; ++cur_cell) {
            if (GetElement()->CellCrystalline(cur_cell)) {
                continue;
            }
            const double phasestate = GetElement()->GetPhasestate(cur_cell);
            if (GetElement()->MonteCarloStep(phasestate)) {
                GetElement()->AddNucleus(GetElement()->CellIndices(cur_cell),
                                         GetPointInTime());
                has_new_nucleus = true;
            }
        }
        //if (has_new_nucleus == true) {
        #pragma omp parallel for
        for (size_t cur_cell = 0; cur_cell < number_cells; ++cur_cell) {
            if (!GetElement()->CellCrystalline(cur_cell)) {
                GetElement()->phasestate_[cur_cell] = 0.0;
            }
        }
        //}

        // add the next nucleation event
        if (GetOffsetNewEvent() > 0.0) {
            GetElement()->AddStateEvent(
                    new NucleationEvent(GetPointInTime() + GetOffsetNewEvent(),
                                        GetOffsetNewEvent(),
                                        GetPointInTime(),
                                        GetElement()));
        }
    }

    bool NucleationEvent::SameEvent(const StateEvent *compare_to) const {
        if (compare_to->GetType() != kNucleation) {
            return false;
        }
        if (compare_to->GetElement() != GetElement()) {
            return false;
        }
        if (fabs(compare_to->GetPointInTime() - GetPointInTime()) > 1e-10) {
            return false;
        }
        return true;
    }

    GrowthEvent::GrowthEvent(const std::size_t cell_index,
                             const complex_index creating_cell,
                             const double point_in_time,
                             const double offset_new_event,
                             const Type growth_type,
                             const Eigen::Vector3i periodic_continuation,
                             Element *element) :
            StateEvent(kGrowth,
                       point_in_time,
                       offset_new_event,
                       element),
            cell_index_(cell_index),
            growth_type_(growth_type),
            creating_cell_(creating_cell),
            periodic_continuation_(periodic_continuation) {
    }

    void GrowthEvent::Execute() {
        if (growth_type_ != kCCG) {
            if (GetElement()->CellCrystalline(cell_index_)) {
                return;
            }  // already crystalline
        }

        if (growth_type_ == kMonteCarlo) {
            // calculate the proceeding of the growthfront from the neighbours
            double percentage_covered = 0.0;
            const double cryst_time =
                    creating_cell_.first->crystallization_time_[creating_cell_.second];
            const double distance =
                    GetElement()->IntegrateCellOverTime(cryst_time,
                                                        GetPointInTime(),
                                                        cell_index_,
                                                        false);
            Eigen::Vector3d dist = GetElement()->CellPosition(cell_index_, true) -
                                   creating_cell_.first->CellPosition(creating_cell_.second,
                                                                      true);
            const double half_dist = 0.5 * dist.norm();
            percentage_covered += distance / half_dist;

            if (GetElement()->MonteCarloStep(percentage_covered)) {
                // set crystalline
                const uint64_t spherulite_id =
                        creating_cell_.first->GetSpheruliteID(creating_cell_.second);
                GetElement()->ChangeCellToSolid(cell_index_,
                                                spherulite_id,
                                                GetPointInTime(),
                                                periodic_continuation_,
                                                false);
            }
        } else if (growth_type_ == kRayTracing) {
            // The raytracing method calculates the point in time when the cell is
            // fully covered by the growth front. This point in time equals the event
            // time => simply change the phase of the cell
            const uint64_t spherulite_id =
                    creating_cell_.first->GetSpheruliteID(creating_cell_.second);
            GetElement()->ChangeCellToSolid(cell_index_,
                                            spherulite_id,
                                            GetPointInTime(),
                                            periodic_continuation_,
                                            false);
        } else if (growth_type_ == kCCG) {
            // The ccg method works more like a nucleation event. Using a fixed
            // timestep the algorithm goes over all cells and adjust the values of
            // the neighbourhodd of crystalline cells
            std::map<complex_index, double> crystal_state_increase;
            std::map<complex_index, double> partial_increase;

            const size_t number_cells = GetElement()->GetTotalNumCells();
            #pragma omp parallel for
            for (size_t cur_cell = 0; cur_cell < number_cells; ++cur_cell) {
                const double cur_state = GetElement()->GetAdditionalCCGState(cur_cell);
                if (cur_state > 0.0) {
                    if (cur_state < 1.1) {
                        const double growth_speed =
                                GetElement()->CrystalGrowthSpeed(GetPointInTime(),
                                                                 cur_cell);
                        // GrowValue is 1/cellSize*timestep*Growthspeed
                        const double grow_value =
                                growth_speed * GetOffsetNewEvent() / GetElement()->GetCellSizeM();

                        IncreaseNeighbourCCG(complex_index(GetElement(), cur_cell),
                                             grow_value,
                                             &crystal_state_increase,
                                             &partial_increase);

                        #pragma omp critical
                        {
                            crystal_state_increase[complex_index(GetElement(), cur_cell)] =
                                    grow_value;
                        }
                    }
                }
            }
            // Calculate the test of the increased cells (Overflow)
            while (!partial_increase.empty()) {
                auto it = partial_increase.begin();
                const Element *cur_el = it->first.first;
                const std::size_t cell_index = it->first.second;
                const double increment = it->second;
                const double growth_speed =
                        cur_el->CrystalGrowthSpeed(GetPointInTime(),
                                                   cell_index);
                const double grow_value =
                        growth_speed * GetOffsetNewEvent() * cur_el->GetCellSizeM() * increment;
                IncreaseNeighbourCCG(it->first,
                                     grow_value,
                                     &crystal_state_increase,
                                     &partial_increase);
                partial_increase.erase(complex_index(cur_el, cell_index));
            }

            // add all increases to the cell values
            for (const auto &it:crystal_state_increase) {
                const Element *cur_el = it.first.first;
                const std::size_t cell_index = it.first.second;
                double increment = it.second;
                const double cur_state = cur_el->GetAdditionalCCGState(cell_index);
                if (cur_state < 1.1) {
                    if (cur_state + increment > 1.1) {
                        increment = 1.1 - cur_state;
                    }
                    const_cast<Element *>(cur_el)->IncreaseAdditionalCCGState(cell_index,
                                                                              increment);
                }
            }
        }

        bool add_next_event = true;
        if (growth_type_ != kCCG) {
            add_next_event = !GetElement()->CellCrystalline(cell_index_);
        }
        if (add_next_event && GetOffsetNewEvent() > 0.0 &&
            !GetElement()->GetFullyCrystalline()) {
            GetElement()->AddStateEvent(new GrowthEvent(cell_index_,
                                                        creating_cell_,
                                                        GetPointInTime() + GetOffsetNewEvent(),
                                                        GetOffsetNewEvent(),
                                                        growth_type_,
                                                        periodic_continuation_,
                                                        GetElement()));
        }
    }

    void GrowthEvent::
    IncreaseNeighbourCCG(const complex_index cur_cell,
                         const double value,
                         std::map<complex_index, double> *crystal_state_increase,
                         std::map<complex_index, double> *partial_increase) {
        std::vector<complex_index> neighbours;
        std::vector<Eigen::Vector3i> periodic_continuations;
        cur_cell.first->NeighbourCells(cur_cell.second,
                                       false,
                                       true,
                                       &neighbours,
                                       &periodic_continuations);


        const Eigen::Vector3i8 old_periodic_continuation =
                cur_cell.first->GetAdditionalCCGPeriodicContinuation(cur_cell.second);

        for (std::size_t cur = 0; cur < neighbours.size(); ++cur) {
            if (neighbours[cur].first == NULL) {
                continue;
            }

            const Element *cur_el = neighbours[cur].first;
            const double crystal_state =
                    cur_el->GetAdditionalCCGState(neighbours[cur].second);
            // state value should not be higher than the caller
            if (crystal_state >= cur_cell.first->GetAdditionalCCGState(cur_cell.second)) {
                continue;
            }

            // increase the state value
            #pragma omp critical
            {
                if (crystal_state_increase->find(neighbours[cur]) ==
                    crystal_state_increase->end()) {
                    (*crystal_state_increase)[neighbours[cur]] = 0.0;
                }

                // 2-Norm
                const double old_increase_value =
                        (*crystal_state_increase)[neighbours[cur]];
                const double increase_value =
                        sqrt(value * value + old_increase_value * old_increase_value);
                (*crystal_state_increase)[neighbours[cur]] = increase_value;

                // check if this is a newly crystalline cell
                if (crystal_state < 0.0 && crystal_state + increase_value >= 0.0) {
                    if (!cur_el->CellCrystalline(neighbours[cur].second)) {
                        // dont forget periodic continuation, which has to be set manually at ccg
                        // because the growth event is not cell based
                        const Eigen::Vector3i cur_periodic_continuation =
                                periodic_continuations[cur];

                        Eigen::Vector3i new_periodic_continuation;
                        for (int32_t dim = 0; dim < 3; ++dim) {
                            if (cur_periodic_continuation[dim] != 0 &&
                                cur_periodic_continuation[dim] != old_periodic_continuation[dim]) {
                                if (old_periodic_continuation[dim] != 0) {
                                    new_periodic_continuation[dim] = 0;
                                } else {
                                    new_periodic_continuation[dim] = cur_periodic_continuation[dim];
                                }
                            } else {
                                new_periodic_continuation[dim] = old_periodic_continuation[dim];
                            }
                        }

                        auto *non_const_element = const_cast<Element *>(cur_el);
                        non_const_element->
                                ChangeCellToSolid(neighbours[cur].second,
                                                  cur_cell.first->GetSpheruliteID(cur_cell.second),
                                                  GetPointInTime(),
                                                  new_periodic_continuation,
                                                  false);
                    }

                    (*partial_increase)[neighbours[cur]] =
                            (crystal_state + increase_value) / increase_value;
                }
            }
        }
    }

    bool GrowthEvent::SameEvent(const StateEvent *compare_to) const {
        if (compare_to->GetType() != kGrowth) {
            return false;
        }
        if (compare_to->GetElement() != GetElement()) {
            return false;
        }
        if (fabs(compare_to->GetPointInTime() - GetPointInTime()) > 1e-10) {
            return false;
        }

        const auto *growth_event =
                reinterpret_cast<const GrowthEvent *>(compare_to);
        if (growth_event->cell_index_ != cell_index_) {
            return false;
        }
        if (growth_event->creating_cell_ != creating_cell_) {
            return false;
        }

        return true;
    }

    WriteResultEvent::WriteResultEvent(const double point_in_time,
                                       const double offset_new_event,
                                       const Simulation &simulation,
                                       Element *element) :
            StateEvent(kWriteResults,
                       point_in_time,
                       offset_new_event,
                       element),
            simulation_(&simulation) {
    }

    bool WriteResultEvent::SameEvent(const StateEvent *compare_to) const {
        if (compare_to->GetType() != kGrowth) {
            return false;
        }
        if (compare_to->GetElement() != GetElement()) {
            return false;
        }
        if (fabs(compare_to->GetPointInTime() - GetPointInTime()) > 1e-10) {
            return false;
        }
        return true;
    }

    void WriteResultEvent::WriteVtkFile(const Simulation &simulation) const {
        //#pragma omp parallel
        //__itt_suppress_pop();
        std::stringstream ss;
        std::map<const std::string, std::vector<double> > numerical_field_values;
        std::map<const std::string, std::vector<Eigen::Vector3d> > vector_field_values;

        ss << simulation.GetOutputFolder();
        if (!simulation.GetVTKOptions()->GetSubfolder().empty()) {
            // create a subfolder
            ss << simulation.GetVTKOptions()->GetSubfolder();
            std::string subfolder = ss.str();
            std::replace(subfolder.begin(), subfolder.end(), '\\', '/');
            if (subfolder.back() != '/') {
                subfolder.push_back('/');
                ss << "/";
            }
            if (!Directory::folderExists(subfolder)) {
                Directory::create_directory(subfolder);
            }
        }
        const long long time_in_ms = llround(GetPointInTime() * 1e3);

        ss << "Element_" << std::setw(7) << std::setfill('0') << GetElement()->GetElementID();
        ss << "_" << std::setw(6) << std::setfill('0') << time_in_ms << ".vtk";

        VTKFile *vtkfile = new VTKFile(ss.str());

        std::stringstream header;
        header << "v" << MAJOR_VERSION << ". " << MINOR_VERSION << "." << UPDATE_VERSION << " ";
        header << "Time " << time_in_ms << "ms origin";
        header << " " << GetElement()->GetCornerPosition(0)[0];
        header << " " << GetElement()->GetCornerPosition(0)[1];
        header << " " << GetElement()->GetCornerPosition(0)[2];
        if (!GetElement()->GetFullyCrystalline()) {
            header << " model 1";
        } else {
            header << " model 0";
        }
        vtkfile->SetHeaderLine(header.str());

        Eigen::Vector3d corner_position(0.0, 0.0, 0.0);
        if (simulation.GetVTKOptions()->GetAbsoluteCoordinates()) {
            corner_position = GetElement()->GetCornerPosition(0);
        }

        vtkfile->SetStructuredGrid(corner_position,
                                   GetElement()->GetSizeM(),
                                  GetElement()->GetNumCells()[0] + 1,
                                  GetElement()->GetNumCells()[1] + 1,
                                  GetElement()->GetNumCells()[2] + 1);

        // get the fieldinformation for the cells
        size_t num_cells = GetElement()->GetTotalNumCells();

        // allocate memory
        if (simulation.GetVTKOptions()->GetTemperature()) {
            numerical_field_values[Element::VTK_TEMPERATUE].resize(num_cells);
            numerical_field_values[Element::VTK_DISENGAGEMENT_TIME].resize(num_cells);
        }
        if (simulation.GetVTKOptions()->GetPhasestate()) {
            numerical_field_values[Element::VTK_PHASESTATE].resize(num_cells);
        }
        if (simulation.GetVTKOptions()->GetSpheruliteid()) {
            numerical_field_values[Element::VTK_SPHERULITEID].resize(num_cells);
        }
        if (simulation.GetVTKOptions()->GetCrystallizationtime()) {
            numerical_field_values[Element::VTK_CRYSTALLIZATION_TIME].resize(num_cells);
        }
        if (simulation.GetVTKOptions()->GetSpherulitediameter()) {
            numerical_field_values[Element::VTK_SPHERULITE_DIAMETER].resize(num_cells);
        }
        if (simulation.GetVTKOptions()->GetRandomcolor()) {
            numerical_field_values[Element::VTK_RANDOM_COLOR].resize(num_cells);
        }
        if (simulation.GetVTKOptions()->GetAspectratio()) {
            numerical_field_values[Element::VTK_ASPECT_RATIO].resize(num_cells);
        }
        if (simulation.GetVTKOptions()->GetEnthropychange()) {
            numerical_field_values[Element::VTK_ENTHROPY_CHANGE].resize(num_cells);
        }
        if (simulation.GetVTKOptions()->GetInduction()) {
            numerical_field_values[Element::VTK_INDUCTION_ATHERMAL].resize(num_cells);
            numerical_field_values[Element::VTK_INDUCTION_FLOW].resize(num_cells);
            numerical_field_values[Element::VTK_INDUCTION_ATHERMAL_FLOW].resize(num_cells);
            numerical_field_values[Element::VTK_INV_ATHERMAL].resize(num_cells);
        }
        if (simulation.GetVTKOptions()->GetCoolingrate()) {
            numerical_field_values[Element::VTK_COOLING_RATE].resize(num_cells);
        }

        if (simulation.GetVTKOptions()->GetLocation()) {
            vector_field_values[Element::VTK_Location].resize(num_cells);
        }
        if (simulation.GetVTKOptions()->GetVelocity()) {
            vector_field_values[Element::VTK_Velocity].resize(num_cells);
        }

        #pragma omp parallel for
        for (size_t cur_cell = 0; cur_cell < num_cells; ++cur_cell) {
            // get the fields for the cell
            std::map<const std::string, double> numerical_fields;
            std::map<const std::string, Eigen::Vector3d> vector_fields;

            GetElement()->GetNumericalFieldsForCell(cur_cell,
                                                    GetPointInTime(),
                                                    &numerical_fields);
            GetElement()->GetVectorFieldsForCell(cur_cell,
                                                 GetPointInTime(),
                                                 &vector_fields);

            // numerical fields
            for (const auto &it_n:numerical_fields) {
                numerical_field_values[it_n.first][cur_cell] = it_n.second;
            }

            // vector fields
            for (const auto &it_v:vector_fields) {
                vector_field_values[it_v.first][cur_cell] = it_v.second;
            }

        }


        // write scalar numerical fields
        for (const auto &cur_numerical:numerical_field_values) {
            vtkfile->AppendFloatCellField(cur_numerical.second, cur_numerical.first);
        }

        // write vector fields
        for (const auto &cur_vector:vector_field_values) {
            vtkfile->AppendVectorCellField(cur_vector.second, cur_vector.first);
        }

        #pragma omp task default(none) firstprivate(vtkfile)
        vtkfile->WriteToFile(true);
        //__itt_suppress_push(__itt_suppress_memory_errors | __itt_suppress_threading_errors);
    }

    void WriteResultEvent::WriteProbeCSVString(const Probe *probe,
                                               const Simulation &simulation,
                                               std::string *header,
                                               std::string *values) const {
        const char separator = simulation.GetCSVSeparator();

        // numerical values
        std::map<const std::string, double> numerical_values;
        probe->NumericalValues(simulation,
                               GetPointInTime(),
                               &numerical_values);

        // convert the numerical values to a string
        std::stringstream header_stream;
        std::stringstream value_stream;

        // point in time
        header_stream << "Time [s]" << separator;
        value_stream << GetPointInTime() << separator;
        auto it_n = numerical_values.begin();
        while (it_n != numerical_values.end()) {
            header_stream << it_n->first;
            value_stream << it_n->second;
            ++it_n;
            if (it_n != numerical_values.end()) {
                header_stream << separator;
                value_stream << separator;
            }
        }

        // vector field values
        std::map<const std::string, Eigen::Vector3d> vector_values;
        probe->VectorValues(simulation,
                            GetPointInTime(),
                            &vector_values);
        if (!vector_values.empty()) {
            header_stream << separator;
            value_stream << separator;

            auto it_v = vector_values.begin();
            while (it_v != vector_values.end()) {
                header_stream << it_v->first << " (x)" << separator;
                header_stream << it_v->first << " (y)" << separator;
                header_stream << it_v->first << " (z)";
                value_stream << it_v->second[0] << separator;
                value_stream << it_v->second[1] << separator;
                value_stream << it_v->second[2];
                ++it_v;
                if (it_v != vector_values.end()) {
                    header_stream << separator;
                    value_stream << separator;
                }
            }
        }
        *header = header_stream.str();
        *values = value_stream.str();
    }

    void WriteResultEvent::WriteProbeCSVFile(const std::string &filename,
                                             const std::string &header,
                                             const std::string &values) const {
        FILE *file;
        if (Directory::file_exists(filename)) {
            file = fopen(filename.c_str(), "a");
            if (file == nullptr) {
                printf("  -> Error: Could not open '%s'\n", filename.c_str());
                return;
            }
        } else {
            file = fopen(filename.c_str(), "w");
            if (file == nullptr) {
                printf("  -> Error: Could not create '%s'\n", filename.c_str());
                return;
            }
            fprintf(file, "%s\n", header.c_str());
        }
        fprintf(file, "%s\n", values.c_str());
        fclose(file);
    }

    void WriteResultEvent::Execute() {
        // write a vtk file of the current element state
        if (simulation_->GetVTKOptions()->GetSaveVTKFiles()) {
            if (GetOffsetNewEvent() < 0.0 ||
                simulation_->GetVTKOptions()->GetAllResultPoints()) {
                WriteVtkFile(*simulation_);
            }
        }

        // probes for the elements and cells
        std::stringstream probe_element_stream[2];
        std::stringstream probe_cell_stream[2];
        std::size_t total_num_probes = simulation_->GetNumProbes();
        for (std::size_t cur_probe = 0; cur_probe < total_num_probes; ++cur_probe) {
            if (simulation_->GetProbe(cur_probe)->GetElementID() ==
                GetElement()->GetElementID()) {
                std::string new_values, new_headings;
                WriteProbeCSVString(simulation_->GetProbe(cur_probe),
                                    *simulation_,
                                    &new_headings,
                                    &new_values);

                const std::string name = simulation_->GetProbe(cur_probe)->GetName();
                if (simulation_->GetProbe(cur_probe)->GetType() == Probe::kElement) {
                    probe_element_stream[0] << name << simulation_->GetCSVSeparator();
                    probe_element_stream[0] << new_headings;
                    probe_element_stream[1] << simulation_->GetCSVSeparator();
                    probe_element_stream[1] << new_values;
                    if (cur_probe != total_num_probes - 1) {
                        probe_element_stream[0] << simulation_->GetCSVSeparator();
                        probe_element_stream[1] << simulation_->GetCSVSeparator();
                    }
                } else {
                    probe_cell_stream[0] << name << simulation_->GetCSVSeparator();
                    probe_cell_stream[0] << new_headings;
                    probe_cell_stream[1] << simulation_->GetCSVSeparator();
                    probe_cell_stream[1] << new_values;
                    if (cur_probe != total_num_probes - 1) {
                        probe_cell_stream[0] << simulation_->GetCSVSeparator();
                        probe_cell_stream[1] << simulation_->GetCSVSeparator();
                    }
                }
            }
        }

        // create the csv-files
        std::stringstream filename;
        if (!probe_element_stream[1].str().empty()) {
            filename << simulation_->GetOutputFolder() << "Element_Probes.csv";
            WriteProbeCSVFile(filename.str(),
                              probe_element_stream[0].str(),
                              probe_element_stream[1].str());
        }
        if (!probe_cell_stream[1].str().empty()) {
            filename.str("");
            filename << simulation_->GetOutputFolder() << "Cell_Probes.csv";
            WriteProbeCSVFile(filename.str(),
                              probe_cell_stream[0].str(),
                              probe_cell_stream[1].str());
        }

        // create the next event
        if (GetOffsetNewEvent() > 0.0) {
            GetElement()->AddStateEvent(
                    new WriteResultEvent(GetPointInTime() + GetOffsetNewEvent(),
                                         GetOffsetNewEvent(),
                                         *simulation_,
                                         GetElement()));
        }
    }

}  // namespace SphaeroSim