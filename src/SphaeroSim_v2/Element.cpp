#include "Element.h"

#ifndef _INCLUDED_STDINT_H_
#include <stdint.h>
#define _INCLUDED_STDINT_H_
#endif

#ifndef _INCLUDED_VECTOR_H_
#include <vector>
#define _INCLUDED_VECTOR_H_
#endif

#ifndef _INCLUDED_MAP_H_
#include <map>
#define _INCLUDED_MAP_H_
#endif

#ifndef _SPHAEROSIM_EIGEN_LIBRARY_H_
#include "EigenLibrary.h"
#endif

#ifndef _SPHAEROSIM_SIMULATION_H_

#include "Simulation.h"

#endif

#ifndef _SPHAEROSIM_EXCEPTION_H_
#include "Exception.h"
#endif

#ifndef _SPHAEROSIM_INTERPOLATION_H_
#include "Interpolation.h"
#endif

#ifndef _SPHAEROSIM_BOUNDARY_CONDITION_H_
#include "BoundaryCondition.h"
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

#ifndef DEINTEGRATOR_H

#include "DEIntegrator.h"

#endif

#ifndef _SPHAEROSIM_MATERIAL_H_
#include "Material.h"
#endif

#ifndef _SPHAEROSIM_HASH_VALUE_H_

#include "HashValue.h"
#include <string>

#endif

#ifdef LIKWID_PERFMON
#include <likwid.h>
#else
#define LIKWID_MARKER_INIT
#define LIKWID_MARKER_THREADINIT
#define LIKWID_MARKER_REGISTER(regionTag)
#define LIKWID_MARKER_START(regionTag)
#define LIKWID_MARKER_STOP(regionTag)
#define LIKWID_MARKER_CLOSE
#endif

namespace SphaeroSim {

    double Element::table_values_sin[fast_table_resolution];
    double Element::table_values_cos[fast_table_resolution];
    double Element::table_values_log[fast_table_resolution];
    bool Element::fast_access_table_initialized = false;

    void Element::InitFastAccessTables() {
        if (!fast_access_table_initialized) {
            for (int32_t iX = 0; iX < Element::fast_table_resolution; ++iX) {
                table_values_cos[iX] = cos(2.0 * M_PI * (double) iX / (double) fast_table_resolution);
                table_values_sin[iX] = sin(2.0 * M_PI * (double) iX / (double) fast_table_resolution);
                table_values_log[iX] = log(1.0 + (double) iX / (double) fast_table_resolution);
            }
            fast_access_table_initialized = true;
        }
    }

    const double Element::TableValue(const FastAccessTables table,
                                     const double value) {
        int iKey;
        double rInput = value;
        if (table == kCos || table == kSin) {
            if (rInput > 2 * M_PI) {
                rInput -= (int) (rInput / (2.0 * M_PI)) * 2.0 * M_PI;
            }
            iKey = (int) (rInput / (2.0 * M_PI) * fast_table_resolution);
            if (table == kCos) {
                return table_values_cos[iKey];
            } else {
                return table_values_sin[iKey];
            }
        } else {
            iKey = (int) ((rInput - 1.0) * (double) fast_table_resolution);
            return table_values_log[iKey];
        }
    }

    Element::Plane::Plane(const Eigen::Vector3d &normal) {
        normal_ = normal;
        distance_ = 0.0;
    }

// returns the a-value for the intersection with a
// straight line g: x = v + a*u
    const double Element::Plane::Intersection(const Eigen::Vector3d &pos,
                                              const Eigen::Vector3d &direction,
                                              bool *parallel) const {
        double val_1 = 0.0, val_2 = 0.0;

        val_1 = normal_.dot(pos);
        val_2 = normal_.dot(direction);

        if (fabs(val_2) < 1e-12) {
            *parallel = true;
            return 0.0;
        }
        return (distance_ - val_1) / val_2;
    }

    Element::AdditionalInformationCCG::
    AdditionalInformationCCG(const int64_t cell_count) {
        crystal_state_.resize(cell_count);
        periodic_continuation_.resize(cell_count);

        for (int64_t cur_cell = 0; cur_cell < cell_count; ++cur_cell) {
            crystal_state_[cur_cell] = -1.0;
            periodic_continuation_[cur_cell] = Eigen::Vector3i8(0, 0, 0);
        }
    }

    Element::Element(const uint64_t element_id,
                     const Eigen::Vector3d minimum_position_m,
                     const double cell_size_m,
                     const Eigen::Vector3s num_cells,
                     const InterpolationType temporal_interpolation_type,
                     const InterpolationType spatial_interpolation_type,
                     const GrowthEvent::Type growth_type,
                     RandomNumberGenerator *rng,
                     Simulation *current_simulation) :
            element_id_(element_id),
            lower_left_front_corner_m_(minimum_position_m),
            cell_size_m_(cell_size_m),
            growth_type_(growth_type) {

        InitFastAccessTables();

        Initialize(num_cells,
                   temporal_interpolation_type,
                   spatial_interpolation_type,
                   rng,
                   current_simulation);
    }

    Element::Element(const uint64_t element_id,
                     const Eigen::Vector3d minimum_position_m,
                     const double cell_size_m,
                     const Eigen::Vector3d size_m,
                     const InterpolationType temporal_interpolation_type,
                     const InterpolationType spatial_interpolation_type,
                     const GrowthEvent::Type growth_type,
                     RandomNumberGenerator *rng,
                     Simulation *current_simulation) :
            element_id_(element_id),
            lower_left_front_corner_m_(minimum_position_m),
            cell_size_m_(cell_size_m),
            growth_type_(growth_type) {
        Initialize(Eigen::Vector3s(static_cast<std::size_t>(size_m[0] / cell_size_m + 0.5),
                                   static_cast<std::size_t>(size_m[1] / cell_size_m + 0.5),
                                   static_cast<std::size_t>(size_m[2] / cell_size_m + 0.5)),
                   temporal_interpolation_type,
                   spatial_interpolation_type,
                   rng,
                   current_simulation);
    }

    Element::~Element() {
        delete temperatures_;
        delete velocities_;
        delete precalculated_temperatures_;
        delete precalculated_velocities_;
        delete precalculated_shear_flags_;
        delete precalculated_elongation_flags_;
        delete additional_information_ccg_;


        RemoveStateEvents(std::numeric_limits<double>::infinity());
    }

    void Element::Initialize(const Eigen::Vector3s num_cells,
                             const InterpolationType temporal_interpolation_type,
                             const InterpolationType spatial_interpolation_type,
                             RandomNumberGenerator *rng,
                             Simulation *current_simulation) {
        rng_ = rng;
        current_simulation_ = current_simulation;
        has_written_final_results_ = false;
        precalculated_velocities_ = nullptr;
        precalculated_temperatures_ = nullptr;
        precalculated_shear_flags_ = nullptr;
        precalculated_elongation_flags_ = nullptr;

        // areal properties
        num_cells_ = num_cells;
        for (std::size_t cur_dim = 0; cur_dim < 3; ++cur_dim) {
            if (num_cells_[cur_dim] == 0) {
                num_cells_[cur_dim] = 1;
            }
        }

        size_m_ =
                Eigen::Vector3d(static_cast<double>(num_cells_[0]) * cell_size_m_,
                                static_cast<double>(num_cells_[1]) * cell_size_m_,
                                static_cast<double>(num_cells_[2]) * cell_size_m_);

        num_crystalline_cells_ = 0;

        // allocate memory
        phasestate_.resize(GetTotalNumCells(), 0.0);
        crystallization_time_.resize(GetTotalNumCells(), 0.0);

        last_nucleation_rates_.resize(GetTotalNumCells());

        temperatures_ =
                new BoundaryCondition<double>(
                        new Interpolation<double>(temporal_interpolation_type),
                        new Interpolation<double>(spatial_interpolation_type));
        velocities_ =
                new BoundaryCondition<Eigen::Vector3d>(
                        new Interpolation<Eigen::Vector3d>(temporal_interpolation_type),
                        new Interpolation<Eigen::Vector3d>(spatial_interpolation_type));

        additional_information_ccg_ = nullptr;
        if (growth_type_ == GrowthEvent::kCCG) {
            additional_information_ccg_ = new AdditionalInformationCCG(GetTotalNumCells());
        }

        PrecalculateIndices();
    }

    const Eigen::Vector3d Element::GetCornerPosition(const std::size_t index) const {
        //                       6--------7
        //                      /|       /|
        //                     2--------3 |                  y
        //                     | |      | |                  |  z
        //                     | 4------|-5                  | /
        //                     |/       |/                   |/
        //                     0--------1                    ------x
        if (index == 0) {
            return lower_left_front_corner_m_;
        } else if (index == 1) {
            return lower_left_front_corner_m_ +
                   Eigen::Vector3d(size_m_[0], 0.0, 0.0);
        } else if (index == 2) {
            return lower_left_front_corner_m_ +
                   Eigen::Vector3d(0.0, size_m_[1], 0.0);
        } else if (index == 3) {
            return lower_left_front_corner_m_ +
                   Eigen::Vector3d(size_m_[0], size_m_[1], 0.0);
        } else if (index == 4) {
            return lower_left_front_corner_m_ +
                   Eigen::Vector3d(0.0, 0.0, size_m_[2]);
        } else if (index == 5) {
            return lower_left_front_corner_m_ +
                   Eigen::Vector3d(size_m_[0], 0.0, size_m_[2]);
        } else if (index == 6) {
            return lower_left_front_corner_m_ +
                   Eigen::Vector3d(0.0, size_m_[1], size_m_[2]);
        } else if (index == 7) {
            return lower_left_front_corner_m_ +
                   Eigen::Vector3d(size_m_[0], size_m_[1], size_m_[2]);
        }

        throw Exception("GetCornerPosition",
                        "Invalid index.");
    }

// returns true if the eight corner velocities describe a pure
// shear flow
    bool Element::
    VelocitiesPureShear(const Eigen::Vector3d corner_velocities[8],
                        const int32_t main_direction,
                        const int32_t shear_direction,
                        const double threshold) const {
        //                       6--------7
        //                      /|       /|
        //                     2--------3 |                  y
        //                     | |      | |                  |  z
        //                     | 4------|-5                  | /
        //                     |/       |/                   |/
        //                     0--------1                    ------x
        int32_t directions_to_check[2];
        int32_t startpoints[4];
        int32_t endpoints[4];

        int32_t counter = 0;
        for (int32_t index = 0; index < 3; ++index) {
            if (index == shear_direction) {
                continue;
            }
            directions_to_check[counter++] = index;
        }

        if (main_direction == 0) {  // x-direction
            int32_t spoints[] = {0, 2, 4, 6};
            memcpy(startpoints, spoints, sizeof(int32_t) * 4);
            int32_t epoints[] = {1, 3, 5, 7};
            memcpy(endpoints, epoints, sizeof(int32_t) * 4);
        } else if (main_direction == 1) {  // y-direction
            int32_t spoints[] = {0, 1, 4, 5};
            memcpy(startpoints, spoints, sizeof(int32_t) * 4);
            int32_t epoints[] = {2, 3, 6, 7};
            memcpy(endpoints, epoints, sizeof(int32_t) * 4);
        } else if (main_direction == 2) {  // z-direction
            int32_t spoints[] = {0, 1, 2, 3};
            memcpy(startpoints, spoints, sizeof(int32_t) * 4);
            int32_t epoints[] = {4, 5, 6, 7};
            memcpy(endpoints, epoints, sizeof(int32_t) * 4);
        }

        // check if no variation in directions_to_check and a siginificant
        // variation in vatiation_direction
        for (int32_t index = 0; index < 4; ++index) {
            const double non_variation_1 =
                    corner_velocities[endpoints[index]][directions_to_check[0]] -
                    corner_velocities[startpoints[index]][directions_to_check[0]];
            const double non_variation_2 =
                    corner_velocities[endpoints[index]][directions_to_check[1]] -
                    corner_velocities[startpoints[index]][directions_to_check[1]];
            const double variation =
                    corner_velocities[endpoints[index]][shear_direction] -
                    corner_velocities[startpoints[index]][shear_direction];

            if (std::abs(non_variation_1) > threshold ||
                std::abs(non_variation_2) > threshold ||
                std::abs(variation) < threshold) {
                return false;
            }  // not a pure flow
        }
        return true;
    }

// returns true if the eight corner velocities describe a pure
// elongational flow
    bool Element::
    VelocitiesPureElongational(const Eigen::Vector3d corner_velocities[8],
                               const int32_t direction,
                               const double threshold) const {
        //                       6--------7
        //                      /|       /|
        //                     2--------3 |                  y
        //                     | |      | |                  |  z
        //                     | 4------|-5                  | /
        //                     |/       |/                   |/
        //                     0--------1                    ------x
        int32_t startpoints[3][4];
        int32_t endpoints[3][4];

        // x-Axis
        int32_t psx[] = {0, 2, 4, 6};
        int32_t pex[] = {1, 3, 5, 7};
        // y-Axis
        int32_t psy[] = {0, 1, 4, 5};
        int32_t pey[] = {2, 3, 6, 7};
        // z-Axis
        int32_t psz[] = {0, 1, 2, 3};
        int32_t pez[] = {4, 5, 6, 7};

        // copy for better access
        memcpy(startpoints[0], psx, sizeof(int32_t) * 4);
        memcpy(startpoints[1], psy, sizeof(int32_t) * 4);
        memcpy(startpoints[2], psz, sizeof(int32_t) * 4);

        memcpy(endpoints[0], pex, sizeof(int32_t) * 4);
        memcpy(endpoints[1], pey, sizeof(int32_t) * 4);
        memcpy(endpoints[2], pez, sizeof(int32_t) * 4);

        // check if an elongational flow exist in direction
        for (unsigned short int cur = 0; cur < 4; ++cur) {
            double variation =
                    corner_velocities[endpoints[direction][cur]][direction] -
                    corner_velocities[startpoints[direction][cur]][direction];
            if (std::abs(variation) < threshold) {
                return false;
            }

            // check if the other two direction satisfy incompressibility
            const int32_t opposite_1 = (direction + 1) % 3;
            double variation_opposite_1 =
                    corner_velocities[endpoints[opposite_1][cur]][opposite_1] -
                    corner_velocities[startpoints[opposite_1][cur]][opposite_1];
            const int32_t opposite_2 = (direction + 2) % 3;
            double variation_opposite_2 =
                    corner_velocities[endpoints[opposite_2][cur]][opposite_2] -
                    corner_velocities[startpoints[opposite_2][cur]][opposite_2];

            // normalize
            variation /= std::abs(GetCornerPosition(startpoints[direction][0])[direction] -
                                  GetCornerPosition(endpoints[direction][0])[direction]);
            variation_opposite_1 /=
                    std::abs(GetCornerPosition(startpoints[opposite_1][0])[opposite_1] -
                             GetCornerPosition(endpoints[opposite_1][0])[opposite_1]);
            variation_opposite_2 /=
                    std::abs(GetCornerPosition(startpoints[opposite_2][0])[opposite_2] -
                             GetCornerPosition(endpoints[opposite_2][0])[opposite_2]);

            if ((std::abs(variation + 2.0 * variation_opposite_1) > threshold) ||
                (std::abs(variation + 2.0 * variation_opposite_2) > threshold)) {
                return false;
            }
        }
        return true;
    }

// precalculate the boundary conditions
    void Element::PrecalculateBoundaryFields() {
        const Simulation::NumericalOptions *options = current_simulation_->GetNumericalOptions();
        std::size_t number_values = options->cache_number_entries() >= 0 ?
                                    options->cache_number_entries() :
                                    temperatures_->NumberValues();

        std::vector<double> point_in_times;
        temperatures_->PointInTimeOfEntries(1e150, &point_in_times);
        //printf("NumberValues:%zu; PointsInTimeOfEntries:%lu\n", number_values, point_in_times.size());
        if (number_values > point_in_times.size()) {
            number_values = point_in_times.size();
        }

        std::vector<bool> neutral_temperature_elements;
        std::vector<bool> neutral_velocity_elements;
        neutral_temperature_elements.resize(number_values);

        const bool consider_flow = current_simulation_->GetCalculateFlowEnergy();
        if (consider_flow) {
            neutral_velocity_elements.resize(number_values);
        }

        if (number_values > 0) {
            // find the neutral elements
//            double start=omp_get_wtime();
            #pragma omp parallel for
            for (size_t cur = 0; cur < number_values; ++cur) {
                // corner values
                double corner_temperatures[8];
                temperatures_->CornerValues(point_in_times[cur], corner_temperatures);

                Eigen::Vector3d corner_velocities[8];
                if (consider_flow) {
                    velocities_->CornerValues(point_in_times[cur], corner_velocities);
                }

                bool neutral_temperature = true;
                bool neutral_velocity = true;
                for (int8_t index = 0; index < 8; ++index) {
                    if (corner_temperatures[index] > 1.0) {
                        neutral_temperature = false;
                    }
                    if (consider_flow) {
                        if (corner_velocities[index].squaredNorm() >
                            options->threshold_pure_flow() * options->threshold_pure_flow()) {
                            neutral_velocity = false;
                        }
                    }
                }
                neutral_temperature_elements[cur] = neutral_temperature;
                if (consider_flow) {
                    neutral_velocity_elements[cur] = neutral_velocity;
                }
            }
//            printf("After first for: %fs\n", omp_get_wtime()-start);
//            start=omp_get_wtime();

            // velocities
            if (consider_flow) {
                precalculated_velocities_ =
                        new PrecalculatedBoundaryCondition<Eigen::Vector3d>(GetTotalNumCells(),
                                                                            number_values,
                                                                            options->timeresolution_cache(),
                                                                            Eigen::Vector3d(0.0, 0.0, 0.0));
                #pragma omp parallel for
                for (std::size_t cur_cell = 0; cur_cell < GetTotalNumCells(); ++cur_cell) {
                    for (std::size_t cur = 0; cur < number_values; ++cur) {
                        const double point_in_time = point_in_times[cur];
                        if (neutral_velocity_elements[cur] == false) {
                            const Eigen::Vector3d velocity = CellVelocity(cur_cell, point_in_time, false);
                            precalculated_velocities_->AssignPrecalculatedValue(cur_cell,
                                                                                cur,
                                                                                point_in_time,
                                                                                velocity);
                        } else {
                            precalculated_velocities_->AssignNeutralValue(cur, point_in_time);
                        }
                    }
                }
            } else {
                precalculated_velocities_ = nullptr;
            }

//            printf("Before precalculation: %fs\n", omp_get_wtime()-start);
//            start=omp_get_wtime();

            // temperatures
            precalculated_temperatures_ =
                    new PrecalculatedBoundaryCondition<double>(GetTotalNumCells(),
                                                               number_values,
                                                               options->timeresolution_cache(),
                                                               0.0);


//            printf("Before second for: %fs\n", omp_get_wtime()-start);
//            start=omp_get_wtime();

            #pragma omp parallel for
            for (size_t cur_cell = 0; cur_cell < GetTotalNumCells(); ++cur_cell) {
                for (std::size_t cur = 0; cur < number_values; ++cur) {
                    const double point_in_time = point_in_times[cur];
                    if (neutral_temperature_elements[cur] == false) {
                        const double temperature = CellTemperature(cur_cell, point_in_time, false);
                        precalculated_temperatures_->AssignPrecalculatedValue(cur_cell,
                                                                              cur,
                                                                              point_in_time,
                                                                              temperature);
                    } else {
                        precalculated_temperatures_->AssignNeutralValue(cur, point_in_time);
                    }
                }
            }

//            printf("After second for: %fs\n", omp_get_wtime()-start);

            if (consider_flow) {
                precalculated_velocities_->CalculateFastAccessIndices();
            }
            precalculated_temperatures_->CalculateFastAccessIndices();
        }

        if (consider_flow) {
            // shear and elongation flags
            precalculated_shear_flags_ = new
                    PrecalculatedBoundaryCondition<bool>(1,
                                                         velocities_->NumberValues(),
                                                         options->timeresolution_cache(),
                                                         false);
            precalculated_elongation_flags_ = new
                    PrecalculatedBoundaryCondition<bool>(1,
                                                         velocities_->NumberValues(),
                                                         options->timeresolution_cache(),
                                                         false);

            const double threshold = options->threshold_pure_flow();
            const size_t number_point_in_time = velocities_->NumberValues();
            for (size_t cur_val = 0; cur_val < number_point_in_time; ++cur_val) {
                Eigen::Vector3d corner_vel[8];
                velocities_->CornerValues(point_in_times[cur_val], corner_vel);

                // check if a velocity exists
                bool velocity_exist = false;
                for (int32_t cur = 0; cur < 8; ++cur) {
                    if (corner_vel[cur].squaredNorm() > threshold * threshold) {
                        velocity_exist = true;
                        break;
                    }
                }

                int32_t shear_counter = 0;
                int32_t elongation_counter = 0;
                if (velocity_exist) {
                    // check shearing
                    for (unsigned short int dir = 0; dir < 3; ++dir) {
                        // shear flow
                        if (VelocitiesPureShear(corner_vel,
                                                dir,
                                                (dir + 1) % 3,
                                                threshold)) {
                            shear_counter += 1;
                        }
                        if (VelocitiesPureShear(corner_vel,
                                                dir,
                                                (dir + 2) % 3,
                                                threshold)) {
                            shear_counter += 1;
                        }
                        // elongational flow
                        if (VelocitiesPureElongational(corner_vel, dir, threshold)) {
                            elongation_counter += 1;
                        }
                    }
                }

                precalculated_elongation_flags_->
                        AssignPrecalculatedValue(0, cur_val, point_in_times[cur_val], !velocity_exist);
                precalculated_shear_flags_->
                        AssignPrecalculatedValue(0, cur_val, point_in_times[cur_val], !velocity_exist);

                if (velocity_exist) {
                    if (shear_counter == 1 && elongation_counter == 0) {
                        // pure shear flow
                        precalculated_shear_flags_->
                                OverwritePrecalculatedValue(0, cur_val, true);
                    } else if (shear_counter == 0 && elongation_counter == 1) {
                        // pure elongational flow
                        precalculated_elongation_flags_->
                                OverwritePrecalculatedValue(0, cur_val, true);
                    }
                }
            }

            precalculated_shear_flags_->CalculateFastAccessIndices();
            precalculated_elongation_flags_->CalculateFastAccessIndices();
        } else {
            precalculated_shear_flags_ = nullptr;
            precalculated_elongation_flags_ = nullptr;
        }
    }

// SetTemperatureField/SetVelocityField: definition of the boundary conditions
// std::pair<double, T>: first entry denotes the point in time
//                       second entry must containt [8] field values
//                       for the element corners
    void Element::
    SetTemperatureField(const std::vector<Element::temperature_point> &field) {
        for (std::size_t cur = 0; cur < field.size(); ++cur)
            temperatures_->AddEntry(field[cur].first, cur, field[cur].second);
    }

    void Element::
    SetVelocityField(const std::vector<Element::velocity_point> &field) {
        for (std::size_t cur = 0; cur < field.size(); ++cur)
            velocities_->AddEntry(field[cur].first, cur, field[cur].second);
    }

// set the neighbouring element information
    void Element::
    SetNeighbourElementInfos(const std::vector<neighbour_info> &infos,
                             const std::map<const uint64_t, Element *> &elements) {
        if (infos.size() != 6) {
            throw Exception("SetNeighbourElementInfos",
                            "The information of the six surrounding elements are required.");
        }
        const uint64_t no_element = UINT64_MAX;
        for (std::size_t cur = 0; cur < 6; ++cur) {
            const uint64_t id = infos[cur].first;
            neighbour_infos_.emplace_back(
                    (id != no_element && infos[cur].second) ? elements.at(id) : nullptr,
                    id != no_element);
        }
    }

// convert the absolut index to its three dimensionale index representation
    void Element::PrecalculateIndices() {
        precalculated_index_vectors_.clear();
        precalculated_index_vectors_.resize(GetTotalNumCells());
        for (std::size_t cur = 0; cur < GetTotalNumCells(); ++cur) {
            Eigen::Vector3s result;
            std::size_t index = cur;
            result[2] = index / (num_cells_[0] * num_cells_[1]);
            index -= result[2] * num_cells_[0] * num_cells_[1];
            result[1] = index / num_cells_[0];
            index -= result[1] * num_cells_[0];
            result[0] = index;
            precalculated_index_vectors_[cur] = result;
        }
    }

// return a cell index from a position
    const Eigen::Vector3s
    Element::CellIndicesFromPosition(const Eigen::Vector3d &position,
                                     const bool force_inside) const {
        Eigen::Vector3s result;
        Eigen::Vector3d coordinate = position - GetCornerPosition(0);

        for (std::size_t cur = 0; cur < 3; ++cur) {
            result[cur] = static_cast<std::size_t>(coordinate[cur] / size_m_[cur]);
            if (force_inside) {
                if (result[cur] >= num_cells_[cur]) {
                    result[cur] = num_cells_[cur] - 1;
                }
            }
        }
        return result;
    }

// return true if the coordinate is within the element boundaries
    const bool Element::
    CoordinateInElement(const Eigen::Vector3d &position) const {
        Eigen::Vector3d min_pos = GetCornerPosition(0);
        Eigen::Vector3d max_pos = min_pos + size_m_;

        for (std::size_t cur = 0; cur < 3; ++cur) {
            if (position[cur] > max_pos[cur] || position[cur] < min_pos[cur]) {
                return false;
            }
        }
        return true;
    }

    const bool Element::
    CoordinateInElement(const Eigen::Vector3s &indices) const {
        for (std::size_t cur = 0; cur < 3; ++cur) {
            if (indices[cur] >= num_cells_[cur]) {
                return false;
            }
        }
        return true;
    }

// return the cell position in the unit meter
    const Eigen::Vector3d
    Element::CellPosition(const std::size_t cell_index,
                          const bool cell_center) const {
        return CellPosition(CellIndices(cell_index), cell_center);
    }

    const Eigen::Vector3d
    Element::CellPosition(const Eigen::Vector3s &index_triple,
                          const bool cell_center) const {
        Eigen::Vector3d result;
        for (std::size_t cur = 0; cur < 3; ++cur) {
            result[cur] = size_m_[cur] / static_cast<double>(num_cells_[cur]);
            result[cur] *= static_cast<double>(index_triple[cur]);
            result[cur] += lower_left_front_corner_m_[cur];
            if (cell_center) {
                result[cur] += cell_size_m_ * 0.5;
            }
        }
        return result;
    }

// return the relative position of a cell in the element
    const Eigen::Vector3d
    Element::RelativeCellLocation(const std::size_t absolut_index) const {
        Eigen::Vector3d cell_pos = CellPosition(absolut_index, true);
        cell_pos -= lower_left_front_corner_m_;
        return Eigen::Vector3d(cell_pos[0] / size_m_[0],
                               cell_pos[1] / size_m_[1],
                               cell_pos[2] / size_m_[2]);
    }

// calculate the temperature for a cell at a specific point in time
    const double Element::
    CellTemperature(const std::size_t cell_index,
                    const double current_time,
                    const bool check_precalculated) const {
        if (!temperatures_->ValuesExist(current_time)) {
            return -1.0;
        }

        if (precalculated_temperatures_ != nullptr && check_precalculated &&
                precalculated_temperatures_->HasPrecalculatedValues(current_time)) {
            return PrecalculatedBoundaryValue<double>(cell_index,
                                                      current_time,
                                                      precalculated_temperatures_,
                                                      temperatures_);
        }

        double neighbour_values[6][8];
        const double *neighbour_values_ptr[6];
        for (std::size_t cur = 0; cur < 6; ++cur) {
            neighbour_values_ptr[cur] = nullptr;
            if (GetNeighbourDataExist(cur)) {
                const BoundaryCondition<double> *temp_field =
                        neighbour_infos_[cur].first->temperatures_;
                if (temp_field->ValuesExist(current_time)) {
                    temp_field->CornerValues(current_time, neighbour_values[cur]);
                    neighbour_values_ptr[cur] = neighbour_values[cur];
                }
            }
        }
        const auto result = BoundaryValue<double>(cell_index,
                                                  current_time,
                                                  temperatures_,
                                                  neighbour_values_ptr);
        return result;
    }

    const double Element::CellCoolingRate(const std::size_t cell_index,
                                          const double current_time) const {
        const double T_prev = CellTemperature(cell_index, current_time - 1e-6, false);
        const double T_next = CellTemperature(cell_index, current_time + 1e-6, false);
        return (T_next - T_prev) / 2e-6;
    }

// calculate the velocity for a cell at a specific point in time
    const Eigen::Vector3d Element::
    CellVelocity(const std::size_t cell_index,
                 const double current_time,
                 const bool check_precalculated) const {
        if (!velocities_->ValuesExist(current_time)) {
            return Eigen::Vector3d(0.0, 0.0, 0.0);
        }
        if (precalculated_velocities_ != nullptr && check_precalculated &&
            precalculated_velocities_->HasPrecalculatedValues(current_time)) {
            return PrecalculatedBoundaryValue<Eigen::Vector3d>(cell_index,
                                                               current_time,
                                                               precalculated_velocities_,
                                                               velocities_);
        }

        Eigen::Vector3d neighbour_values[6][8];
        const Eigen::Vector3d *neighbour_values_ptr[6];
        for (std::size_t cur = 0; cur < 6; ++cur) {
            neighbour_values_ptr[cur] = nullptr;
            if (GetNeighbourDataExist(cur)) {
                const BoundaryCondition<Eigen::Vector3d> *vel_field =
                        neighbour_infos_[cur].first->velocities_;
                if (vel_field->ValuesExist(current_time)) {
                    vel_field->CornerValues(current_time, neighbour_values[cur]);
                    neighbour_values_ptr[cur] = neighbour_values[cur];
                }
            }
        }

        const Eigen::Vector3d result =
                BoundaryValue<Eigen::Vector3d>(cell_index,
                                               current_time,
                                               velocities_,
                                               neighbour_values_ptr);
        return result;
    }

// solve the differential equations for the deformation tensor for a single timestep
    void Element::solveDifferentialEquations(Eigen::Matrix3d &velGradient,
                                             double &dt,
                                             double zero1,
                                             double zero2,
                                             double zero3,
                                             double prev1_1,
                                             double prev2_1,
                                             double prev3_1,
                                             bool fHasTwoSteps,
                                             double prev1_2,
                                             double prev2_2,
                                             double prev3_2,
                                             double &dst1,
                                             double &dst2,
                                             double &dst3) const {
        double A, B, C, D, E;
        double temp;
        double denominator_1, denominator_2, denominator_3;

        if (!fHasTwoSteps) {
            CALCULATION:
            // calculate prefactors
            temp = 1.0 - dt * velGradient(0, 0);
            if (std::abs(temp) < 1e-8) {
                dt *= 0.5;
            }

            denominator_1 = 1.0 / temp;
            A = velGradient(0, 1) * dt * denominator_1;
            B = velGradient(0, 2) * dt * denominator_1;
            C = zero1 * denominator_1;

            temp = 1.0 - dt * velGradient(1, 0) * A - dt * velGradient(1, 1);
            if (std::abs(temp) < 1e-8) {
                dt *= 0.5;
                goto CALCULATION; // restart
            }

            denominator_2 = 1.0 / temp;
            D = (dt * velGradient(1, 0) * B + dt * velGradient(1, 2)) * denominator_2;
            E = (dt * velGradient(1, 0) * C + zero2) * denominator_2;

            temp = 1.0 - dt * velGradient(2, 0) * (A * D + B) - dt * velGradient(2, 1) * D - dt * velGradient(2, 2);
            if (std::abs(temp) < 1e-8) {
                dt *= 0.5;
                goto CALCULATION; // restart
            }
            denominator_3 = 1.0 / temp;

            dst3 = denominator_3 * (dt * velGradient(2, 0) * (A * E + C) + dt * velGradient(2, 1) * E + zero3);
            dst2 = D * dst3 + E;
            dst1 = A * dst2 + B * dst3 + C;
        } else {
            dst3 = 2 * dt * (velGradient(2, 0) * prev1_1 + velGradient(2, 1) * prev2_1 + velGradient(2, 2) * prev3_1) +
                   prev3_2;
            dst2 = 2 * dt * (velGradient(1, 0) * prev1_1 + velGradient(1, 1) * prev2_1 + velGradient(1, 2) * prev3_1) +
                   prev2_2;
            dst1 = 2 * dt * (velGradient(0, 0) * prev1_1 + velGradient(0, 1) * prev2_1 + velGradient(0, 2) * prev3_1) +
                   prev1_2;
        }
    }

// get the deformation tensor for the timeintervall t1 to t2 (not accurate yet. Use integrate over time for accurate calculatione (not supported yet))
// t1 > t2
    Eigen::Matrix3d Element::DeformationTensor(const std::size_t cell_index,
                                               const double time_1,
                                               const double time_2,
                                               bool *consider_pure_elongational) const {
        const Simulation::NumericalOptions *options = current_simulation_->GetNumericalOptions();

        double rTimeStep = options->timestep_fdm();
        bool fFirstIteration = true;
        long lTimeFactor = static_cast<long>(options->timeresolution_fdm());
        long lTime1 = (long) (time_1 * lTimeFactor);
        long lTime2 = (long) (time_2 * lTimeFactor);
        long lCurTime = lTime2;
        long lTimeStep = (long) (rTimeStep * lTimeFactor);
        double rInvTimeFactor = (double) 1.0 / lTimeFactor;

        *consider_pure_elongational = false;

        Eigen::Vector3d cornerVel[8];
        Eigen::Matrix3d velGradTensor;
        Eigen::Matrix3d curTensor, nextTensor;
        Eigen::Matrix3d prevTensor = Eigen::Matrix3d::Zero();

        prevTensor(0, 0) = prevTensor(1, 1) = prevTensor(2, 2) = 1.0; // E(t, t) = 1
        curTensor = prevTensor;
        nextTensor = prevTensor;

        while (lCurTime < lTime1) {
            double rCurTime = (double) lCurTime * rInvTimeFactor;

            // build velocity gradient
            for (unsigned short int iX = 0; iX < 3; ++iX) {
                for (unsigned short int iY = 0; iY < 3; ++iY) {
                    velGradTensor(iX, iY) = CellVelocityGradient(cell_index, iX, iY, rCurTime);
                    if (std::abs(velGradTensor(iX, iY)) < 1e-10) {
                        velGradTensor(iX, iY) = 0.0;
                    }
                }
            }

            solveDifferentialEquations(velGradTensor, rTimeStep, 1.0, 0.0, 0.0, curTensor(0, 0), curTensor(1, 0),
                                       curTensor(2, 0), !fFirstIteration, prevTensor(0, 0), prevTensor(1, 0),
                                       prevTensor(2, 0), nextTensor(0, 0), nextTensor(1, 0), nextTensor(2, 0));
            solveDifferentialEquations(velGradTensor, rTimeStep, 0.0, 1.0, 0.0, curTensor(0, 1), curTensor(1, 1),
                                       curTensor(2, 1), !fFirstIteration, prevTensor(0, 1), prevTensor(1, 1),
                                       prevTensor(2, 1), nextTensor(0, 1), nextTensor(1, 1), nextTensor(2, 1));
            solveDifferentialEquations(velGradTensor, rTimeStep, 0.0, 0.0, 1.0, curTensor(0, 2), curTensor(1, 2),
                                       curTensor(2, 2), !fFirstIteration, prevTensor(0, 2), prevTensor(1, 2),
                                       prevTensor(2, 2), nextTensor(0, 2), nextTensor(1, 2), nextTensor(2, 2));

            // check incompressibility
            double rDet = nextTensor(0, 0) * nextTensor(1, 1) * nextTensor(2, 2) +
                          nextTensor(0, 1) * nextTensor(1, 2) * nextTensor(2, 0) +
                          nextTensor(0, 2) * nextTensor(1, 0) * nextTensor(2, 1) -
                          nextTensor(0, 2) * nextTensor(1, 1) * nextTensor(2, 0) -
                          nextTensor(0, 1) * nextTensor(1, 0) * nextTensor(2, 2) -
                          nextTensor(0, 0) * nextTensor(1, 2) * nextTensor(2, 1);
            if (std::abs(rDet - 1.0) > options->threshold_incompressibility()) {
                // try again with smaller timestep
                lTimeStep /= 10;
                if (lTimeStep == 0) {
                    if (nextTensor(0, 0) < 0.0 || nextTensor(1, 1) < 0.0 || nextTensor(2, 2) < 0.0) {
                        // instable -> consider as purely elongational flow
                        *consider_pure_elongational = true;
                        break;
                    }
                    // force inkompressible deformation tensor
                    //makeIncompressible(curTensor);
                } else {
                    continue;
                }
            }
            fFirstIteration = false;

            prevTensor = curTensor;
            curTensor = nextTensor;
            if (curTensor(0, 0) > 1e150 || curTensor(0, 0) < -1e150 ||
                curTensor(1, 1) > 1e150 || curTensor(1, 1) < -1e150 ||
                curTensor(2, 2) > 1e150 || curTensor(2, 2) < -1e150) {
                // numerically instable because of double limitation
                // consider the flow as a purely elongational flow
                *consider_pure_elongational = true;
                break;
            }
            lCurTime += lTimeStep;
        }
        return curTensor;
    }

// calculate the dimensionless induction time
    const double Element::DimensionlessInduction(const double point_in_time,
                                                 const std::size_t cell_index,
                                                 const bool flow,
                                                 const bool thermal) const {
        if (flow && !thermal) {
            // quiescent
            const double quiescent_rate =
                    current_simulation_->GetNucleationModel()->Rate(point_in_time,
                                                                    cell_index,
                                                                    true,
                                                                    false,
                                                                    this);

            // flowinduced
            const double flow_rate =
                    current_simulation_->GetNucleationModel()->Rate(point_in_time,
                                                                    cell_index,
                                                                    false,
                                                                    false,
                                                                    this);

            if (std::abs(flow_rate) == 0.0) {
                return -1.0;
            }
            return quiescent_rate / flow_rate;
        } else if (!flow && thermal) {
            // thermal
            const double thermal_rate =
                    current_simulation_->GetNucleationModel()->Rate(point_in_time,
                                                                    cell_index,
                                                                    true,
                                                                    true,
                                                                    this);
            // athermal
            const double athermal_rate =
                    current_simulation_->GetNucleationModel()->Rate(point_in_time,
                                                                    cell_index,
                                                                    true,
                                                                    false,
                                                                    this);
            if (std::abs(athermal_rate) == 0.0) {
                return -1.0;
            }
            return thermal_rate / athermal_rate;
        } else if (flow && thermal) {
            // thermal
            const double thermal_rate =
                    current_simulation_->GetNucleationModel()->Rate(point_in_time,
                                                                    cell_index,
                                                                    false,
                                                                    true,
                                                                    this);
            // athermal
            const double athermal_rate =
                    current_simulation_->GetNucleationModel()->Rate(point_in_time,
                                                                    cell_index,
                                                                    false,
                                                                    false,
                                                                    this);
            if (std::abs(athermal_rate) == 0.0) {
                return -1.0;
            }
            return thermal_rate / athermal_rate;
        } else {
            return -1.0;
        }
    }

// check if the flow can be represented as pure shearing or
// elongational flow
    bool Element::FlowPureShearElongational(const double start_time,
                                            const double end_time,
                                            bool *consider_elongational,
                                            bool *consider_shearing,
                                            bool *no_flow) const {
        // start index
        std::size_t start_index_1;
        std::size_t start_index_2;
        precalculated_shear_flags_->Indices(start_time,
                                            &start_index_1,
                                            &start_index_2);
        const std::size_t start_index = start_index_1 < start_index_2 ?
                                        start_index_1 :
                                        start_index_2;
        if (start_index >= precalculated_shear_flags_->NumberValuesPerCell()) {
            throw Exception("Internal error",
                            "Wrong number of precalculated flags");
        }

        // end index
        std::size_t end_index_1;
        std::size_t end_index_2;
        precalculated_shear_flags_->Indices(end_time,
                                            &end_index_1,
                                            &end_index_2);
        std::size_t end_index = end_index_1 > end_index_2 ?
                                end_index_1 :
                                end_index_2;
        if (end_index >= precalculated_shear_flags_->NumberValuesPerCell()) {
            end_index = start_index;
        }

        *consider_elongational = true;
        *consider_shearing = true;
        *no_flow = true;
        for (std::size_t cur = start_index; cur <= end_index; ++cur) {
            if (!precalculated_elongation_flags_->Value(0, cur)) {
                *consider_elongational = false;
            }
            if (!precalculated_shear_flags_->Value(0, cur)) {
                *consider_shearing = false;
            }
            if (!precalculated_shear_flags_->NeutralValue(cur)) {
                *no_flow = false;
            }
            if (!*consider_shearing &&
                !*consider_elongational &&
                !*no_flow) {
                return false;
            }
        }
        return true;
    }

// calculat the deformation of a cell in a pure shearflow
    const double Element::CellPointShearing(const std::size_t cell_index,
                                            const double start_time,
                                            const double end_time) const {
        // find the shearing direction
        double max_gradient = -1e10;
        int32_t main_direction = -1;
        int32_t shear_direction = -1;
        for (int32_t time_index = 0; time_index < 2; ++time_index) {
            double cur_time;
            if (time_index == 0) {
                cur_time = start_time;
            } else {
                cur_time = end_time;
            }
            for (unsigned short int direction = 0; direction < 3; ++direction) {
                double gradient = CellVelocityGradient(cell_index,
                                                       direction,
                                                       (direction + 1) % 3,
                                                       cur_time);
                if (std::abs(gradient) > max_gradient) {
                    main_direction = direction;
                    shear_direction = (direction + 1) % 3;
                    max_gradient = std::abs(gradient);
                }
                gradient = CellVelocityGradient(cell_index,
                                                direction,
                                                (direction + 2) % 3,
                                                cur_time);
                if (std::abs(gradient) > max_gradient) {
                    main_direction = direction;
                    shear_direction = (direction + 2) % 3;
                    max_gradient = (gradient);
                }
            }
        }

        class IntegralShearRate {
        public:
            IntegralShearRate(const std::size_t cell_index,
                              const int32_t main_direction,
                              const int32_t shear_direction,
                              const Element *el) {
                cell_index_ = cell_index;
                main_direction_ = main_direction;
                shear_direction_ = shear_direction;
                el_ = el;
            }

            double operator()(const double point_in_time) const {
                const double shearrate =
                        el_->CellVelocityGradient(cell_index_,
                                                  main_direction_,
                                                  shear_direction_,
                                                  point_in_time);
                return std::abs(shearrate);
            }

        private:
            std::size_t cell_index_;
            int32_t main_direction_;
            int32_t shear_direction_;
            const Element *el_;
        } shear_integral(cell_index, main_direction, shear_direction, this);

        // integrate
        const double integral_error = current_simulation_->GetIntegralError();
        const double shear_deformation =
                DEIntegrator<IntegralShearRate>::Integrate(shear_integral,
                                                           start_time,
                                                           end_time,
                                                           integral_error);

        return shear_deformation;
    }

// calculate the ln  of the point elongation in case of a purely elongational flow
    const double Element::CellPointElongationLn(const std::size_t cell_index,
                                                const double start_time,
                                                const double end_time) const {

        class IntegralElongationRate {
        public:
            IntegralElongationRate(const std::size_t cell_index,
                                   const std::size_t component,
                                   const Element *el) :
                    cell_index_(cell_index),
                    component_(component),
                    el_(el) {
            }

            double operator()(const double point_in_time) const {
                const double elongational_rate =
                        el_->CellVelocityGradient(cell_index_,
                                                  component_,
                                                  component_,
                                                  point_in_time);
                return std::abs(elongational_rate);
            }

        private:
            const std::size_t cell_index_;
            const std::size_t component_;
            const Element *el_;

            IntegralElongationRate &operator=(const IntegralElongationRate &) {
                return *this;
            }
        };
        const double integral_error = current_simulation_->GetIntegralError();

        // compute the elongation in x, y and z direction
        double elongation[3];
        for (std::size_t direction = 0; direction < 3; ++direction)
            elongation[direction] =
                    DEIntegrator<IntegralElongationRate>::
                    Integrate(IntegralElongationRate(cell_index, direction, this),
                              start_time,
                              end_time,
                              integral_error);

        if (elongation[0] > elongation[1]) {
            if (elongation[0] > elongation[2]) {
                return elongation[0];
            }
        } else {
            if (elongation[1] > elongation[2]) {
                return elongation[1];
            }
        }
        return elongation[2];
    }

// calcualte the flow integral using the reptations theory
    const double Element::CellFlowIntegral(const std::size_t cell_index,
                                           const double current_time) const {
        if (!GetCurrentSimulation()->GetCalculateFlowEnergy()) {
            return 0.0;
        }
        if (current_time < 1e-10) {
            return 0.0;
        }

        // integrate(ln(|E*u|), |u| = 1) (theta integration)
        class InnerIntegrationTheta {
        public:
            InnerIntegrationTheta(const Eigen::Matrix3d &deformation_tensor,
                                  const double phi) :
                    deformation_tensor_(deformation_tensor),
                    phi_(phi) {
            }

            InnerIntegrationTheta(const InnerIntegrationTheta &src) :
                    deformation_tensor_(src.deformation_tensor_),
                    phi_(src.phi_) {
            }

            double operator()(const double theta) const {
                //      sin(Theta)cos(Phi)
                // u =  sin(Theta)sin(Phi)
                //          cos(Theta)
                const Eigen::Vector3d u(
                        Element::TableValue(kSin, theta) * Element::TableValue(kCos, phi_),
                        Element::TableValue(kSin, theta) * Element::TableValue(kSin, phi_),
                        Element::TableValue(kCos, theta));

                const Eigen::Vector3d proj = deformation_tensor_ * u;
                const double squared_proj = proj.squaredNorm();
                return 0.5 * fast_log(squared_proj) * Element::TableValue(kSin, theta);
            }

        private:
            inline double fast_log(double val) const {
                uint64_t raw_val = *(uint64_t *) &val;
                int exp_2 = (int) ((raw_val >> 52l) & 0x7FF) - 1023;
                uint64_t mant_int = raw_val & 0xFFFFFFFFFFFFF;
                return exp_2 * 0.69314718 + Element::TableValue(kLog, 1.0 + (double) mant_int * 2.22044604925e-16);
            }

            const double phi_;
            const Eigen::Matrix3d deformation_tensor_;

            InnerIntegrationTheta &operator=(const InnerIntegrationTheta &) {
                return *this;
            }
        };

        // integrate(ln(|E*u|), |u| = 1) (phi integration)
        class InnerIntegration {
        public:
            InnerIntegration(const Eigen::Matrix3d &deformation_tensor,
                             const double integral_error) :
                    integral_error_(integral_error),
                    deformation_tensor_(deformation_tensor) {
            }

            InnerIntegration(const InnerIntegration &src) :
                    integral_error_(src.integral_error_),
                    deformation_tensor_(src.deformation_tensor_) {
            }

            double operator()(const double phi) const {
                const double inner_integral =
                        DEIntegrator<InnerIntegrationTheta>::
                        Integrate(InnerIntegrationTheta(deformation_tensor_, phi),
                                  0.0,
                                  M_PI,
                                  integral_error_);
                return inner_integral;
            }

        private:
            const double integral_error_;
            const Eigen::Matrix3d deformation_tensor_;

            InnerIntegration &operator=(const InnerIntegration &) {
                return *this;
            }
        };
        // integrate(memory_function * innner_integral, 0, t)
        class OuterIntegration {
        public:
            OuterIntegration(const std::size_t cell_index,
                             const double scaling,
                             const double end_time,
                             const double integral_error,
                             const Element *element) :
                    cell_index_(cell_index),
                    scaling_(scaling),
                    end_time_(end_time),
                    integral_error_(integral_error),
                    el_(element) {
            }

            OuterIntegration(const OuterIntegration &src) = default;

            double operator()(const double point_in_time) const {
                const double memory_value = MemoryFunction(end_time_,
                                                           point_in_time);
                bool consider_elongational = false;
                bool consider_shearing = false;
                bool no_flow = false;

                // check if the flow can be represented as pure shearing or
                // elongational flow
                el_->FlowPureShearElongational(point_in_time,
                                               end_time_,
                                               &consider_elongational,
                                               &consider_shearing,
                                               &no_flow);
                if (no_flow) {
                    return 0.0;
                }

                const Simulation::NumericalOptions *options =
                        el_->GetCurrentSimulation()->GetNumericalOptions();
                if (options->consider_pure_elongational()) {
                    consider_elongational = true;
                }
                if (options->consider_pure_shear()) {
                    consider_shearing = true;
                }

                // if no pure shearing or elongational the deformation tensor has to
                // be calculated via finite difference methode
                Eigen::Matrix3d deformation_tensor;
                if (!consider_elongational && !consider_shearing) {
                    deformation_tensor = el_->DeformationTensor(cell_index_,
                                                                end_time_,
                                                                point_in_time,
                                                                &consider_elongational);
                }

                double inner_integral = 0.0;
                if (consider_elongational) {
                    // consider the flow as purely elongational
                    // for detailed information about the formular see
                    // "Microrheological Modeling of Flow-Induced Crystallization"
                    // from Coppola and Grizzuti
                    const double elongation_ln =
                            el_->CellPointElongationLn(cell_index_, point_in_time, end_time_);

                    if (elongation_ln > 20.0) {
                        inner_integral = elongation_ln - 1.0;  // atan(x)/x can be neglected
                    } else {
                        const double lambda = exp(elongation_ln);
                        if (std::abs(lambda - 1.0) < 1e-9) {
                            inner_integral = 0.0;
                        } else {
                            const double lambda_sqrt = lambda < 1.0 ? 0.0 : sqrt(lambda * lambda * lambda - 1);
                            inner_integral = elongation_ln +
                                             atan(lambda_sqrt) / lambda_sqrt -
                                             1.0;
                        }
                    }
                } else if (consider_shearing) {
                    // consider the flow as purely elongational
                    // for detailed information about the formular see
                    // "Microrheological Modeling of Flow-Induced Crystallization"
                    // from Coppola and Grizzuti
                    const double shearing_deformation =
                            el_->CellPointShearing(cell_index_, point_in_time, end_time_);

                    class ShearingOperator {
                    public:
                        ShearingOperator(const double shearing_deformation) {
                            gamma_ = shearing_deformation;
                            gamma_2_ = shearing_deformation *
                                       shearing_deformation;
                        }

                        inline double operator()(const double x) const {
                            const double x_2 = x * x;
                            const double sqrt_argument =
                                    x_2 * x_2 * (gamma_2_ * gamma_2_ + 4.0 * gamma_2_) -
                                    2.0 * gamma_2_ * x_2 +
                                    1.0;
                            const double sqrt_term = std::sqrt(sqrt_argument);
                            return fast_log(0.5 * (1.0 + gamma_2_ * x_2 + sqrt_term));
                        }

                    private:
                        inline double fast_log(double val) const {
                            uint64_t raw_val = *(uint64_t *) &val;
                            int exp_2 = (int) ((raw_val >> 52l) & 0x7FF) - 1023;
                            uint64_t mant_int = raw_val & 0xFFFFFFFFFFFFF;
                            return exp_2 * 0.69314718 +
                                   Element::TableValue(kLog, 1.0 + (double) mant_int * 2.22044604925e-16);
                        }

                        inline float fast_sqrt(float x) const {
                            unsigned int i = *(unsigned int *) &x;
                            // adjust bias
                            i += 127 << 23;
                            // approximation of square root
                            i >>= 1;
                            return *(float *) &i;
                        }

                        double gamma_;
                        double gamma_2_;
                    } shearing_operator(shearing_deformation);

                    if (std::abs(shearing_deformation) < 1e-9) {
                        inner_integral = 0.0;
                    } else {
                        inner_integral = 0.5 *
                                         DEIntegrator<ShearingOperator>::Integrate(shearing_operator,
                                                                                   0.0,
                                                                                   1.0,
                                                                                   integral_error_);
                    }
                } else {
                    inner_integral = 1.0 / (4.0 * M_PI) *
                                     DEIntegrator<InnerIntegration>::
                                     Integrate(InnerIntegration(deformation_tensor,
                                                                integral_error_),
                                               0.0,
                                               2.0 * M_PI,
                                               integral_error_);
                }
                return scaling_ * memory_value * inner_integral;
            }

        private:
            const double MemoryFunction(const double time_1,
                                        const double time_2) const {
                const double temperature = el_->CellTemperature(cell_index_,
                                                                time_2);
                const double relaxation_time =
                        el_->GetCurrentSimulation()->GetMaterial()->
                                GetDisengagementTime(temperature);
                // d/dt(8/Pi²*sum(1/p²exp(-p²(t-t')/Td), p:odd))
                int iP = 1;
                int iCounter = 1;
                double rResult = 0.0;
                double rValue = 0.0;
                for (;;) {
                    rValue = exp(-iP * iP * (time_1 - time_2) / relaxation_time);
                    if (rValue < 1e-6 || std::abs(time_1 - time_2) < 1e-8) {
                        break;
                    }
                    ++iCounter;
                    iP = 2 * iCounter + 1;
                    rResult += rValue;

                    if (iCounter >= 100) {
                        break;
                    }
                }
                return 8.0 * rResult / (M_PI * M_PI * relaxation_time);
            }

            inline double fast_log(double val) const {
                uint64_t raw_val = *(uint64_t *) &val;
                int exp_2 = (int) ((raw_val >> 52l) & 0x7FF) - 1023;
                uint64_t mant_int = raw_val & 0xFFFFFFFFFFFFF;
                return exp_2 * 0.69314718 + Element::TableValue(kLog, 1.0 + (double) mant_int * 2.22044604925e-16);
            }

            const double scaling_;
            const double end_time_;
            const double integral_error_;
            const std::size_t cell_index_;
            const Element *el_;

            OuterIntegration &operator=(const OuterIntegration &) {
                return *this;
            }
        };

        const double integral_error = current_simulation_->GetIntegralError();
        const double start_point = temperatures_->FirstPointInTime();

        double integral =
                DEIntegrator<OuterIntegration>::
                Integrate(OuterIntegration(cell_index, 1.0, current_time, integral_error, this),
                          start_point,
                          current_time,
                          integral_error);

        if (std::abs(integral) < integral_error) {
            integral = 0.0;
        }
        return integral;
    }

// calcualte the flow energy using the reptations theory
    const double Element::CellFlowEnergy(const std::size_t cell_index,
                                         const double current_time,
                                         double *enthropy) const {
        if (!GetCurrentSimulation()->GetCalculateFlowEnergy()) {
            return 0.0;
        }
        if (current_time < 1e-10) {
            return 0.0;
        }
        /*
        union {
          int64_t key_value;
          int32_t single_values[2];
        } key;

        key.single_values[0] = static_cast<int32_t>(cell_index);
        key.single_values[1] =
          static_cast<int32_t>(current_time *
                               current_simulation_->GetNumericalOptions()->timeresolution_cache());
        */
        const double integral = CellFlowIntegral(cell_index,
                                                 current_time);

        const Material *material = GetCurrentSimulation()->GetMaterial();
        const double rDensity = material->GetDensity();
        const double rMolecularWeightEntangle = material->GetMolWeightEnt() / 1000.0;

        const double rAvogadro = 6e23;
        const double rBoltzmann = 1.380650424e-23;
        const double rConcentration = rDensity * rAvogadro / rMolecularWeightEntangle;

        const double temperature = CellTemperature(cell_index, current_time);
        const double prefactor = 3.0 * rConcentration * rBoltzmann * temperature;

        if (enthropy != nullptr) {
            *enthropy = -integral * 3.0 * rBoltzmann;
        }

        return prefactor * integral;
    }

    const double Element::RelativeCrystallizationDegree() const {
        const std::size_t total_num_cells = GetTotalNumCells();
        return static_cast<double>(num_crystalline_cells_) /
               static_cast<double>(total_num_cells);
    }

// get a complete set of information for the current state of a cell
    void Element::
    GetNumericalFieldsForCell(const std::size_t &index,
                              const double &current_point_in_time,
                              std::map<const std::string,
                                      double> *field_values) const {
        // temperature
        if (current_simulation_->GetVTKOptions()->GetTemperature()) {
            (*field_values)[VTK_TEMPERATUE] =
                    CellTemperature(index, current_point_in_time) - 273.15;

            // disengagement time
            (*field_values)[VTK_DISENGAGEMENT_TIME] =
                    current_simulation_->GetMaterial()->GetDisengagementTime(
                            (*field_values)[VTK_TEMPERATUE] + 273.15);
        }

        // phase state
        if (current_simulation_->GetVTKOptions()->GetPhasestate()) {
            (*field_values)[VTK_PHASESTATE] = phasestate_[index];
        }

        // spherulite id
        if (current_simulation_->GetVTKOptions()->GetSpheruliteid()) {
            if (CellCrystalline(index)) {
                //new option checking boundry layer
                if (CellInBorder(index)) {
                    if (current_simulation_->GetVTKOptions()->Getcheckboundarylayer()) {
                        (*field_values)[VTK_SPHERULITEID] = 1.0;
                    } else {
                        (*field_values)[VTK_SPHERULITEID] = static_cast<double>(GetSpheruliteID(index));
                    }
                } else {
                    (*field_values)[VTK_SPHERULITEID] = static_cast<double>(GetSpheruliteID(index));
                }
            } else {
                (*field_values)[VTK_SPHERULITEID] = 0.0;
            }
        }

//		  if (current_simulation_->GetVTKOptions()->GetSpheruliteid())
//		  {
//			  if (CellInBorder(index))
//			  {
//				  if (current_simulation_->GetVTKOptions()->Getcheckboundarylayer())
//					  (*field_values)[VTK_SPHERULITEID] = 1.0;
//				  else
//					  (*field_values)[VTK_SPHERULITEID] = static_cast<double>(GetSpheruliteID(index));
//			  }
//
//			  if (current_simulation_->GetVTKOptions()->Getcheckboundarylayer())
//			  {
//				  if (CellInBorder(index))
//					  (*field_values)[VTK_SPHERULITEID] = 1.0;
//			  }
//			  else
//				  (*field_values)[VTK_SPHERULITEID] = static_cast<double>(GetSpheruliteID(index));
//		  }
//		  
//	  }
//	  else {
//     (*field_values)[VTK_SPHERULITEID] = 0.0;
//    }
//  }


        // crystallization time
        if (current_simulation_->GetVTKOptions()->GetCrystallizationtime()) {
            (*field_values)[VTK_CRYSTALLIZATION_TIME] = crystallization_time_[index];
        }

        const SpheruliteInfo *spherulite = nullptr;
        if (CellCrystalline(index)) {
            spherulite = current_simulation_->getSpherulite_infos()->at(GetSpheruliteID(index));
        }

        // spherulite diameter
        if (current_simulation_->GetVTKOptions()->GetSpherulitediameter()) {
            if (CellCrystalline(index)) {
                (*field_values)[VTK_SPHERULITE_DIAMETER] =
                        current_simulation_->SpheruliteDiameterNoDataracePossible(spherulite) * 1e6;
            } else {
                (*field_values)[VTK_SPHERULITE_DIAMETER] = 0.0;
            }
        }

        // random color per spherulite
        if (current_simulation_->GetVTKOptions()->GetRandomcolor()) {
            if (CellCrystalline(index)) {
                (*field_values)[VTK_RANDOM_COLOR] = current_simulation_->SpheruliteRandomColor(spherulite);
            } else {
                (*field_values)[VTK_RANDOM_COLOR] = -1.0;
            }
        }
        // spherulite aspect ratio
        if (current_simulation_->GetVTKOptions()->GetAspectratio()) {
            if (CellCrystalline(index)) {
                (*field_values)[VTK_ASPECT_RATIO] =
                        current_simulation_->SpheruliteAspectRatioNoDatraracePossible(spherulite);
            } else {
                (*field_values)[VTK_ASPECT_RATIO] = 0.0;
            }
        }

        // flow integral
        if (current_simulation_->GetVTKOptions()->GetEnthropychange()) {
            const double rBoltzmann = 1.380650424e-23;
            (*field_values)[VTK_ENTHROPY_CHANGE] =
                    -3.0 * rBoltzmann * CellFlowIntegral(index, current_point_in_time);
        }

        if (current_simulation_->GetVTKOptions()->GetInduction()) {
            // dimensionless induction time for flow
            (*field_values)[VTK_INDUCTION_FLOW] =
                    DimensionlessInduction(current_point_in_time, index, true, false);

            // dimensionless induction time for athermal nucleation
            (*field_values)[VTK_INDUCTION_ATHERMAL] =
                    DimensionlessInduction(current_point_in_time, index, false, true);

            // dimensionless induction time for athermal+flow nucleation
            (*field_values)[VTK_INDUCTION_ATHERMAL_FLOW] =
                    DimensionlessInduction(current_point_in_time, index, true, true);

            // inv. athermal nucleation
            if ((*field_values)[VTK_INDUCTION_ATHERMAL] > 0.0) {
                (*field_values)[VTK_INV_ATHERMAL] = 1.0 / (*field_values)[VTK_INDUCTION_ATHERMAL];
            } else {
                (*field_values)[VTK_INV_ATHERMAL] = 0.0;
            }
        }

        // cooling rate
        if (current_simulation_->GetVTKOptions()->GetCoolingrate()) {
            (*field_values)[VTK_COOLING_RATE] =
                    CellCoolingRate(index, current_point_in_time);
        }
    }

    void Element::
    GetVectorFieldsForCell(const std::size_t index,
                           const double current_point_in_time,
                           std::map<const std::string,
                                   Eigen::Vector3d> *field_values) const {
        // location
        if (current_simulation_->GetVTKOptions()->GetLocation()) {
            (*field_values)[VTK_Location] = CellPosition(index, true) * 1000.0;
        }

        // velocity
        if (current_simulation_->GetVTKOptions()->GetVelocity()) {
            (*field_values)[VTK_Velocity] =
                    CellVelocity(index, current_point_in_time) * 1000.0;
        }
    }

    void Element::
    GetNumericalFieldsForElement(const double current_point_in_time,
                                 const Eigen::Vector3s starting_cell,
                                 const Eigen::Vector3s num_cells,
                                 std::map<const std::string,
                                         double> *field_values) const {
        std::size_t cells_in_area = 0;
        std::size_t crystalline_cells_in_area = 0;
        double summed_temperature = 0.0;

        for (std::size_t cur_x = 0; cur_x < num_cells[0]; ++cur_x) {
            for (std::size_t cur_y = 0; cur_y < num_cells[1]; ++cur_y) {
                for (std::size_t cur_z = 0; cur_z < num_cells[2]; ++cur_z) {
                    Eigen::Vector3s cur_index = starting_cell +
                                                Eigen::Vector3s(cur_x, cur_y, cur_z);
                    if (cur_index[0] < num_cells_[0] &&
                        cur_index[1] < num_cells_[1] &&
                        cur_index[2] < num_cells_[2]) {

                        // valid index in the element
                        summed_temperature += CellTemperature(AbsoluteCellIndex(cur_index),
                                                              current_point_in_time);

                        // check if crystalline
                        if (CellCrystalline(AbsoluteCellIndex(cur_index))) {
                            ++crystalline_cells_in_area;
                        }
                        ++cells_in_area;
                    }
                }
            }
        }

        if (cells_in_area > 0) {
            std::size_t num_spherulites = 0;
            const std::pair<const double, const double> average_diameter =
                    current_simulation_->AverageSpheruliteDiameterForElement(this,
                                                                             starting_cell,
                                                                             num_cells,
                                                                             &num_spherulites);

            // insert the fields
            (*field_values)["Average Temp. [C]"] =
                    summed_temperature / cells_in_area - 273.15;
            (*field_values)["Cryst. degree [-]"] =
                    static_cast<double>(crystalline_cells_in_area) /
                    static_cast<double>(cells_in_area);
            (*field_values)["Average Sph. Dia. [um]"] =
                    average_diameter.first * 1e6;
            (*field_values)["Sigma Sph. Dia. [um]"] =
                    average_diameter.second * 1e6;
            (*field_values)["Num. Spherulites [-]"] =
                    static_cast<double>(num_spherulites);
            (*field_values)["Spherulite density [(1/m)^3]"] =
                    static_cast<double>(num_spherulites) /
                    (GetCellSizeM() * GetCellSizeM() * GetCellSizeM() * static_cast<double>(GetTotalNumCells()));
        }
    }

    void Element::
    GetVectorFieldsForElement(const double current_point_in_time,
                              const Eigen::Vector3s starting_cell,
                              const Eigen::Vector3s num_cells,
                              std::map<const std::string,
                                      Eigen::Vector3d> *field_values) const {
        std::size_t cells_in_area = 0;
        Eigen::Vector3d summed_velocity(0.0, 0.0, 0.0);
        Eigen::Vector3d location(0.0, 0.0, 0.0);

        for (std::size_t cur_x = 0; cur_x < num_cells[0]; ++cur_x) {
            for (std::size_t cur_y = 0; cur_y < num_cells[1]; ++cur_y) {
                for (std::size_t cur_z = 0; cur_z < num_cells[2]; ++cur_z) {
                    Eigen::Vector3s cur_index = starting_cell +
                                                Eigen::Vector3s(cur_x, cur_y, cur_z);
                    if (cur_index[0] < num_cells_[0] &&
                        cur_index[1] < num_cells_[1] &&
                        cur_index[2] < num_cells_[2]) {

                        // valid index in the element
                        summed_velocity += CellVelocity(AbsoluteCellIndex(cur_index),
                                                        current_point_in_time);
                        location += CellPosition(cur_index, true);
                        ++cells_in_area;
                    }
                }
            }
        }

        if (cells_in_area > 0) {
            const double factor = 1.0 / static_cast<double>(cells_in_area);
            (*field_values)["Location [m]"] = location * factor;
            (*field_values)["Average Vel. [m/s]"] = summed_velocity * factor;
        }
    }

// add a single nucleus at a specific nucleation
    void Element::AddNucleus(const Eigen::Vector3s &location,
                             const double point_in_time) {
        if (GetFullyCrystalline()) {
            return;
        }
        const std::size_t absolute_index = AbsoluteCellIndex(location);
        const uint64_t spherulite_id =
                current_simulation_->NewSpherulite(absolute_index,
                                                   point_in_time,
                                                   this);
        ChangeCellToSolid(absolute_index,
                          spherulite_id,
                          point_in_time,
                          Eigen::Vector3i(0, 0, 0),
                          true);
    }

// add a single nucleus at a random location
    void Element::AddRandomNucleus(const double point_in_time) {
        if (GetFullyCrystalline()) {
            return;
        }

        Eigen::Vector3s location;
        do {
            location = Eigen::Vector3s(rng_->rand() % num_cells_[0],
                                       rng_->rand() % num_cells_[1],
                                       rng_->rand() % num_cells_[2]);
        }
        while (CellCrystalline(AbsoluteCellIndex(location)));
        AddNucleus(location,
                   point_in_time);
    }

    const bool Element::EventExists(StateEvent *event_to_check) const {
        if (!current_simulation_->GetSaveMemory()) {
            return false;
        }
//TODO: Ask what this method does
        auto it = state_events_.upper_bound(event_to_check);
        auto rit = state_events_.rbegin();
        if (it != state_events_.end()) {
            const long distance = std::distance(state_events_.begin(), it);
            for (long dist = 0; dist < distance; ++dist)
                ++rit;
        }

        while (rit != state_events_.rend()) {
            if ((*rit)->SameEvent(event_to_check)) {
                return true;
            }
            if ((*rit)->GetPointInTime() - event_to_check->GetPointInTime() < -1e-6) {
                return false;
            }
            ++rit;
        }
        return false;
    }

    void Element::ChangeCellToSolid(const std::size_t cell_index,
                                    const uint64_t spherulite_id,
                                    const double point_in_time,
                                    const Eigen::Vector3i periodic_continuation,
                                    const bool nucleation) {
        phasestate_[cell_index] = static_cast<double>(spherulite_id);
        crystallization_time_[cell_index] = point_in_time;
        num_crystalline_cells_ += 1;
        current_simulation_->AddCellToSpherulite(spherulite_id,
                                                 complex_index(this, cell_index),
                                                 periodic_continuation);

        // the CCG method is not creating growth events for each cell but only one
        // global one
        if (growth_type_ == GrowthEvent::kCCG) {
            if (nucleation) {
                additional_information_ccg_->crystal_state_[cell_index] = 1e-6;
            } else {
                additional_information_ccg_->periodic_continuation_[cell_index] =
                        Eigen::Vector3i8(static_cast<int8_t>(periodic_continuation[0]),
                                         static_cast<int8_t>(periodic_continuation[1]),
                                         static_cast<int8_t>(periodic_continuation[2]));
            }
            return;
        }

        const bool moore = (growth_type_ == GrowthEvent::kRayTracing);
        std::vector<Eigen::Vector3i> periodic_continuations;
        std::vector<complex_index> indices;
        NeighbourCells(cell_index,
                       moore,
                       true,
                       &indices,
                       &periodic_continuations);

        complex_index spherulite_center =
                current_simulation_->SpheruliteCenter(spherulite_id);

        // add a growth event to the neighbourcells
        #pragma omp parallel
        {
            //double start = omp_get_wtime();
            LIKWID_MARKER_START("ChangeCellToSolid");
            #pragma omp for schedule(dynamic) nowait
            for (size_t cur = 0; cur < indices.size(); ++cur) {
                const Element *cur_el = indices[cur].first;
                if (cur_el != nullptr) {
                    if (!cur_el->CellCrystalline(indices[cur].second)) {
                        // get the periodic continuations that were used for this cell
                        Eigen::Vector3i new_periodic_continuation = periodic_continuation;
                        const Eigen::Vector3i cell_periodic_continuation =
                                periodic_continuations[cur];

                        for (int8_t dim = 0; dim < 3; ++dim) {
                            if (cell_periodic_continuation[dim] != 0 &&
                                cell_periodic_continuation[dim] != periodic_continuation[dim]) {
                                if (periodic_continuation[dim] != 0) {
                                    new_periodic_continuation[dim] = 0;
                                } else {
                                    new_periodic_continuation[dim] = cell_periodic_continuation[dim];
                                }
                            } else {
                                new_periodic_continuation[dim] = periodic_continuation[dim];
                            }
                        }

                        double new_point_in_time;
                        double delta_time = -1.0;
                        if (growth_type_ == GrowthEvent::kMonteCarlo) { // monte carlo method
                            delta_time = current_simulation_->Timestep();
                            new_point_in_time = point_in_time + delta_time;
                        } else if (growth_type_ == GrowthEvent::kRayTracing) {  // ray-tracing method
                            new_point_in_time =
                                    cur_el->RaytracingMethod(indices[cur].second,
                                                             spherulite_center,
                                                             new_periodic_continuation);
                        } else {
                            throw Exception("Internal error",
                                            "Not implemented yet");
                        }
                        GrowthEvent *new_event = nullptr;
                        if (new_point_in_time < 9e11) {
                            new_event =
                                    new GrowthEvent(indices[cur].second,
                                                    complex_index(this, cell_index),
                                                    new_point_in_time,
                                                    delta_time,
                                                    growth_type_,
                                                    new_periodic_continuation,
                                                    current_simulation_->
                                                            GetElement(cur_el->GetElementID()));
                        }
                        #pragma omp critical
                        {
                            if (new_event != nullptr) {
                                if (EventExists(new_event)) {
                                    delete new_event;
                                    new_event = nullptr;
                                } else {
                                    AddStateEvent(new_event);
                                }
                            }
                        }
                    }
                }
            }
            LIKWID_MARKER_STOP("ChangeCellToSolid");
            //printf("T%d: %fs\n", omp_get_thread_num(), omp_get_wtime()-start);
            #pragma omp barrier
        };
        if (indices.size() != 26) {
            printf("#Cells:%zu\n", indices.size());
        }
    }

// return the cell ids of the neighbouring cells
    void Element::NeighbourCells(const std::size_t cell_index,
                                 const bool moore_environment,
                                 const bool check_periodic_continuation,
                                 std::vector<complex_index> *indices,
                                 std::vector<Eigen::Vector3i> *periodic_continuation) const {
        const Eigen::Vector3s cell_position = CellIndices(cell_index);
        const Eigen::Vector3i64 cell_position_i64 = Eigen::Vector3i64(
                static_cast<int64_t>(cell_position[0]),
                static_cast<int64_t>(cell_position[1]),
                static_cast<int64_t>(cell_position[2]));

        // van Neumann neighbourhood
        const Eigen::Vector3i64 offsets_neumann[] = {
                Eigen::Vector3i64(1, 0, 0),
                Eigen::Vector3i64(-1, 0, 0),
                Eigen::Vector3i64(0, 1, 0),
                Eigen::Vector3i64(0, -1, 0),
                Eigen::Vector3i64(0, 0, 1),
                Eigen::Vector3i64(0, 0, -1),
        };
        const std::size_t num_entries_neumann = sizeof(offsets_neumann) /
                                                sizeof(Eigen::Vector3i64);
        for (std::size_t cur = 0; cur < num_entries_neumann; ++cur) {
            const Eigen::Vector3i64 cur_index = cell_position_i64 +
                                                offsets_neumann[cur];
            Eigen::Vector3i periodic_continuation_vector;
            indices->push_back(ComplexIndex(cur_index[0],
                                            cur_index[1],
                                            cur_index[2],
                                            check_periodic_continuation,
                                            &periodic_continuation_vector));
            if (periodic_continuation != nullptr) {
                periodic_continuation->push_back(periodic_continuation_vector);
            }
        }

        if (moore_environment) {
            // additional offset for the moore neighbourhood
            const Eigen::Vector3i64 offsets_moore[] = {
                    // diagonal
                    Eigen::Vector3i64(1, 1, 0),
                    Eigen::Vector3i64(1, -1, 0),
                    Eigen::Vector3i64(-1, 1, 0),
                    Eigen::Vector3i64(-1, -1, 0),

                    // z + 1
                    Eigen::Vector3i64(1, 1, 1),
                    Eigen::Vector3i64(1, -1, 1),
                    Eigen::Vector3i64(-1, 1, 1),
                    Eigen::Vector3i64(-1, -1, 1),
                    Eigen::Vector3i64(1, 0, 1),
                    Eigen::Vector3i64(-1, 0, 1),
                    Eigen::Vector3i64(0, 1, 1),
                    Eigen::Vector3i64(0, -1, 1),

                    // z - 1
                    Eigen::Vector3i64(1, 1, -1),
                    Eigen::Vector3i64(1, -1, -1),
                    Eigen::Vector3i64(-1, 1, -1),
                    Eigen::Vector3i64(-1, -1, -1),
                    Eigen::Vector3i64(1, 0, -1),
                    Eigen::Vector3i64(-1, 0, -1),
                    Eigen::Vector3i64(0, 1, -1),
                    Eigen::Vector3i64(0, -1, -1)
            };
            const std::size_t num_entries_moore = sizeof(offsets_moore) /
                                                  sizeof(Eigen::Vector3i64);
            for (std::size_t cur = 0; cur < num_entries_moore; ++cur) {
                const Eigen::Vector3i64 cur_index = cell_position_i64 +
                                                    offsets_moore[cur];
                Eigen::Vector3i periodic_continuation_vector;
                indices->push_back(ComplexIndex(cur_index[0],
                                                cur_index[1],
                                                cur_index[2],
                                                check_periodic_continuation,
                                                &periodic_continuation_vector));
                if (periodic_continuation != nullptr) {
                    periodic_continuation->push_back(periodic_continuation_vector);
                }
            }
        }
    }

// return the complex index. In contrast to normal indices, complex indices
// allow index values beyond the current element. ComplexIndex checks if there
// are neighbouring elements to which the cell belongs. If not the perdiodic
// continuation flags are checked. If 'check_periodic_continuation' flag is
// false the function returns the closest index in the current element if a
// perdiodic continuation would be applyable
    const complex_index
    Element::ComplexIndex(const int64_t index_x,
                          const int64_t index_y,
                          const int64_t index_z,
                          const bool check_periodic_continuation,
                          Eigen::Vector3i *periodic_continuation) const {
        if (periodic_continuation != nullptr) {
            *periodic_continuation = Eigen::Vector3i(0, 0, 0);
        }

        if (index_x >= 0 && index_x < static_cast<int64_t>(num_cells_[0]) &&
            index_y >= 0 && index_y < static_cast<int64_t>(num_cells_[1]) &&
            index_z >= 0 && index_z < static_cast<int64_t>(num_cells_[2])) {
            return complex_index(this, AbsoluteCellIndex(Eigen::Vector3s(index_x,
                                                                         index_y,
                                                                         index_z)));
        }

        Eigen::Vector3i64 cur_index(index_x, index_y, index_z);
        const Element *owner = this;

        for (std::size_t cur_neighbour = 0; cur_neighbour < 6; ++cur_neighbour) {
            const std::size_t cur_dim = cur_neighbour / 2;
            const bool out_of_bounds = ((cur_neighbour % 2) == 0) ?
                                       cur_index[cur_dim] >= static_cast<int64_t>(owner->num_cells_[cur_dim]) :
                                       cur_index[cur_dim] < 0;
            const int64_t factor = ((cur_neighbour % 2) == 0) ? -1 : 1;

            // neighbour index order: (+/- X_AXIS; +/- Y_AXIS; +/- Z_AXIS)
            if (out_of_bounds) {
                if (!owner->neighbour_infos_[cur_neighbour].second) {  // no slip condition
                    return complex_index(NULL, std::numeric_limits<uint64_t>::max());
                }

                if (owner->neighbour_infos_[cur_neighbour].first == NULL) {
                    // a neighbouring element should exist, but no data exists ->
                    //   check periodic continuation
                    if (check_periodic_continuation) {
                        Eigen::Vector3i direction(0, 0, 0);
                        direction[cur_dim] = static_cast<int>(factor);
                        if (periodic_continuation != nullptr) {
                            (*periodic_continuation)[cur_dim] = direction[cur_dim];
                        }

                        const Element *old_owner = owner;
                        owner = owner->PeriodicBoundaryNeighbour(direction);
                        if (owner == nullptr) {  // no slip condition
                            return complex_index(NULL, std::numeric_limits<uint64_t>::max());
                        }

                        if (factor == -1) { // >= num_cells
                            cur_index[cur_dim] = cur_index[cur_dim] -
                                                 old_owner->num_cells_[cur_dim];
                            if (cur_index[cur_dim] >=
                                static_cast<int64_t>(owner->num_cells_[cur_dim])) {
                                cur_index[cur_dim] = owner->num_cells_[cur_dim] - 1;
                            }
                        } else { // < 0
                            cur_index[cur_dim] += owner->num_cells_[cur_dim];
                            if (cur_index[cur_dim] < 0) {
                                cur_index[cur_dim] = 0;
                            }
                        }
                    } else {
                        cur_index[cur_dim] = (factor == -1) ?
                                             owner->num_cells_[cur_dim] - 1 :
                                             0;  // nearest cell
                    }
                } else {
                    // the cell belongs to the neighbouring element
                    const Element *old_owner = owner;
                    owner = owner->neighbour_infos_[cur_neighbour].first;
                    if (factor == -1) { // >= num_cells
                        cur_index[cur_dim] = cur_index[cur_dim] - old_owner->num_cells_[cur_dim];
                        if (cur_index[cur_dim] >=
                            static_cast<int64_t>(owner->num_cells_[cur_dim])) {
                            cur_index[cur_dim] = owner->num_cells_[cur_dim] - 1;
                        }
                    } else { // < 0
                        cur_index[cur_dim] += owner->num_cells_[cur_dim];
                        if (cur_index[cur_dim] < 0) {
                            cur_index[cur_dim] = 0;
                        }
                    }
                }
            }
        }
        return complex_index(owner,
                             owner->AbsoluteCellIndex(Eigen::Vector3s(cur_index[0],
                                                                      cur_index[1],
                                                                      cur_index[2])));
    }

    const Element *Element::
    PeriodicBoundaryNeighbour(const Eigen::Vector3i direction) const {
        // neighbour index order: (+/- X_AXIS; +/- Y_AXIS; +/- Z_AXIS)
        std::size_t neighbour_index;
        if (direction[0] == 1) {
            neighbour_index = 0;
        } else if (direction[0] == -1) {
            neighbour_index = 1;
        } else if (direction[1] == 1) {
            neighbour_index = 2;
        } else if (direction[1] == -1) {
            neighbour_index = 3;
        } else if (direction[2] == 1) {
            neighbour_index = 4;
        } else if (direction[2] == -1) {
            neighbour_index = 5;
        } else {
            throw Exception("Internal error",
                            "Wrong neighbourhood direction");
        }

        const Element *cur_element = this;
        while (cur_element->neighbour_infos_[neighbour_index].first != NULL) {
            // check no slip condition
            if (!cur_element->neighbour_infos_[neighbour_index].second) {
                return nullptr;
            }
            cur_element = cur_element->neighbour_infos_[neighbour_index].first;
        }
        return cur_element;
    }

// execute a monte-carlo step for a cell
    bool Element::MonteCarloStep(const double probability) const {
        const double comparison = rng_->rand() /
                                  static_cast<double>(std::numeric_limits<unsigned int>::max());
        return comparison < probability;
    }

// retrieve the center position of the element
    const Eigen::Vector3d Element::CenterPosition() const {
        return lower_left_front_corner_m_ + 0.5 * size_m_;
    }

// calculate the crystal growth speed at a given point in time
    const double Element::CrystalGrowthSpeed(const double point_in_time,
                                             const std::size_t cell_index) const {
        return current_simulation_->GetGrowthModel()->Speed(point_in_time,
                                                            cell_index,
                                                            this);
    }

// calculate the nucleation rate at a given point in time
    const double Element::NucleationRate(const double point_in_time,
                                         const std::size_t cell_index) const {
        return current_simulation_->GetNucleationModel()->Rate(point_in_time,
                                                               cell_index,
                                                               !current_simulation_->GetCalculateFlowEnergy(),
                                                               false,
                                                               this);
    }

// integrate nucleation rate or growth speed over time
    const double Element::IntegrateCellOverTime(const double start_time,
                                                const double end_time,
                                                const std::size_t cell_index,
                                                const bool nucleation) const {
        class Integration {
        public:
            Integration(const std::size_t cell_index,
                        const double scaling,
                        const bool nucleation,
                        const Element *element) :
                    cell_index_(cell_index),
                    scaling_(scaling),
                    nucleation_(nucleation),
                    el_(element) {
            }

            Integration(const Integration &src) = default;

            double operator()(const double point_in_time) const {
                if (nucleation_) {
                    return scaling_ * el_->NucleationRate(point_in_time,
                                                          cell_index_);
                } else {
                    return scaling_ * el_->CrystalGrowthSpeed(point_in_time,
                                                              cell_index_);
                }
            }

        private:
            const double scaling_;
            const std::size_t cell_index_;
            const bool nucleation_;
            const Element *el_;

            Integration &operator=(const Integration &) {
                return *this;
            }
        };

        const double scaling =
                (!nucleation) ?
                current_simulation_->GetGrowthModel()->ScalingFactor() :
                current_simulation_->GetNucleationModel()->ScalingFactor();

        const double integral =
                DEIntegrator<Integration>::Integrate(Integration(cell_index,
                                                                 scaling,
                                                                 nucleation,
                                                                 this),
                                                     start_time,
                                                     end_time,
                                                     current_simulation_->GetIntegralError());

        return integral / scaling;
    }

// return the time delta needed to reach a certain value
    const double Element::
    IntegrateCellToValue(const std::size_t cell_index,
                         const double dest_value,
                         const double start_time,
                         const bool nucleation) const {
        if (dest_value < 0.0) {
            throw Exception("Internal error",
                            "Destination value must be larger than zero.");
        }

        const double integration_step = current_simulation_->Timestep();

        double current_value = 0.0;
        double current_delta = 0.0;
        double current_time = start_time;
        int64_t no_improvement_counter = 0;
        while (current_value < dest_value) {
            current_delta = IntegrateCellOverTime(current_time,
                                                  current_time + integration_step,
                                                  cell_index,
                                                  nucleation);
            if (fabs(current_delta) < current_value * 1e-6) {
                ++no_improvement_counter;
            } else {
                current_value += current_delta;
                no_improvement_counter = 0;
            }
            if (no_improvement_counter == 10000) {
                return 1e20;
            }  // will take very long
            current_time += integration_step;
        }

        // subtract the overshoot
        const double overshoot = current_value - dest_value;
        const double percentage = overshoot / current_delta;
        current_time -= integration_step * percentage;
        return current_time - start_time;
    }

// execute a periodic continuation on a position
    const Eigen::Vector3d
    Element::PeriodicCorrection(const Eigen::Vector3i &periodic_continuation,
                                const Element *min_neighbours[3],
                                const Element *max_neighbours[3],
                                const Eigen::Vector3d &position,
                                const Eigen::Vector3d &end_position) const {
        Eigen::Vector3d new_position = position;
        for (int32_t cur = 0; cur < 3; ++cur) {
            if (periodic_continuation[cur] != 0) {
                // neighbour index order: (+/- X_AXIS; +/- Y_AXIS; +/- Z_AXIS)
                if (periodic_continuation[cur] == 1) { // continuation to positive direction
                    if (position[cur] < end_position[cur] - 1e-12) {
                        new_position[cur] = max_neighbours[cur]->GetCornerPosition(7)[cur] +
                                            position[cur] - min_neighbours[cur]->GetCornerPosition(0)[cur];
                    }
                } else {
                    if (position[cur] > end_position[cur] + 1e-12) {
                        new_position[cur] = min_neighbours[cur]->GetCornerPosition(0)[cur] -
                                            (max_neighbours[cur]->GetCornerPosition(7)[cur] - position[cur]);
                    }
                }
            }
        }
        return new_position;
    }

// get the minimum/maximum neighbour of the simulation area
    void Element::MinMaxNeighbour(const Element *min_neighbours[3],
                                  const Element *max_neighbours[3]) const {
        for (size_t cur = 0; cur < 3; ++cur) {
            min_neighbours[cur] = this;
            while (min_neighbours[cur]->GetNeighbourElement(1 + cur * 2) != nullptr) {
                min_neighbours[cur] = min_neighbours[cur]->GetNeighbourElement(1 + cur * 2);
            }

            max_neighbours[cur] = this;
            while (max_neighbours[cur]->GetNeighbourElement(cur * 2) != nullptr) {
                max_neighbours[cur] = max_neighbours[cur]->GetNeighbourElement(cur * 2);
            }
        }
    }

// calculate the point in time when the cell at absolute_index is fully covered
// by the growth front of the spherulite at spherulite_center
    const double
    Element::RaytracingMethod(const std::size_t absolute_index,
                              const complex_index spherulite_center,
                              const Eigen::Vector3i periodic_continuation) const {
        Eigen::Vector3d start_pos =
                spherulite_center.first->CellPosition(spherulite_center.second, true);
        Eigen::Vector3d end_center_pos = CellPosition(absolute_index, true);

        const Element *min_neighbours[3];
        const Element *max_neighbours[3];
        MinMaxNeighbour(min_neighbours, max_neighbours);

        // adjust the startposition due to periodic continuation
        start_pos = PeriodicCorrection(periodic_continuation,
                                       min_neighbours,
                                       max_neighbours,
                                       start_pos,
                                       end_center_pos);

        // direction of the ray
        Eigen::Vector3d direction = end_center_pos - start_pos;
        direction.normalize();

        // increment of the cell index
        const Eigen::Vector3i64 increment(direction[0] < -1e-20 ? -1 : 1,
                                          direction[1] < -1e-20 ? -1 : 1,
                                          direction[2] < -1e-20 ? -1 : 1);

        // Determine how far we can travel along the ray before we hit a voxel boundary.
        Eigen::Vector3d max_dist(0.5 * cell_size_m_ * static_cast<double>(increment[0]),
                                 0.5 * cell_size_m_ * static_cast<double>(increment[1]),
                                 0.5 * cell_size_m_ * static_cast<double>(increment[2]));
        Eigen::Vector3d delta(cell_size_m_ * static_cast<double>(increment[0]),
                              cell_size_m_ * static_cast<double>(increment[1]),
                              cell_size_m_ * static_cast<double>(increment[2]));
        for (std::size_t cur = 0; cur < 3; ++cur) {
            if (fabs(direction[cur]) < 1e-20) {
                max_dist[cur] = 1e20;
                delta[cur] = 1e20;
            } else {
                max_dist[cur] /= direction[cur];
                delta[cur] /= direction[cur];
            }
        }

        const double start_time =
                spherulite_center.first->crystallization_time_[spherulite_center.second];

        double current_time = start_time;
        complex_index end_position(this, absolute_index);
        complex_index cur_position = spherulite_center;
        complex_index next_position = spherulite_center;
        do {
            cur_position = next_position;

            Eigen::Vector3d cell_position =
                    cur_position.first->CellPosition(cur_position.second, true);

            // constraints
            if (cur_position != end_position) {
                //if (cur_position.first->CellCrystalline(cur_position.second) == false) {
                //  return 9e11;  // not allowed
                //}
                // if (cur_position.first->GetSpheruliteID(cur_position.second) !=
                //     spherulite_center.first->GetSpheruliteID(spherulite_center.second))
                //   return 9e11;  // not allowed
            }

            // move because of the periodic continuation,
            cell_position = PeriodicCorrection(periodic_continuation,
                                               min_neighbours,
                                               max_neighbours,
                                               cell_position,
                                               end_center_pos);

            // min position
            cell_position[0] -= cell_size_m_ * 0.5;
            cell_position[1] -= cell_size_m_ * 0.5;
            cell_position[2] -= cell_size_m_ * 0.5;

            double delta_distance = DeltaDistance(cell_position,
                                                  start_pos,
                                                  direction);
            if (delta_distance == 0.0) {
                printf("ERROR: delta_distance == 0.0 (Element.cpp, Line 2241)\n");
            }

            if (cur_position == spherulite_center) {
                delta_distance *= 0.5;
            }  // only half the distance in the spherulite
            // center cell

            if (fabs(delta_distance) > 1e-9) {
                const double time_delta =
                        cur_position.first->IntegrateCellToValue(cur_position.second,
                                                                 delta_distance,
                                                                 current_time,
                                                                 false);
                current_time += time_delta;
                if (current_time > 9e11) {
                    printf("INFO: current_time > 9e11 (Element.cpp, Line 2259)\n");
                    return 9e11;  // will take a long time...
                }
            }

            Eigen::Vector3s cur_index_unsigned =
                    cur_position.first->CellIndices(cur_position.second);
            Eigen::Vector3i64 cur_index(cur_index_unsigned[0],
                                        cur_index_unsigned[1],
                                        cur_index_unsigned[2]);
            // Do the next step.
            if (fabs(max_dist[0] - max_dist[1]) < 1e-12 &&
                fabs(max_dist[0] - max_dist[2]) < 1e-12) {
                cur_index += increment;
                max_dist += delta;
            } else if (fabs(max_dist[0] - max_dist[1]) < 1e-12 &&
                       max_dist[0] < max_dist[2]) {
                cur_index[0] += increment[0];
                cur_index[1] += increment[1];
                max_dist[0] += delta[0];
                max_dist[1] += delta[1];
            } else if (fabs(max_dist[0] - max_dist[2]) < 1e-12 &&
                       max_dist[0] < max_dist[1]) {
                cur_index[0] += increment[0];
                cur_index[2] += increment[2];
                max_dist[0] += delta[0];
                max_dist[2] += delta[2];
            } else if (fabs(max_dist[1] - max_dist[2]) < 1e-12 &&
                       max_dist[1] < max_dist[0]) {
                cur_index[1] += increment[1];
                cur_index[2] += increment[2];
                max_dist[1] += delta[1];
                max_dist[2] += delta[2];
            } else if (fabs(max_dist[0]) < fabs(max_dist[1]) &&
                       max_dist[0] < max_dist[2]) {
                // tMax.X is the lowest, an YZ cell boundary plane is nearest.
                cur_index[0] += increment[0];
                max_dist[0] += delta[0];
            } else if (max_dist[1] < max_dist[2]) {
                // tMax.Y is the lowest, an XZ cell boundary plane is nearest.
                cur_index[1] += increment[1];
                max_dist[1] += delta[1];
            } else {
                // tMax.Z is the lowest, an XY cell boundary plane is nearest.
                cur_index[2] += increment[2];
                max_dist[2] += delta[2];
            }
            next_position = cur_position.first->ComplexIndex(cur_index[0],
                                                             cur_index[1],
                                                             cur_index[2],
                                                             true);
        }
        while (cur_position != end_position);
        if (current_time - start_time < 1e-12) {
            current_time = start_time;
        }

        return current_time;
    }

// calculate the intersection with the cell walls and the distance between it
    const double Element::DeltaDistance(const Eigen::Vector3d &min_cell_position,
                                        const Eigen::Vector3d &spherulite_center,
                                        const Eigen::Vector3d &direction) const {

        Eigen::Vector3d min_pos, max_pos;
        std::vector<Eigen::Vector3d> results;
        Plane planes[] = {Plane(Eigen::Vector3d(-1.0, 0.0, 0.0)),
                          Plane(Eigen::Vector3d(0.0, -1.0, 0.0)),
                          Plane(Eigen::Vector3d(0.0, 0.0, -1.0))};


        min_pos = min_cell_position;
        max_pos = min_pos;
        max_pos[0] += cell_size_m_;
        max_pos[1] += cell_size_m_;
        max_pos[2] += cell_size_m_;

        // a straight line g is given by: x = v + a*u
        // v is position vector, it is given by v = spherulite_center
        // u is the directional vector, it is given by u = direction
        // check the six planes of the cell box
        // check min
        for (int32_t cur = 0; cur < 3; ++cur) {
            // calculate intersections
            planes[cur].set_distance_from_point(min_pos);
            StoreIntersectionResults(planes[cur],
                                     direction,
                                     spherulite_center,
                                     min_pos,
                                     max_pos,
                                     &results);
            if (results.size() == 2) {
                break;
            }
            planes[cur].set_distance_from_point(max_pos);
            StoreIntersectionResults(planes[cur],
                                     direction,
                                     spherulite_center,
                                     min_pos,
                                     max_pos,
                                     &results);
            if (results.size() == 2) {
                break;
            }
        }
        // check if there are two intersection point
        if (results.size() < 2) {
            return 0.0;
        }

        // check if all point are in the cell
        if (results[0][0] < min_pos[0] - 1e-12 ||
            results[0][1] < min_pos[1] - 1e-12 ||
            results[0][2] < min_pos[2] - 1e-12 ||
            results[0][0] > max_pos[0] + 1e-12 ||
            results[0][1] > max_pos[1] + 1e-12 ||
            results[0][2] > max_pos[2] + 1e-12) {
            return 0.0;
        }
        if (results[1][0] < min_pos[0] - 1e-12 ||
            results[1][1] < min_pos[1] - 1e-12 ||
            results[1][2] < min_pos[2] - 1e-12 ||
            results[1][0] > max_pos[0] + 1e-12 ||
            results[1][1] > max_pos[1] + 1e-12 ||
            results[1][2] > max_pos[2] + 1e-12) {
            return 0.0;
        }

        Eigen::Vector3d distance = results[0] - results[1];
        return distance.norm();
    }

    void Element::
    StoreIntersectionResults(const Plane &plane,
                             const Eigen::Vector3d &direction,
                             const Eigen::Vector3d &center,
                             const Eigen::Vector3d &min_pos,
                             const Eigen::Vector3d &max_pos,
                             std::vector<Eigen::Vector3d> *results) const {
        bool parallel = false;
        double factor = plane.Intersection(center,
                                           direction,
                                           &parallel);
        if (parallel) {
            return;
        }

        Eigen::Vector3d possible_intersection = center + direction * factor;
        bool store = possible_intersection[0] <= max_pos[0] + 1e-12 &&
                     possible_intersection[1] <= max_pos[1] + 1e-12 &&
                     possible_intersection[2] <= max_pos[2] + 1e-12 &&
                     possible_intersection[0] + 1e-12 >= min_pos[0] &&
                     possible_intersection[1] + 1e-12 >= min_pos[1] &&
                     possible_intersection[2] + 1e-12 >= min_pos[2];
        if (store) {
            // check if the point is already in the list
            bool already_in_list = false;
            for (std::size_t cur = 0; cur < results->size(); ++cur) {
                Eigen::Vector3d res = results->at(cur) - possible_intersection;
                double dist = res.dot(res);
                if (dist < 1e-30) {
                    already_in_list = true;
                    break;
                }
            }
            if (!already_in_list) {
                results->push_back(possible_intersection);
            }
        }
    }

    void Element::WriteFinalResults(const double last_executed_event) {
        // execute an additional 'WriteResults' event for the element
        WriteResultEvent result_event(last_executed_event,
                                      -1.0,  // no automatic generation of new events
                                      *current_simulation_,
                                      this);
        result_event.Execute();
        has_written_final_results_ = true;
    }

    const double Element::
    CellVelocityGradient(const std::size_t cell_index,
                         const std::size_t component,
                         const std::size_t axis,
                         const double current_time) const {
        const Eigen::Vector3s indices = CellIndices(cell_index);
        Eigen::Vector3s prev = indices;

        double delta = 0.0;
        if (prev[axis] != 0) {
            prev[axis] -= 1;
            delta += GetCellSizeM();
        }
        Eigen::Vector3s next = indices;
        if (next[axis] != GetNumCells()[axis] - 1) {
            next[axis] += 1;
            delta += GetCellSizeM();
        }

        if (std::abs(delta) < 1e-9) {
            return 0.0;
        }

        const std::size_t abs_prev = AbsoluteCellIndex(prev);
        const std::size_t abs_next = AbsoluteCellIndex(next);

        const double delta_component =
                CellVelocity(abs_next, current_time)[component] -
                CellVelocity(abs_prev, current_time)[component];
        return delta_component / delta;
    }

    bool Element::CellInBorder(const size_t index) const {
        std::vector<complex_index> neighbour_indices;
        NeighbourCells(index, false, true, &neighbour_indices);

        // never in border if not crystalline
        if (!CellCrystalline(index)) {
            return false;
        }

        auto sph_id = GetSpheruliteID(index);

        // always in border if at fixed side of the simulation area
        for (auto &neighbour_indice : neighbour_indices) {
            if (neighbour_indice.first == NULL) {
                return true;
            }
        }

        // check positive direction neighbours if they belong to a different
        // spherulite
        int indices_to_check[] = {0, 2, 4};
        for (int &cur_id : indices_to_check) {
            if (neighbour_indices[cur_id].first->
                    CellCrystalline(neighbour_indices[cur_id].second)) {
                auto neigh_id =
                        neighbour_indices[cur_id].first->GetSpheruliteID(neighbour_indices[cur_id].second);
                if (neigh_id != sph_id) {
                    return true;
                }
            }
        }

        return false;
    }

// hashkey for a cell at a given point in time
    const uint64_t Element::CellHashkey(const std::size_t absolute_index,
                                        const double point_in_time,
                                        std::vector<double> *history_temperature_times,
                                        std::vector<double> *history_velocity_times) const {
        uint64_t hash_value = 0;
        const double cache_convert_temperature =
                GetCurrentSimulation()->GetCacheConvertTemperature();
        const double cache_convert_velocity =
                GetCurrentSimulation()->GetCacheConvertVelocity();

        // temperatures
        if (history_temperature_times->empty()) {
            temperatures_->PointInTimeOfEntries(point_in_time,
                                                history_temperature_times);
        }
        for (auto &history_temperature_time : *history_temperature_times) {
            const double temperature =
                    CellTemperature(absolute_index, history_temperature_time);
            const auto temperature_data =
                    static_cast<uint64_t>(temperature * cache_convert_temperature + 0.5);

            hash_value = HashValue::I()->Recalculate(hash_value, temperature_data);
        }
        const double temperature = CellTemperature(absolute_index, point_in_time);
        const auto temperature_data =
                static_cast<uint64_t>(temperature * cache_convert_temperature + 0.5);
        hash_value = HashValue::I()->Recalculate(hash_value, temperature_data);

        if (GetCurrentSimulation()->GetCalculateFlowEnergy()) {
            // velocities
            if (history_velocity_times->empty()) {
                velocities_->PointInTimeOfEntries(point_in_time, history_velocity_times);
            }
            for (auto &history_velocity_time : *history_velocity_times) {
                const Eigen::Vector3d vel =
                        CellVelocity(absolute_index, history_velocity_time);
                const uint64_t vel_data_x =
                        static_cast<uint64_t>(vel[0] * cache_convert_velocity + 0.5);
                const uint64_t vel_data_y =
                        static_cast<uint64_t>(vel[1] * cache_convert_velocity + 0.5);
                const uint64_t vel_data_z =
                        static_cast<uint64_t>(vel[2] * cache_convert_velocity + 0.5);
                hash_value = HashValue::I()->Recalculate(hash_value, vel_data_x);
                hash_value = HashValue::I()->Recalculate(hash_value, vel_data_y);
                hash_value = HashValue::I()->Recalculate(hash_value, vel_data_z);
            }
            const Eigen::Vector3d vel = CellVelocity(absolute_index, point_in_time);
            const uint64_t vel_data_x =
                    static_cast<uint64_t>(vel[0] * cache_convert_velocity + 0.5);
            const uint64_t vel_data_y =
                    static_cast<uint64_t>(vel[1] * cache_convert_velocity + 0.5);
            const uint64_t vel_data_z =
                    static_cast<uint64_t>(vel[2] * cache_convert_velocity + 0.5);
            hash_value = HashValue::I()->Recalculate(hash_value, vel_data_x);
            hash_value = HashValue::I()->Recalculate(hash_value, vel_data_y);
            hash_value = HashValue::I()->Recalculate(hash_value, vel_data_z);
        }
        return hash_value;
    }

}  // namespace SphaeroSim