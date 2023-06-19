#ifndef _SPHAEROSIM_ELEMENT_H_
#define _SPHAEROSIM_ELEMENT_H_

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

#ifndef _INCLUDED_SET_H_

#include <set>

#define _INCLUDED_SET_H_
#endif

#ifndef _INCLUDED_UTILITY_H_

#include <utility>

#define _INCLUDED_UTILITY_H_
#endif

#ifndef _SPHAEROSIM_EIGEN_LIBRARY_H_

#include "EigenLibrary.h"

#endif

#ifndef _SPHAEROSIM_STATE_EVENT_H_

#include "StateEvent.h"

#endif

#ifndef _SPHAEROSIM_INTERPOLATION_H_

#include "Interpolation.h"

#endif

#ifndef _SPHAEROSIM_CACHE_H_

#include "Cache.h"

#endif

#ifndef _SPHAEROSIM_BOUNDARY_CONDITION_H_

#include "BoundaryCondition.h"

#endif

namespace SphaeroSim {

    template<class T>
    class BoundaryCondition;

    class RandomNumberGenerator;

    class Simulation;

// the '_m' at the end of a variable denotes the SI unit meter
    class Element {
    public:
        enum FastAccessTables {
            kSin,
            kCos,
            kLog
        };

        Element(const uint64_t element_id,
                const Eigen::Vector3d minimum_position_m,
                const double cell_size_m,
                const Eigen::Vector3s num_cells,
                const InterpolationType temporal_interpolation_type,
                const InterpolationType spatial_interpolation_type,
                const GrowthEvent::Type growth_type,
                RandomNumberGenerator *rng,
                Simulation *current_simulation);

        Element(const uint64_t element_id,
                const Eigen::Vector3d minimum_position_m,
                const double cell_size_m,
                const Eigen::Vector3d size_m,
                const InterpolationType temporal_interpolation_type,
                const InterpolationType spatial_interpolation_type,
                const GrowthEvent::Type growth_type,
                RandomNumberGenerator *rng,
                Simulation *current_simulation);

        virtual ~Element();

        // fast access tables
        static const double TableValue(const FastAccessTables table,
                                       const double value);

        // precalculate the boundary conditions
        void PrecalculateBoundaryFields();

        // read operations
        const Eigen::Vector3d GetCornerPosition(const std::size_t index) const;

        const Eigen::Vector3d CenterPosition() const;

        const Eigen::Vector3s
        CellIndicesFromPosition(const Eigen::Vector3d &position,
                                const bool force_inside) const;

        const bool CoordinateInElement(const Eigen::Vector3d &position) const;

        const bool CoordinateInElement(const Eigen::Vector3s &indices) const;

        const Eigen::Vector3d PeriodicCorrection(const Eigen::Vector3i &periodic_continuation,
                                                 const Element *min_neighbours[3],
                                                 const Element *max_neighbours[3],
                                                 const Eigen::Vector3d &position,
                                                 const Eigen::Vector3d &end_position) const;

        void NeighbourCells(const std::size_t cell_index,
                            const bool moore_environment,
                            const bool check_periodic_continuation,
                            std::vector<complex_index> *indices,
                            std::vector<Eigen::Vector3i> *periodic_continuation = NULL) const;

        const Eigen::Vector3d CellPosition(const std::size_t cell_index,
                                           const bool cell_center) const;

        const Eigen::Vector3d CellPosition(const Eigen::Vector3s &index_triple,
                                           const bool cell_center) const;

        const double RelativeCrystallizationDegree() const;

        void MinMaxNeighbour(const Element *min_neighbours[3],
                             const Element *max_neighbours[3]) const;

        const uint64_t CellHashkey(const std::size_t absolute_index,
                                   const double point_in_time,
                                   std::vector<double> *history_temperature_times,
                                   std::vector<double> *history_velocity_times) const;

        // add a single nucleus at a specific nucleation
        void AddNucleus(const Eigen::Vector3s &location,
                        const double point_in_time);

        // add a single nucleus at a random location
        void AddRandomNucleus(const double point_in_time);

        // get a complete set of information for state of a cell
        void GetNumericalFieldsForCell(const std::size_t &index,
                                       const double &current_point_in_time,
                                       std::map<const std::string,
                                               double> *field_values) const;

        void GetVectorFieldsForCell(const std::size_t index,
                                    const double current_point_in_time,
                                    std::map<const std::string,
                                            Eigen::Vector3d> *field_values) const;

        void GetNumericalFieldsForCellOLD(const std::size_t &index,
                                          const double &current_point_in_time,
                                          std::map<const std::string,
                                                  double> *field_values) const;

        void GetVectorFieldsForCellOLD(const std::size_t index,
                                       const double current_point_in_time,
                                       std::map<const std::string,
                                               Eigen::Vector3d> *field_values) const;

        void GetNumericalFieldsForElement(const double current_point_in_time,
                                          const Eigen::Vector3s starting_cell,
                                          const Eigen::Vector3s num_cells,
                                          std::map<const std::string,
                                                  double> *field_values) const;

        void GetVectorFieldsForElement(const double current_point_in_time,
                                       const Eigen::Vector3s starting_cell,
                                       const Eigen::Vector3s num_cells,
                                       std::map<const std::string,
                                               Eigen::Vector3d> *field_values) const;

        const double CellTemperature(const std::size_t cell_index,
                                     const double current_time,
                                     const bool check_precalculated = true) const;

        const double CellCoolingRate(const std::size_t cell_index,
                                     const double current_time) const;

        const double CellFlowEnergy(const std::size_t cell_index,
                                    const double current_time,
                                    double *enthropy) const;

        const double CellFlowIntegral(const std::size_t cell_index,
                                      const double current_time) const;

        const Eigen::Vector3d CellVelocity(const std::size_t cell_index,
                                           const double current_time,
                                           const bool check_precalculated = true) const;

        const double CellVelocityGradient(const std::size_t cell_index,
                                          const std::size_t component,
                                          const std::size_t axis,
                                          const double current_time) const;

        // neighbouring element indices
        // first entry: element id of the neighbouring element
        //              a value of UINT64_MAX denotes
        //              no slip conditions
        // second entry: true if boundary data exist for the neighbour-element
        //               false otherwise
        typedef std::pair<uint64_t, bool> neighbour_info;

        // index order: (+/- X_AXIS; +/- Y_AXIS; +/- Z_AXIS)
        void
        SetNeighbourElementInfos(const std::vector<neighbour_info> &infos,
                                 const std::map<const uint64_t, Element *> &elements);

        // point for the boundary conditions
        // first entry: point in time
        // second entry: 8 corner values
        typedef std::pair<const double, const std::vector<double> > temperature_point;
        typedef std::pair<const double,
                const std::vector<Eigen::Vector3d> > velocity_point;

        void SetTemperatureField(const std::vector<temperature_point> &field);

        void SetVelocityField(const std::vector<velocity_point> &field);

        void WriteFinalResults(const double last_executed_event);

        // remove all state events whose point in time is before 'time'
        inline void RemoveStateEvents(const double time) {
            NullEvent tmp_event(time, this);
            std::multiset<StateEvent *>::iterator remove =
                    state_events_.lower_bound(&tmp_event);
            std::multiset<StateEvent *>::iterator it = state_events_.begin();
            while (it != remove) {
                delete (*it);
                ++it;
            }
            state_events_.erase(state_events_.begin(), remove);
        }

        // remove a stateevent
        inline std::multiset<StateEvent*>::iterator RemoveStateEvent(const StateEvent *event_to_remove) {
            std::multiset<StateEvent *>::iterator it = state_events_.begin();
            while ((*it) != event_to_remove) {
                ++it;
            }
            delete (*it);
            return state_events_.erase(it);
        }

        inline std::multiset<StateEvent*>::iterator RemoveStateEvent(const std::multiset<StateEvent*>::iterator &event_to_remove) {
            delete (*event_to_remove);
            return state_events_.erase(event_to_remove);
        }

        inline void AddStateEvent(StateEvent *new_event) {
            state_events_.insert(new_event);
        }

        // Getter
        inline const std::size_t GetTotalNumCells() const {
            return num_cells_[0] * num_cells_[1] * num_cells_[2];
        }

        inline const Eigen::Vector3d GetSizeM() const {
            return size_m_;
        }

        inline const uint64_t GetElementID() const {
            return element_id_;
        }

        inline const bool GetNeighbourDataExist(const std::size_t index) const {
            return neighbour_infos_[index].second &&
                   neighbour_infos_[index].first != NULL;
        }

        inline std::multiset<StateEvent *>::iterator GetStateEventIterator() {
            return state_events_.begin();
        }

        inline std::multiset<StateEvent *>::iterator GetStateEventEnd() const {
            return state_events_.end();
        }

        inline std::multiset<StateEvent *>::iterator
        GetStateEventLowerBound(const double time) {
            NullEvent tmp(time, this);
            return state_events_.lower_bound(&tmp);
        }

        inline const bool GetFullyCrystalline() const {
            if (num_crystalline_cells_ > GetTotalNumCells()) {
                throw Exception("Internal error",
                                "Too many crystalline cells.");
            }
            return num_crystalline_cells_ == GetTotalNumCells();
        }

        inline const Eigen::Vector3s GetNumCells() const {
            return num_cells_;
        }

        inline const double GetVolume() const {
            return size_m_[0] * size_m_[1] * size_m_[2];
        }

        inline const uint64_t GetSpheruliteID(const std::size_t index) const {
            if (!CellCrystalline(index)) {
                throw Exception("Internal error",
                                "Cell not crystalline.");
            }
            return static_cast<uint64_t>(phasestate_[index] + 0.5);
        }

        inline const Element *GetNeighbourElement(const std::size_t index) const {
            return neighbour_infos_[index].first;
        }

        inline bool GetHasWrittenFinalResults() const {
            return has_written_final_results_;
        }

        inline const double GetAdditionalCCGState(const std::size_t index) const {
            return additional_information_ccg_->crystal_state_[index];
        }

        inline const Eigen::Vector3i8
        GetAdditionalCCGPeriodicContinuation(const std::size_t index) const {
            return additional_information_ccg_->periodic_continuation_[index];
        }

        inline void IncreaseAdditionalCCGState(const std::size_t index,
                                               const double increment) const {
            additional_information_ccg_->crystal_state_[index] += increment;
        }

        inline const double GetCellSizeM() const {
            return cell_size_m_;
        }

        inline const GrowthEvent::Type GetCrystalGrowthType() const {
            return growth_type_;
        }

        inline const double GetPhasestate(const std::size_t index) const {
            return phasestate_[index];
        }

        inline const Simulation *GetCurrentSimulation() const {
            return current_simulation_;
        }

        inline const Eigen::Vector3s
        CellIndices(const std::size_t absolute_index) const {
            return precalculated_index_vectors_[absolute_index];
        }

        // convert the three dimensionale index representation to an absolut index
        inline const std::size_t
        AbsoluteCellIndex(const Eigen::Vector3s &index_triple) const {
            return index_triple[2] * (num_cells_[0] * num_cells_[1]) +
                   index_triple[1] * num_cells_[0] +
                   index_triple[0];
        }

        friend GrowthEvent;
        friend NucleationEvent;
    private:
        class Plane {
        public:
            explicit Plane(const Eigen::Vector3d &normal);

            inline void set_distance_from_point(const Eigen::Vector3d &point) {
                distance_ = normal_.dot(point);
            }

            // returns the a-value for the intersection with a
            // straight line g: x = v + a*u
            const double Intersection(const Eigen::Vector3d &pos,
                                      const Eigen::Vector3d &direction,
                                      bool *parallel) const;

        private:
            double distance_;
            Eigen::Vector3d normal_;
        };

        // additional fields that are only needed for the CCG growth method
        class AdditionalInformationCCG {
        public:
            explicit AdditionalInformationCCG(const int64_t cell_count);

            std::vector<double> crystal_state_;
            std::vector<Eigen::Vector3i8> periodic_continuation_;
        };

        // precalculated boundary condition
        template<class T>
        class PrecalculatedBoundaryCondition {
        public:
            PrecalculatedBoundaryCondition(const std::size_t total_num_cells,
                                           const std::size_t total_num_values_per_cell,
                                           const double time_resolution,
                                           const T neutral_element) :
                    total_num_cells_(total_num_cells),
                    time_resolution_(time_resolution),
                    neutral_element_(neutral_element) {

                data_indices=(bool*) malloc(sizeof(bool)*total_num_cells);
                for (int j = 0; j < total_num_cells; ++j) {
                    data_indices[j]= false;
                }

                cell_values_.resize(total_num_cells_);
                for (std::size_t cur = 0; cur < total_num_cells_; ++cur)
                    cell_values_[cur].reserve(total_num_values_per_cell);
                point_in_times_.resize(total_num_values_per_cell);
                lastPrecalculatedValue=-1;
            }

            void AssignNeutralValue(const std::size_t value_index,
                                    const double point_in_time) {
                point_in_times_[value_index] = point_in_time;
                data_indices[value_index]= false;
            }

            void AssignPrecalculatedValue(const std::size_t cell_index,
                                          const std::size_t value_index,
                                          const double point_in_time,
                                          const T &value) {
                point_in_times_[value_index] = point_in_time;
                lastPrecalculatedValue=std::max(lastPrecalculatedValue, point_in_time);

                data_indices[value_index] = true;
                //TODO:Is it resized?
                cell_values_[cell_index].push_back(value);
            }

            void OverwritePrecalculatedValue(const std::size_t cell_index,
                                             const std::size_t value_index,
                                             const T &value) {
                if (!data_indices[value_index]) {
                    throw Exception("PrecalculatedBoundaryCondition",
                                    "Overwrite not supported");
                }
                cell_values_[cell_index][value_index] = value;
            }

            inline bool HasPrecalculatedValues(const double point_in_time) const {
                return point_in_time<=lastPrecalculatedValue;
            }

            void CalculateFastAccessIndices() {
                auto last_time = static_cast<int32_t>(LastPointinTime() * time_resolution_);
                fast_access_indices_.resize(last_time, -1);
                auto start_time = static_cast<int32_t>(point_in_times_[0] * time_resolution_);
                int32_t cur_index = 0;
                for (int32_t cur_time = start_time; cur_time < last_time; ++cur_time) {
                    const double point_in_time = static_cast<double>(cur_time) / time_resolution_;
                    if (point_in_time >= point_in_times_[cur_index + 1]) {
                        cur_index += 1;
                    }
                    fast_access_indices_[cur_time] = cur_index;
                }
            }

            inline const T neutral_element() const {
                return neutral_element_;
            }

            inline bool Indices(const double point_in_time,
                                std::size_t *index_1,
                                std::size_t *index_2) const {
                auto access = static_cast<int32_t>(point_in_time * time_resolution_);
                if (access == -1) {
                    return false;
                }
                if (access >= fast_access_indices_.size()) {
                    if (index_1 != nullptr && index_2 != nullptr) {
                        *index_1 = fast_access_indices_.back();
                        *index_2 = *index_1;
                    }
                    return true;
                }
                if (index_1 != nullptr && index_2 != nullptr) {
                    *index_1 = fast_access_indices_[access];
                    if (fast_access_indices_[access] == -1) {
                        *index_2 = *index_1;
                    } else {
                        *index_2 = *index_1 + 1;
                    }
                }
                return true;
            }

            inline const std::size_t NumberValuesPerCell() const {
                return point_in_times_.size();
            }

            inline const T Value(const std::size_t cell_index,
                                 const std::size_t value_index) const {
                if (!data_indices[value_index]) {
                    return neutral_element_;
                }
                return cell_values_[cell_index][value_index];
            }

            inline const T *ValuePtr(const std::size_t cell_index,
                                     const std::size_t value_index) const {
                if (!data_indices[value_index]) {
                    return &neutral_element_;
                }
                return &cell_values_[cell_index][value_index];
            }

            inline const double PointInTime(const std::size_t value_index) const {
                return point_in_times_[value_index];
            }

            inline const double LastPointinTime() const {
                return point_in_times_[point_in_times_.size() - 1];
            }

            inline const bool NeutralValue(const std::size_t value_index) const {
                return !data_indices[value_index];
            }

        private:
            const std::size_t total_num_cells_;
            const double time_resolution_;
            const T neutral_element_;

            bool *data_indices;
            std::vector<double> point_in_times_;
            std::vector<std::vector<T> > cell_values_;
            std::vector<int32_t> fast_access_indices_;
            double lastPrecalculatedValue;

            // NOT IMPLEMENTED
            PrecalculatedBoundaryCondition<T> &operator=(const PrecalculatedBoundaryCondition<T> &);
        };

        // adjust the phasestate of a cell
        inline void AdjustCellPhasestate(const std::size_t cell_index,
                                         const double adjustment) {
            phasestate_[cell_index] += adjustment;
        }

        // return true if a cell is crystalline
        inline bool CellCrystalline(const std::size_t index) const {
            return phasestate_[index] > 1.5;
        }

        // calculate the crystal growth speed at a given point in time
        const double CrystalGrowthSpeed(double point_in_time,
                                        std::size_t cell_index) const;

        // calculate the nucleation rate at a given point in time
        const double NucleationRate(double point_in_time,
                                    std::size_t cell_index) const;

        // integrate growth speed or nucleation rate over time
        const double IntegrateCellOverTime(double start_time,
                                           double end_time,
                                           std::size_t cell_index,
                                           bool nucleation) const;

        const Eigen::Vector3d
        RelativeCellLocation(const std::size_t absolut_index) const;

        const complex_index
        ComplexIndex(const int64_t index_x,
                     const int64_t index_y,
                     const int64_t index_z,
                     const bool check_periodic_continuation,
                     Eigen::Vector3i *periodic_continuation = NULL) const;

        const Element *
        PeriodicBoundaryNeighbour(const Eigen::Vector3i direction) const;

        const double RaytracingMethod(const std::size_t absolute_index,
                                      const complex_index spherulite_center,
                                      const Eigen::Vector3i periodic_continuation) const;

        const double
        DeltaDistance(const Eigen::Vector3d &min_cell_position,
                      const Eigen::Vector3d &direction,
                      const Eigen::Vector3d &spherulite_center_pos) const;

        void StoreIntersectionResults(const Plane &plane,
                                      const Eigen::Vector3d &direction,
                                      const Eigen::Vector3d &center,
                                      const Eigen::Vector3d &min_pos,
                                      const Eigen::Vector3d &max_pos,
                                      std::vector<Eigen::Vector3d> *results) const;

        const bool EventExists(StateEvent *event_to_check) const;

        bool VelocitiesPureElongational(const Eigen::Vector3d corner_velocities[8],
                                        const int32_t direction,
                                        const double threshold) const;

        bool VelocitiesPureShear(const Eigen::Vector3d corner_velocities[8],
                                 const int32_t main_direction,
                                 const int32_t shear_direction,
                                 const double threshold) const;

        bool CellInBorder(const size_t index) const;

        void PrecalculateIndices();

        void Initialize(const Eigen::Vector3s num_cells,
                        const InterpolationType temporal_interpolation_type,
                        const InterpolationType spatial_interpolation_type,
                        RandomNumberGenerator *rng,
                        Simulation *current_simulation);

        void ChangeCellToSolid(const std::size_t cell_index,
                               const uint64_t spherulite_id,
                               const double point_in_time,
                               const Eigen::Vector3i periodic_continuation,
                               const bool nucleation);

        bool MonteCarloStep(const double probability) const;

        const double IntegrateCellToValue(const std::size_t cell_index,
                                          const double dest_value,
                                          const double start_time,
                                          const bool nucleation) const;

        Eigen::Matrix3d DeformationTensor(const std::size_t cell_index,
                                          const double time_1,
                                          const double time_2,
                                          bool *consider_pure_elongational) const;

        const double CellPointElongationLn(const std::size_t cell_index,
                                           const double start_time,
                                           const double end_time) const;

        const double CellPointShearing(const std::size_t cell_index,
                                       const double start_time,
                                       const double end_time) const;

        bool FlowPureShearElongational(const double start_time,
                                       const double end_time,
                                       bool *consider_elongational,
                                       bool *consider_shearing,
                                       bool *no_flow) const;

        const double DimensionlessInduction(const double point_in_time,
                                            const std::size_t cell_index,
                                            const bool flow,
                                            const bool thermal) const;

        void solveDifferentialEquations(Eigen::Matrix3d &velGradient,
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
                                        double &dst3) const;

        template<class T>
        const T BoundaryValue(const std::size_t cell_index,
                              const double current_time,
                              const BoundaryCondition<T> *boundary_condition,
                              const T *neighbour_values[6]) const {
            T corner_values[8];
            boundary_condition->CornerValues(current_time, corner_values);

            const Eigen::Vector3d relative_offset = RelativeCellLocation(cell_index);

            // get the values from the neighbours
            const Eigen::Vector3d neighbour_sizes[] = {
                    neighbour_values[0] != NULL ?
                    neighbour_infos_[0].first->GetSizeM() : Eigen::Vector3d(0.0, 0.0, 0.0),
                    neighbour_values[1] != NULL ?
                    neighbour_infos_[1].first->GetSizeM() : Eigen::Vector3d(0.0, 0.0, 0.0),
                    neighbour_values[2] != NULL ?
                    neighbour_infos_[2].first->GetSizeM() : Eigen::Vector3d(0.0, 0.0, 0.0),
                    neighbour_values[3] != NULL ?
                    neighbour_infos_[3].first->GetSizeM() : Eigen::Vector3d(0.0, 0.0, 0.0),
                    neighbour_values[4] != NULL ?
                    neighbour_infos_[4].first->GetSizeM() : Eigen::Vector3d(0.0, 0.0, 0.0),
                    neighbour_values[5] != NULL ?
                    neighbour_infos_[5].first->GetSizeM() : Eigen::Vector3d(0.0, 0.0, 0.0),
            };

            return boundary_condition->SpatialInterpolation(relative_offset,
                                                            corner_values,
                                                            GetSizeM(),
                                                            neighbour_values,
                                                            neighbour_sizes);
        }

        template<class T>
        const T PrecalculatedBoundaryValue(const std::size_t cell_index,
                                           const double current_time,
                                           PrecalculatedBoundaryCondition<T> *precalc,
                                           BoundaryCondition<T> *boundary_cond) const {
            // retrieve the two indices in the precalculated value storage
            std::size_t index_1, index_2;
            if (!precalc->Indices(current_time, &index_1, &index_2)) {
                return precalc->neutral_element();
            }
            if (index_1 == index_2) {
                return precalc->Value(cell_index, index_1);
            }

            typedef typename Interpolation<T>::ref_point rp;
            typedef typename Interpolation<T>::ref_point_ptr rp_ptr;

            // set interpolation-points
            rp p1(precalc->PointInTime(index_1),
                  precalc->Value(cell_index, index_1));
            rp p2(precalc->PointInTime(index_2),
                  precalc->Value(cell_index, index_2));
            rp_ptr before(index_1 > 0 ?
                          precalc->PointInTime(index_1 - 1) :
                          0.0,
                          index_1 > 0 ?
                          precalc->ValuePtr(cell_index, index_1 - 1) :
                          nullptr);
            rp_ptr after(index_2 < velocities_->NumberValues() - 1 ?
                         precalc->PointInTime(index_2 + 1) :
                         0.0,
                         index_2 < velocities_->NumberValues() - 1 ?
                         precalc->ValuePtr(cell_index, index_2 + 1) :
                         nullptr);

            // interpolation for the correct point in time
            Interpolation<T> inter(boundary_cond->TemporalInterpolationType());
            return inter.ExecuteInterpolation(p1, p2, before, after, current_time);
        }

        // fast access tables
        static const int32_t fast_table_resolution = 100000;
        static bool fast_access_table_initialized;
        static double table_values_sin[fast_table_resolution];
        static double table_values_cos[fast_table_resolution];
        static double table_values_log[fast_table_resolution];

        static void InitFastAccessTables();

        // element properties
        const uint64_t element_id_;

        // first entry: pointer to the element
        // second entry: true if the neighbour element should exist (first entry can still
        //               be zero in case the data does not exist)
        //               false if the neighbour element does not exist and has only 0.0-values
        //               for temperature and velocity
        std::vector<std::pair<const Element *, const bool> > neighbour_infos_;

        // areal properties
        const Eigen::Vector3d lower_left_front_corner_m_;
        const double cell_size_m_;
        const GrowthEvent::Type growth_type_;

        Eigen::Vector3d size_m_;
        Eigen::Vector3s num_cells_;
        BoundaryCondition<double> *temperatures_;
        BoundaryCondition<Eigen::Vector3d> *velocities_;
        bool has_written_final_results_;

        std::vector<double> phasestate_;
        std::vector<double> crystallization_time_;
        AdditionalInformationCCG *additional_information_ccg_;
        std::size_t num_crystalline_cells_;

        RandomNumberGenerator *rng_;
        Simulation *current_simulation_;

        // precomputed values and caching
        std::vector<Eigen::Vector3s> precalculated_index_vectors_;
        PrecalculatedBoundaryCondition<double> *precalculated_temperatures_;
        PrecalculatedBoundaryCondition<Eigen::Vector3d> *precalculated_velocities_;
        PrecalculatedBoundaryCondition<bool> *precalculated_shear_flags_;
        PrecalculatedBoundaryCondition<bool> *precalculated_elongation_flags_;

        // state events
        std::multiset<StateEvent *, StateEventCompare> state_events_;

        // cached values
        std::vector<double> last_nucleation_rates_;

        // NOT IMPLEMENTED
        Element &operator=(const Element &);

        //constants
    public:
        inline constexpr static const char *VTK_TEMPERATUE = "Temperature [C]";
        inline constexpr static const char *VTK_PHASESTATE = "Phasestate [-]";
        inline constexpr static const char *VTK_DISENGAGEMENT_TIME = "Disengagement time [s]";
        inline constexpr static const char *VTK_SPHERULITEID = "spheru_id";
        inline constexpr static const char *VTK_CRYSTALLIZATION_TIME = "Crystallization time [s]";
        inline constexpr static const char *VTK_SPHERULITE_DIAMETER = "Spherulite diameter [um]";
        inline constexpr static const char *VTK_RANDOM_COLOR = "Random color [-]";
        inline constexpr static const char *VTK_ASPECT_RATIO = "Aspect ratio [-]";
        inline constexpr static const char *VTK_ENTHROPY_CHANGE = "Enthropy_change [J/K]";
        inline constexpr static const char *VTK_INDUCTION_FLOW = "Induction flow [-]";
        inline constexpr static const char *VTK_INDUCTION_ATHERMAL = "Induction athermal [-]";
        inline constexpr static const char *VTK_INV_ATHERMAL = "Inv. athermal [-]";
        inline constexpr static const char *VTK_INDUCTION_ATHERMAL_FLOW = "Induction athermal+flow [-]";
        inline constexpr static const char *VTK_COOLING_RATE = "Cooling rate [K/s]";
        inline constexpr static const char *VTK_Location = "Location [mm]";
        inline constexpr static const char *VTK_Velocity = "Velocity [mm/s]";
    };
}  // namespace SphaeroSim

#endif  // _SPHAEROSIM_ELEMENT_H_