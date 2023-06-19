#ifndef _SPHAEROSIM_BOUNDARY_CONDITION_H_
#define _SPHAEROSIM_BOUNDARY_CONDITION_H_

#ifndef _INCLUDED_SET_H_
#include <set>
#define _INCLUDED_SET_H_
#endif

#ifndef _INCLUDED_VECTOR_H_
#include <vector>
#define _INCLUDED_VECTOR_H_
#endif

#ifndef _INCLUDED_UTILITY_H_
#include <utility>
#define _INCLUDED_UTILITY_H_
#endif

#ifndef _INCLUDED_STDINT_H_
#include <stdint.h>
#define _INCLUDED_STDINT_H_
#endif

#ifndef _SPHAEROSIM_INTERPOLATION_H_
#include "Interpolation.h"
#endif

#ifndef _SPHAEROSIM_EXCEPTION_H_
#include "Exception.h"
#endif

namespace SphaeroSim {

    template<class T>
    class BoundaryCondition {
    public:
        BoundaryCondition(Interpolation <T> *temporal_interpolation,
                          Interpolation <T> *spatial_interpolation) {
            temporal_interpolation_ = temporal_interpolation;
            spatial_interpolation_ = spatial_interpolation;
        }

        ~BoundaryCondition() {
            typename std::set<Entry *>::iterator it = entries_.begin();
            while (it != entries_.end()) {
                delete (*it);
                ++it;
            }
            delete temporal_interpolation_;
            delete spatial_interpolation_;
        }

        // retrieve all point of times of the entries < time
        void PointInTimeOfEntries(const double time,
                                  std::vector<double> *dst) const {
            typename std::set<Entry *, EntryComperator>::const_iterator it = entries_.begin();
            while ((*it)->GetPointInTime() < time) {
                dst->push_back((*it)->GetPointInTime());
                ++it;
                if (it == entries_.end()) {
                    break;
                }
            }
        }

        // retrieve the first point in time for which data exists
        double FirstPointInTime() const {
            return (*entries_.begin())->GetPointInTime();
        }

        // add an entry of the corner values
        inline void AddEntry(const double time,
                             const std::size_t index,
                             const std::vector<T> &values) {
            entries_.insert(new Entry(time, index, values));
        }

        // return true if corner values exist for the point in time
        inline bool ValuesExist(const double time) const {
            if (entries_.size() == 0) {
                throw Exception("BoundaryCondition",
                                "No values stored in the boundary conditions.");
            }
            return time >= (*entries_.begin())->GetPointInTime();
        }

        // retrieve the indices of the values before and after a point in time
        void EntryIndicesBeforeAfter(const double current_time,
                                     std::size_t *index_1,
                                     std::size_t *index_2) const {
            if (entries_.size() == 1) {
                *index_1 = 0;
                *index_2 = 0;
                return;
            }

            // get the right interval
            Entry tmp_entry(current_time, 0, std::vector<T>());
            typename std::set<Entry *>::iterator first_entry = entries_.upper_bound(&tmp_entry);
            if (first_entry == entries_.end()) {
                *index_1 = entries_.size() - 1;
                *index_2 = entries_.size() - 1;
                return;
            }

            typename std::set<Entry *>::iterator second_entry = first_entry;
            if (first_entry != entries_.begin()) {
                --first_entry;
            }
            if (second_entry == entries_.end()) {
                second_entry = first_entry;
            }

            *index_1 = (*first_entry)->GetIndex();
            *index_2 = (*second_entry)->GetIndex();
        }

        // get the corner values for a specific point in time
        void CornerValues(const double time,
                          T results[8]) const {
            if (entries_.size() == 0) {
                throw Exception("BoundaryCondition",
                                "No values stored in the boundary conditions.");
            }
            if (time < (*entries_.begin())->GetPointInTime()) {
                throw Exception("BoundaryConditions",
                                "No values for the requested point in time");
            }

            if (time > (*entries_.rbegin())->GetPointInTime()) {
                for (std::size_t cur = 0; cur < 8; ++cur)
                    results[cur] = (*entries_.rbegin())->GetCornerValue(cur);
                return;
            }
            if (entries_.size() == 1) {
                for (std::size_t cur = 0; cur < 8; ++cur)
                    results[cur] = (*entries_.begin())->GetCornerValue(cur);
                return;
            }

            // get the right interval
            Entry tmp_entry(time, 0, std::vector<T>());
            typename std::set<Entry *>::iterator first_entry = entries_.upper_bound(&tmp_entry);
            if (first_entry == entries_.end()) {
                for (std::size_t cur = 0; cur < 8; ++cur)
                    results[cur] = (*entries_.rbegin())->GetCornerValue(cur);
                return;
            }

            typename std::set<Entry *>::iterator second_entry = first_entry;
            if (first_entry != entries_.begin()) {
                --first_entry;
            }
            if (second_entry == entries_.end()) {
                second_entry = first_entry;
            }

            if (first_entry == second_entry) {
                for (std::size_t cur = 0; cur < 8; ++cur)
                    results[cur] = (*first_entry)->GetCornerValue(cur);
                return;
            }

            typename std::set<Entry *>::iterator entry_before = entries_.end();
            if (first_entry != entries_.begin()) {
                entry_before = first_entry;
                --entry_before;
            }

            typename std::set<Entry *>::iterator entry_after = second_entry;
            entry_after++;
            if (entry_after == entries_.end()) {
                entry_after = entries_.begin();
            }

            // interpolate the values
            for (std::size_t cur = 0; cur < 8; ++cur) {
                typename Interpolation<T>::ref_point p1((*first_entry)->GetPointInTime(),
                                                        (*first_entry)->GetCornerValue(cur));
                typename Interpolation<T>::ref_point p2((*second_entry)->GetPointInTime(),
                                                        (*second_entry)->GetCornerValue(cur));

                typename Interpolation<T>::ref_point_ptr before(
                        entry_before != entries_.end() ?
                        (*first_entry)->GetPointInTime() - (*entry_before)->GetPointInTime() :
                        0.0,
                        entry_before != entries_.end() ?
                        (*entry_before)->GetCornerValuePtr(cur) :
                        NULL);
                typename Interpolation<T>::ref_point_ptr after(
                        entry_after != entries_.begin() ?
                        (*entry_after)->GetPointInTime() - (*second_entry)->GetPointInTime() :
                        0.0,
                        entry_after != entries_.begin() ?
                        (*entry_after)->GetCornerValuePtr(cur) :
                        NULL);

                results[cur] = temporal_interpolation_->ExecuteInterpolation(p1,
                                                                             p2,
                                                                             before,
                                                                             after,
                                                                             time);
            }
        }

        // get spatial interpolation for a specific coordinate
        // neighbour_values:
        // first_entry: relative size of the element compared to the current one
        //              if this value is < 0.0, no neighbouring element is specified
        // second_entry: corner values of the neighbouring element
        const T
        SpatialInterpolation(const Eigen::Vector3d relative_offset,
                             const T corner_values[8],
                             const Eigen::Vector3d size_m,
                             const T *neighbours[6],
                             const Eigen::Vector3d neighbour_size_m[6]) const {
            Eigen::Vector3d rel_neigh_size[6];
            for (std::size_t cur = 0; cur < 6; ++cur)
                rel_neigh_size[cur] = Eigen::Vector3d(
                        neighbour_size_m[cur][0] / size_m[0],
                        neighbour_size_m[cur][1] / size_m[1],
                        neighbour_size_m[cur][2] / size_m[2]);
            //                       6--------7
            //                      /|       /|
            //                     2--------3 |                  y
            //                     | |      | |                  |  z
            //                     | 4------|-5                  | /
            //                     |/       |/                   |/
            //                     0--------1                    ------x
            // the neighbour elements (+/- X_AXIS; +/- Y_AXIS; +/- Z_AXIS)
            // easy access
            typedef typename Interpolation<T>::ref_point point;
            typedef typename Interpolation<T>::ref_point_ptr point_ptr;
            // interpolation for +/- Z_AXIS
            // 0 -> 4
            const T val04 = spatial_interpolation_->ExecuteInterpolation(
                    point(0.0, corner_values[0]),
                    point(1.0, corner_values[4]),
                    point_ptr(neighbours[5] != NULL ? rel_neigh_size[5][2] : 0.0,
                              neighbours[5] != NULL ? &neighbours[5][0] : NULL),
                    point_ptr(neighbours[4] != NULL ? rel_neigh_size[4][2] : 0.0,
                              neighbours[4] != NULL ? &neighbours[4][4] : NULL),
                    relative_offset[2]);

            // 1 -> 5
            const T val15 = spatial_interpolation_->ExecuteInterpolation(
                    point(0.0, corner_values[1]),
                    point(1.0, corner_values[5]),
                    point_ptr(neighbours[5] != NULL ? rel_neigh_size[5][2] : 0.0,
                              neighbours[5] != NULL ? &neighbours[5][1] : NULL),
                    point_ptr(neighbours[4] != NULL ? rel_neigh_size[4][2] : 0.0,
                              neighbours[4] != NULL ? &neighbours[4][5] : NULL),
                    relative_offset[2]);

            // 2 -> 6
            const T val26 = spatial_interpolation_->ExecuteInterpolation(
                    point(0.0, corner_values[2]),
                    point(1.0, corner_values[6]),
                    point_ptr(neighbours[5] != NULL ? rel_neigh_size[5][2] : 0.0,
                              neighbours[5] != NULL ? &neighbours[5][2] : NULL),
                    point_ptr(neighbours[4] != NULL ? rel_neigh_size[4][2] : 0.0,
                              neighbours[4] != NULL ? &neighbours[4][6] : NULL),
                    relative_offset[2]);

            // 3 -> 7
            const T val37 = spatial_interpolation_->ExecuteInterpolation(
                    point(0.0, corner_values[3]),
                    point(1.0, corner_values[7]),
                    point_ptr(neighbours[5] != NULL ? rel_neigh_size[5][2] : 0.0,
                              neighbours[5] != NULL ? &neighbours[5][3] : NULL),
                    point_ptr(neighbours[4] != NULL ? rel_neigh_size[4][2] : 0.0,
                              neighbours[4] != NULL ? &neighbours[4][7] : NULL),
                    relative_offset[2]);

            // interpolation for neighbours +/- Y_AXIS
            T interMinusY_04 = val04;
            T interMinusY_15 = val15;
            if (neighbours[3] != NULL) {
                interMinusY_04 = neighbours[3][0] * (1.0 - relative_offset[2]) +
                                 neighbours[3][4] * relative_offset[2];
                interMinusY_15 = neighbours[3][1] * (1.0 - relative_offset[2]) +
                                 neighbours[3][5] * relative_offset[2];
            }

            T interPlusY_04 = val26;
            T interPlusY_15 = val37;
            if (neighbours[2] != NULL) {
                interPlusY_04 = neighbours[2][2] * (1.0 - relative_offset[2]) +
                                neighbours[2][6] * relative_offset[2];
                interPlusY_15 = neighbours[2][3] * (1.0 - relative_offset[2]) +
                                neighbours[2][7] * relative_offset[2];
            }

            // 04 -> 26
            const T val04_26 = spatial_interpolation_->ExecuteInterpolation(
                    point(0.0, val04),
                    point(1.0, val26),
                    point_ptr(neighbours[3] != NULL ? rel_neigh_size[3][1] : 0.0,
                              neighbours[3] != NULL ? &interMinusY_04 : NULL),
                    point_ptr(neighbours[2] != NULL ? rel_neigh_size[2][1] : 0.0,
                              neighbours[2] != NULL ? &interPlusY_04 : NULL),
                    relative_offset[1]);

            // 15 -> 37
            const T val15_37 = spatial_interpolation_->ExecuteInterpolation(
                    point(0.0, val15),
                    point(1.0, val37),
                    point_ptr(neighbours[3] != NULL ? rel_neigh_size[3][1] : 0.0,
                              neighbours[3] != NULL ? &interMinusY_15 : NULL),
                    point_ptr(neighbours[2] != NULL ? rel_neigh_size[2][1] : 0.0,
                              neighbours[2] != NULL ? &interPlusY_15 : NULL),
                    relative_offset[1]);

            // interpolation for neighbours +/- X_AXIS
            T interMinusX = val04_26;
            if (neighbours[1] != NULL) {
                T tmp1 = neighbours[1][0] * (1.0 - relative_offset[2]) +
                         neighbours[1][4] * relative_offset[2]; // 0 -> 4
                T tmp2 = neighbours[1][2] * (1.0 - relative_offset[2]) +
                         neighbours[1][6] * relative_offset[2]; // 2 -> 6
                interMinusX = tmp1 * (1.0 - relative_offset[1]) + tmp2 * relative_offset[1];
            }
            T interPlusX = val15_37;
            if (neighbours[0] != NULL) {
                T tmp1 = neighbours[0][1] * (1.0 - relative_offset[2]) +
                         neighbours[0][5] * relative_offset[2]; // 1 -> 5
                T tmp2 = neighbours[0][3] * (1.0 - relative_offset[2]) +
                         neighbours[0][7] * relative_offset[2]; // 3 -> 7
                interPlusX = tmp1 * (1.0 - relative_offset[1]) + tmp2 * relative_offset[1];
            }

            // final interpolation
            const T result = spatial_interpolation_->ExecuteInterpolation(
                    point(0.0, val04_26),
                    point(1.0, val15_37),
                    point_ptr(neighbours[1] != NULL ? rel_neigh_size[1][0] : 0.0,
                              neighbours[1] != NULL ? &interMinusX : NULL),
                    point_ptr(neighbours[0] != NULL ? rel_neigh_size[0][0] : 0.0,
                              neighbours[0] != NULL ? &interPlusX : NULL),
                    relative_offset[0]);
            return result;
        }

        // getter
        inline const InterpolationType TemporalInterpolationType() const {
            return temporal_interpolation_->GetType();
        }

        inline const std::size_t NumberValues() const {
            return entries_.size();
        }

    private:
        class Entry {
        public:
            Entry(const double time,
                  const std::size_t index,
                  const std::vector<T> &values) {
                point_in_time_ = time;
                corner_values_ = values;
                index_ = index;
            }

            // Getter
            inline const double GetPointInTime() const {
                return point_in_time_;
            }

            inline const T GetCornerValue(const std::size_t index) const {
                return corner_values_[index];
            }

            inline const T *GetCornerValuePtr(const std::size_t index) const {
                return &corner_values_[index];
            }

            inline const std::vector<T> *GetCornerValues() const {
                return &corner_values_;
            }

            inline const std::size_t GetIndex() const {
                return index_;
            }

        private:
            double point_in_time_;
            std::vector<T> corner_values_;
            std::size_t index_;
        };

        class EntryComperator {
        public:
            bool operator()(const Entry *entry_1,
                            const Entry *entry_2) const {
                return entry_1->GetPointInTime() < entry_2->GetPointInTime();
            }
        };

        std::set<Entry *, EntryComperator> entries_;

        Interpolation <T> *temporal_interpolation_;
        Interpolation <T> *spatial_interpolation_;
    };

}  // namespace SphaeroSim

#endif  // _SPHAEROSIM_BOUNDARY_CONDITION_H_