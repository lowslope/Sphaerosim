#include "SpheruliteInfo.h"

#ifndef _INCLUDED_STDINT_H_
#include <stdint.h>
#define _INCLUDED_STDINT_H_
#endif

#ifndef _INCLUDED_VECTOR_H_
#include <vector>
#define _INCLUDED_VECTOR_H_
#endif

#ifndef _SPHAEROSIM_STATE_EVENT_H_
#include "StateEvent.h"
#endif

#ifndef _SPHAEROSIM_ELEMENT_H_
#include "Element.h"
#endif

namespace SphaeroSim {

    SpheruliteInfo::SpheruliteInfo(const uint64_t spherulite_id,
                                   const std::size_t start_location,
                                   const double creation_time,
                                   const int32_t random_color,
                                   const Element *owner) :
            start_location_(start_location),
            creation_time_(creation_time),
            spherulite_id_(spherulite_id),
            random_color_(random_color),
            owner_(owner) {

        bb_min_ = Eigen::Vector3d(1e100, 1e100, 1e100);
        bb_max_ = Eigen::Vector3d(-1e100, -1e100, -1e100);
        number_cells_in_spherulite_ = 1;

        // sqrt(3/2) * CellSize
        spherulite_diameter_ = 1.2247448 * owner->GetCellSizeM();
    }

    const double SpheruliteInfo::AspectRatio(const double cellsize_m) const {
        const double side_1 = std::abs(bb_max_[0] - bb_min_[0]);
        const double side_2 = std::abs(bb_max_[1] - bb_min_[1]);
        const double side_3 = std::abs(bb_max_[2] - bb_min_[2]);

        double max_side = 0.0;
        if (side_1 > side_2) {
            if (side_1 > side_3) {
                max_side = side_1;
            } else {
                max_side = side_3;
            }
        } else {
            if (side_2 > side_3) {
                max_side = side_2;
            } else {
                max_side = side_3;
            }
        }

        double min_side = 0.0;
        if (side_1 < side_2) {
            if (side_1 < side_3) {
                min_side = side_1;
            } else {
                min_side = side_3;
            }
        } else {
            if (side_2 < side_3) {
                min_side = side_2;
            } else {
                min_side = side_3;
            }
        }

        if (min_side < cellsize_m) {
            min_side = cellsize_m;
        }
        if (max_side < cellsize_m) {
            max_side = cellsize_m;
        }
        return max_side / min_side;
    }

    const double SpheruliteInfo::EquivalentDiameter() const {
        const double volume = number_cells_in_spherulite_ *
                              GetOwnerElement()->GetCellSizeM() *
                              GetOwnerElement()->GetCellSizeM() *
                              GetOwnerElement()->GetCellSizeM();

        const double radius = pow(3.0 * volume / (4.0 * M_PI), 1.0 / 3.0);
        return 2.0 * radius;
    }

    void SpheruliteInfo::AddCellToSpherulite(const complex_index &new_cell,
                                             const Eigen::Vector3i &periodic_continuation) {
        const Eigen::Vector3d spherulite_centre =
                owner_->CellPosition(start_location_, true);
        Eigen::Vector3d cell_centre =
                new_cell.first->CellPosition(new_cell.second, true);

        const Eigen::Vector3s sph_indices = owner_->CellIndices(start_location_);
        const Eigen::Vector3s cell_indices = new_cell.first->CellIndices(new_cell.second);

        const Element *min_neighbour[3];
        const Element *max_neighbour[3];
        new_cell.first->MinMaxNeighbour(min_neighbour, max_neighbour);
        cell_centre = new_cell.first->PeriodicCorrection(-1 * periodic_continuation,
                                                         min_neighbour,
                                                         max_neighbour,
                                                         cell_centre,
                                                         spherulite_centre);
        #pragma omp critical
        {
            for (unsigned short int cur = 0; cur < 3; ++cur) {
                if (cell_centre[cur] < bb_min_[cur]) {
                    bb_min_[cur] = cell_centre[cur];
                }
                if (cell_centre[cur] > bb_max_[cur]) {
                    bb_max_[cur] = cell_centre[cur];
                }
            }

            for (unsigned short int cur = 0; cur < 3; ++cur) {
                if (bb_max_[cur] - bb_min_[cur] > spherulite_diameter_) {
                    spherulite_diameter_ = bb_max_[cur] - bb_min_[cur];
                }
            }

            if (spherulite_diameter_ < owner_->GetCellSizeM()) {
                spherulite_diameter_ = owner_->GetCellSizeM();
            }

            number_cells_in_spherulite_ += 1;
        }

        /*
        const double new_radius_sqr =
          (cell_centre - spherulite_centre).squaredNorm();
        const double old_radius_sqr =
          0.25 * spherulite_diameter_ * spherulite_diameter_;
        if (new_radius_sqr > old_radius_sqr)
          spherulite_diameter_ = 2.0 * sqrt(new_radius_sqr);
          */
    }
}  // namespace SphaeroSim
