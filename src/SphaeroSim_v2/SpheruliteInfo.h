#ifndef _SPHAEROSIM_SPHERULITE_INFO_H_
#define _SPHAEROSIM_SPHERULITE_INFO_H_

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

    class Element;

    class SpheruliteInfo {
    public:
        SpheruliteInfo(const uint64_t spherulite_id,
                       const std::size_t start_location,
                       const double creation_time,
                       const int32_t random_color,
                       const Element *owner);

        void AddCellToSpherulite(const complex_index &new_cell,
                                 const Eigen::Vector3i &periodic_continuation);

        const double AspectRatio(const double cellsize_m) const;

        const double EquivalentDiameter() const;

        // Getter
        inline const std::size_t GetStartLocation() const {
            return start_location_;
        }

        inline const Element *GetOwnerElement() const {
            return owner_;
        }

        inline const double GetDiameter() const {
            return spherulite_diameter_;
        }

        inline const int32_t GetRandomColor() const {
            return random_color_;
        }

        inline const double GetCreationTime() const {
            return creation_time_;
        }

        inline const uint64_t GetID() const {
            return spherulite_id_;
        }

        inline const Eigen::Vector3d GetStartLocationAbsolute() const {
            return owner_->CellPosition(start_location_, true);
        }

        inline const Eigen::Vector3d GetBoxCenter() const {
            return 0.5 * (bb_max_ + bb_min_);
        }

        inline const Eigen::Vector3d GetBoxSize() const {
            return bb_max_ - bb_min_;
        }

        inline const int32_t GetNumberCellsInSpherulite() const {
            return number_cells_in_spherulite_;
        }

    private:
        const std::size_t start_location_;
        const double creation_time_;
        const uint64_t spherulite_id_;
        const Element *owner_;
        const int32_t random_color_;

        double spherulite_diameter_;
        Eigen::Vector3d bb_min_;  // bounding box minimum
        Eigen::Vector3d bb_max_;  // bounding box maximum
        int32_t number_cells_in_spherulite_;

        // not implemented
        SpheruliteInfo &operator=(const SpheruliteInfo &);
    };

}  // namespace SphaeroSim

#endif