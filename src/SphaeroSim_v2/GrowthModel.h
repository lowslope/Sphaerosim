#ifndef _SPHAEROSIM_GROWTH_MODEL_H_
#define _SPHAEROSIM_GROWTH_MODEL_H_

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

namespace SphaeroSim {

    class Element;

    class Material;

    class GrowthModel {
    public:
        virtual const double Speed(const double point_in_time,
                                   const std::size_t cell_index,
                                   const Element *element) const = 0;

        virtual const double ScalingFactor() const = 0;

        inline const GrowthEvent::Type GetGrowthType() const {
            return growth_type_;
        }

    protected:
        explicit GrowthModel(const GrowthEvent::Type growth_type);

    private:
        const GrowthEvent::Type growth_type_;

        // not implemented
        GrowthModel &operator=(const GrowthModel &);
    };

    class ConstantGrowthModel : public GrowthModel {
    public:
        ConstantGrowthModel(const GrowthEvent::Type growth_type,
                            const double constant_speed);

        const double Speed(const double point_in_time,
                           const std::size_t cell_index,
                           const Element *element) const;

        const double ScalingFactor() const;

    private:
        const double constant_speed_;

        // not implemented
        ConstantGrowthModel &operator=(const ConstantGrowthModel &);
    };

    class HDLGrowthModel : public GrowthModel {
    public:
        HDLGrowthModel(const GrowthEvent::Type growth_type,
                       const Material *material,
                       const double G0,
                       const double U0,
                       const double K_g);

        const double Speed(const double point_in_time,
                           const std::size_t cell_index,
                           const Element *element) const;

        const double ScalingFactor() const;

    private:
        const double GrowthSpeed(const Material *material,
                                 const double temperature) const;

        const double G0_;
        const double U0_;
        const double K_g_;
        double scaling_factor_;

        // not implemented
        HDLGrowthModel &operator=(const HDLGrowthModel &);
    };

}  // namespace SphaeroSim

#endif  // _SPHAEROSIM_GROWTH_MODEL_H_