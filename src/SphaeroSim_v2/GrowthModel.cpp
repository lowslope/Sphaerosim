#include "GrowthModel.h"

#ifndef _INCLUDED_STDINT_H_
#include <stdint.h>
#define _INCLUDED_STDINT_H_
#endif

#ifndef _INCLUDED_VECTOR_H_
#include <vector>
#define _INCLUDED_VECTOR_H_
#endif

#ifndef _SPHAEROSIM_EXCEPTION_H_

#include "Exception.h"

#endif

#ifndef _SPHAEROSIM_ELEMENT_H_

#include "Element.h"

#endif

#ifndef _SPHAEROSIM_MATERIAL_H_

#include "Material.h"

#endif

#ifndef _SPHAEROSIM_SIMULATION_H_

#include "Simulation.h"

#endif

namespace SphaeroSim {

    GrowthModel::GrowthModel(const GrowthEvent::Type growth_type) :
            growth_type_(growth_type) {
    }

    ConstantGrowthModel::ConstantGrowthModel(const GrowthEvent::Type growth_type,
                                             const double constant_speed) :
            GrowthModel(growth_type),
            constant_speed_(constant_speed) {

        if (constant_speed_ <= 0.0) {
            throw Exception("GrowthModel",
                            "The growth speed should be larger than zero.");
        }
    }

    const double ConstantGrowthModel::Speed(const double point_in_time,
                                            const std::size_t cell_index,
                                            const Element *element) const {
        point_in_time;  // not needed
        cell_index;  // not needed
        element;  // not needed
        return constant_speed_;
    }

    const double ConstantGrowthModel::ScalingFactor() const {
        return 1.0 / constant_speed_;
    }

    HDLGrowthModel::HDLGrowthModel(const GrowthEvent::Type growth_type,
                                   const Material *material,
                                   const double G0,
                                   const double U0,
                                   const double K_g) :
            GrowthModel(growth_type),
            G0_(G0),
            U0_(U0),
            K_g_(K_g) {

        // calculate the scaling factor
        double current_temperature = 273.15;
        double maximum_temperature = current_temperature;
        double maximum_growth_speed = GrowthSpeed(material,
                                                  current_temperature);
        while (current_temperature < 573.15) {
            const double current_growth_speed = GrowthSpeed(material,
                                                            current_temperature);
            if (current_growth_speed > maximum_growth_speed) {
                maximum_growth_speed = current_growth_speed;
                maximum_temperature = current_temperature;
            }
            current_temperature += 5.0;
        }
        scaling_factor_ = 1.0 / maximum_growth_speed;
    }

    const double HDLGrowthModel::GrowthSpeed(const Material *material,
                                             const double temperature) const {
        const double T_m0 = material->GetTm0();
        const double T_inf = material->GetTinf();

        if (temperature > T_m0) {
            return 0.0;
        }
        if (temperature < T_inf) {
            return 0.0;
        }
        // calculate growth speed
        const double rF = 2.0 * temperature / (T_m0 + temperature);
        //double exp1 = exp(-U0_/(8.314472 * (temperature - T_inf)));
        //double exp2 = exp(-K_g_*T_m0*T_m0/(rF * temperature * (T_m0 - temperature)));
        //double rGrowthSpeed = G0_*exp1*exp2; // growth speed for this cell
        const double exp_value = exp(-K_g_ * T_m0 * T_m0 / (rF * temperature * (T_m0 - temperature)) +
                                     -U0_ / (8.314472 * (temperature - T_inf)));
        const double rGrowthSpeed = G0_ * exp_value;
        return rGrowthSpeed;  // growth speed for this cell
    }

    const double HDLGrowthModel::Speed(const double point_in_time,
                                       const std::size_t cell_index,
                                       const Element *element) const {
        return GrowthSpeed(element->GetCurrentSimulation()->GetMaterial(),
                           element->CellTemperature(cell_index, point_in_time));
    }

    const double HDLGrowthModel::ScalingFactor() const {
        return scaling_factor_;
    }
}  // namespace SphaeroSim
