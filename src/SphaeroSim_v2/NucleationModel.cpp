#include "NucleationModel.h"

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

    NucleationModel::NucleationModel() = default;

    ConstantNucleationModel::ConstantNucleationModel(const double constant_rate) :
            constant_rate_(constant_rate) {

        if (constant_rate_ <= 0.0) {
            throw Exception("NucleationModel",
                            "The nucleation rate should be larger than zero.");
        }
    }

    const double ConstantNucleationModel::Rate(const double point_in_time,
                                               const std::size_t cell_index,
                                               const bool quiescent_only,
                                               const bool thermal_only,
                                               const Element *element) const {
        point_in_time;  // not needed
        cell_index;  // not needed
        element;  // not needed
        quiescent_only;  // no needed
        thermal_only;  // no needed
        return constant_rate_;
    }

    const double ConstantNucleationModel::ScalingFactor() const {
        return 1.0 / constant_rate_;
    }

    HDLNucleationModel::HDLNucleationModel(const Material *material,
                                           const double C,
                                           const double U,
                                           const double K,
                                           const double n,
                                           const bool enthalpy_correction,
                                           const bool athermal,
                                           const double C_ath) :
            C_(C),
            U_(U),
            K_(K),
            n_(n),
            athermal_(athermal),
            C_ath_(AthermalConstant(C_ath)),
            enthalpy_correction_(enthalpy_correction) {

        // calculate the scaling factor
        double current_temperature = 273.15;
        double maximum_temperature = current_temperature;
        double maximum_nucleation_rate = CalculateNucleationRate(material,
                                                                 current_temperature,
                                                                 0.0, 0.0, 0.0);
        while (current_temperature < 573.15) {
            const double current_nucleation_rate =
                    CalculateNucleationRate(material,
                                            current_temperature,
                                            0.0, 0.0, 0.0);
            if (current_nucleation_rate > maximum_nucleation_rate) {
                maximum_nucleation_rate = current_nucleation_rate;
                maximum_temperature = current_temperature;
            }
            current_temperature += 5.0;
        }
        scaling_factor_ = 1.0 / maximum_nucleation_rate;
    }

    const double HDLNucleationModel::Rate(const double point_in_time,
                                          const std::size_t cell_index,
                                          const bool quiescent_only,
                                          const bool thermal_only,
                                          const Element *element) const {
        double coolingrate = 0.0;
        if (!thermal_only) {
            coolingrate = element->CellCoolingRate(cell_index, point_in_time);
        }
        if (quiescent_only) {
            return CalculateNucleationRate(element->GetCurrentSimulation()->GetMaterial(),
                                           element->CellTemperature(cell_index, point_in_time),
                                           coolingrate,
                                           0.0,
                                           0.0);
        } else {
            double flow_enthropy = 0.0;
            double flow_energy = element->CellFlowEnergy(cell_index, point_in_time, &flow_enthropy);
            return CalculateNucleationRate(element->GetCurrentSimulation()->GetMaterial(),
                                           element->CellTemperature(cell_index, point_in_time),
                                           coolingrate,
                                           flow_energy,
                                           flow_enthropy);
        }
    }

    const double HDLNucleationModel::ScalingFactor() const {
        return scaling_factor_;
    }

    const double HDLNucleationModel::AthermalConstant(const double volume) const {
        if (volume <= 0.0) {
            return 0.0;
        }

        const double planck = 6.626E-34;
        const double boltzmann = 1.38E-23;
        const double i = n_ + 1.0;
        const double factor_1 = planck / (2 * i) * i * i / (i - 1.0);
        const double factor_2 = K_ * boltzmann * (i - 1.0) / volume;
        return C_ * factor_1 * pow(factor_2, (i + 1.0) / i);
    }

    const double HDLNucleationModel::
    CalculateNucleationRate(const Material *material,
                            const double temperature,
                            const double cooling_rate,
                            const double flow_energy,
                            const double flow_enthropy) const {
        const double boltzmann = 1.380650424e-23;
        const double gas_constant = 8.314472;
        const double Tm0 = material->GetTm0();

        double rDeltaG = 0.0;
        double rDeltaG_N = 0.0;
        if (temperature > material->GetTm0()) {
            rDeltaG = 0.0;
        } else if (temperature < 320.0) {
            rDeltaG = 1e12;
        } else {
            rDeltaG =
                    material->GetEnthalpy() * (1.0 - temperature / material->GetTm0());
            if (enthalpy_correction_) {
                rDeltaG *= 2 * temperature / (temperature + Tm0);
            }
        }
        rDeltaG += flow_energy;

        if (std::abs(n_ - 1.0) > 0.0001) {
            rDeltaG_N = pow(rDeltaG, n_);
        } else {
            rDeltaG_N = rDeltaG;
        }

        // thermal
        const double rAmp = C_ * boltzmann * temperature * rDeltaG;
        const double rArgExp1 = -U_ / (gas_constant * temperature);
        const double rArgExp2 = -K_ / (temperature * rDeltaG_N);
        const double rExp_1 = exp(rArgExp1);
        const double rExp_2 = exp(rArgExp2);
        const double thermal_nucleation = rAmp * rExp_1 * rExp_2;

        // athermal
        double athermal_nucleation = 0.0;
        if (athermal_) {
            const double rDeltaG_N2 = pow(rDeltaG, n_ + 2.0);
            double rDeltaG_q_dT = 0.0;
            if (enthalpy_correction_) {
                rDeltaG_q_dT = -2.0 * material->GetEnthalpy() *
                               (temperature * temperature + 2.0 * temperature * Tm0 - Tm0 * Tm0) /
                               (Tm0 * (Tm0 + temperature) * (Tm0 + temperature));
            } else {
                rDeltaG_q_dT = -material->GetEnthalpy() / Tm0;
            }
            const double rDeltaG_dT = rDeltaG_q_dT - flow_enthropy;

            if (fabs(rDeltaG_N2) > 1e-30) {
                athermal_nucleation = C_ath_ / rDeltaG_N2 * rExp_2;
                athermal_nucleation *= rDeltaG_dT * cooling_rate;
            }
        }
        if (thermal_nucleation + athermal_nucleation < 0.0) {
            return 0.0;
        }
        return thermal_nucleation + athermal_nucleation;
    }

}  // namespace SphaeroSim