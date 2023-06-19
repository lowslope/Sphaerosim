#ifndef _SPHAEROSIM_NUCLEATION_MODEL_H_
#define _SPHAEROSIM_NUCLEATION_MODEL_H_

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

    class NucleationModel {
    public:
        virtual const double Rate(const double point_in_time,
                                  const std::size_t cell_index,
                                  const bool quiescent_only,
                                  const bool thermal_only,
                                  const Element *element) const = 0;

        virtual const double ScalingFactor() const = 0;

    protected:
        NucleationModel();

    private:
        // not implemented
        NucleationModel &operator=(const NucleationModel &);
    };

    class ConstantNucleationModel : public NucleationModel {
    public:
        ConstantNucleationModel(const double constant_rate);

        const double Rate(const double point_in_time,
                          const std::size_t cell_index,
                          const bool quiescent_only,
                          const bool thermal_only,
                          const Element *element) const;

        const double ScalingFactor() const override;

    private:
        const double constant_rate_;

        // not implemented
        ConstantNucleationModel &operator=(const ConstantNucleationModel &);
    };

    class HDLNucleationModel : public NucleationModel {
    public:
        HDLNucleationModel(const Material *material,
                           const double C,
                           const double U,
                           const double K,
                           const double n,
                           const bool enthalpy_correction,
                           const bool athermal,
                           const double C_ath);

        const double Rate(const double point_in_time,
                          const std::size_t cell_index,
                          const bool quiescent_only,
                          const bool thermal_only,
                          const Element *element) const;

        const double ScalingFactor() const override;

        const double AthermalConstant(const double volume) const;

    private:
        const double C_;
        const double U_;
        const double K_;
        const double n_;
        const double C_ath_;
        const bool athermal_;
        const bool enthalpy_correction_;

        double scaling_factor_;

        const double CalculateNucleationRate(const Material *material,
                                             const double temperature,
                                             const double cooling_rate,
                                             const double flow_energy,
                                             const double flow_enthropy) const;

        // not implemented
        HDLNucleationModel &operator=(const HDLNucleationModel &);
    };

}  // namespace SphaeroSim

#endif  // _SPHAEROSIM_NUCLEATION_MODEL_H_