#ifndef _SPHAEROSIM_DISENGAGEMENT_TIME_H_
#define _SPHAEROSIM_DISENGAGEMENT_TIME_H_

#define _USE_MATH_DEFINES

#ifndef _INCLUDED_MATH_H_

#include <math.h>

#define _INCLUDED_MATH_H_
#endif

#undef _USE_MATH_DEFINES

namespace SphaeroSim {

    class DisengagementTime {
    public:
        virtual ~DisengagementTime() {}

        // interface
        virtual const double Calculate(const double temperature,
                                       const double density,
                                       const double mol_weight_ent) = 0;
    };

// constant disengagement time
    class ConstantDisengagementTime : public DisengagementTime {
    public:
        ConstantDisengagementTime(const double constant_value) :
                constant_value_(constant_value) {
        }

        // interface
        const double Calculate(const double temperature,
                               const double density,
                               const double mol_weight_ent) {
            temperature;
            density;
            mol_weight_ent;
            return constant_value_;
        }

    private:
        const double constant_value_;

        ConstantDisengagementTime();

        ConstantDisengagementTime &operator=(const ConstantDisengagementTime &);
    };

// disengagement time calculate with a WLF approach
    class WLFDisengagementTime : public DisengagementTime {
    public:
        WLFDisengagementTime(const double D1,
                             const double D2,
                             const double A1,
                             const double A2) :
                D1_(D1),
                D2_(D2),
                A1_(A1),
                A2_(A2) {
        }

        // interface
        const double Calculate(const double temperature,
                               const double density,
                               const double mol_weight_ent) {
            double delta_T = temperature - D2_;
            double val_1 = -A1_ * delta_T;
            double val_2 = A2_ + delta_T;
            double zero_shear = D1_ * exp(val_1 / val_2);

            double relaxation_time = 20.0 / (M_PI * M_PI) *
                                     mol_weight_ent / (8.314472 * temperature * density) *
                                     zero_shear;
            return relaxation_time;
        }

    private:
        const double D1_;
        const double D2_;
        const double A1_;
        const double A2_;

        WLFDisengagementTime();

        WLFDisengagementTime &operator=(const WLFDisengagementTime &);
    };

}  // namespace SphaeroSim

#endif  // _SPHAEROSIM_DISENGAGEMENT_TIME_H_