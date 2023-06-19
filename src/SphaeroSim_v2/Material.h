#ifndef _SPHAEROSIM_MATERIAL_H_
#define _SPHAEROSIM_MATERIAL_H_

#ifndef _INCLUDED_STRING_H_
#include <string>
#define _INCLUDED_STRING_H_
#endif

#ifndef _SPHAEROSIM_DISENGAGEMENT_TIME_H_

#include "DisengagementTime.h"

#endif

namespace SphaeroSim {

    class Material {
    public:
        Material(const std::string &name,
                 const double density,
                 const double T_m0,
                 const double T_inf,
                 const double mol_weight_ent,
                 const double enthalpy,
                 DisengagementTime *disengagement_time) :
                name_(name),
                density_(density),
                T_m0_(T_m0),
                T_inf_(T_inf),
                mol_weight_ent_(mol_weight_ent),
                enthalpy_(enthalpy) {
            disengagement_time_ = disengagement_time;
        }

        virtual ~Material() {
            if (disengagement_time_ != NULL) {
                delete disengagement_time_;
            }
        }

        // getter
        inline const std::string GetName() const {
            return name_;
        }

        inline const double GetDensity() const {
            return density_;
        }

        inline const double GetTm0() const {
            return T_m0_;
        }

        inline const double GetTinf() const {
            return T_inf_;
        }

        inline const double GetEnthalpy() const {
            return enthalpy_;
        }

        inline const double GetMolWeightEnt() const {
            return mol_weight_ent_;
        }

        inline const double GetDisengagementTime(const double temperature) const {
            return disengagement_time_->Calculate(temperature,
                                                  GetDensity(),
                                                  GetMolWeightEnt() / 1000.0);
        }

        inline const double GetDisengagementTimeDT(const double temperature) const {
            temperature;
            return 0.0;
        }

    private:
        const std::string name_;
        const double density_;
        const double T_m0_;
        const double T_inf_;
        const double enthalpy_;
        const double mol_weight_ent_;
        DisengagementTime *disengagement_time_;

        Material();

        Material &operator=(const Material &);
    };

}

#endif  // _SPHAEROSIM_MATERIAL_H_
