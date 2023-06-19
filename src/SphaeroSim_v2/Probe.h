#ifndef _SPHAEROSIM_PROBE_H_
#define _SPHAEROSIM_PROBE_H_

#ifndef _INCLUDED_STDINT_H_

#include <stdint.h>

#define _INCLUDED_STDINT_H_
#endif

#ifndef _INCLUDED_MAP_H_

#include <map>

#define _INCLUDED_MAP_H_
#endif

#ifndef _SPHAEROSIM_EIGEN_LIBRARY_H_

#include "EigenLibrary.h"

#endif

namespace SphaeroSim {

    class Simulation;

    class Probe {
    public:
        enum Type {
            kCell,
            kElement
        };

        virtual ~Probe() {}

        virtual void
        NumericalValues(const Simulation &simulation,
                        const double point_in_time,
                        std::map<const std::string, double> *dst) const = 0;

        virtual void
        VectorValues(const Simulation &simulation,
                     const double point_in_time,
                     std::map<const std::string, Eigen::Vector3d> *dst) const = 0;

        // Getter
        inline const Type GetType() const {
            return type_;
        }

        inline const std::size_t GetElementID() const {
            return element_id_;
        }

        inline const std::string GetName() const {
            return name_;
        }

    protected:
        Probe(const Type type,
              const std::size_t element_id,
              const std::string &name);

    private:
        const Type type_;
        const std::size_t element_id_;
        const std::string name_;

        // NOT IMPLEMENTED
        Probe &operator=(const Probe &);
    };

    class ProbeCell : public Probe {
    public:
        ProbeCell(const std::size_t element_index,
                  const std::size_t cell_index,
                  const std::string &name);

        void NumericalValues(const Simulation &simulation,
                             const double point_in_time,
                             std::map<const std::string, double> *dst) const;

        void VectorValues(const Simulation &simulation,
                          const double point_in_time,
                          std::map<const std::string, Eigen::Vector3d> *dst) const;

    private:
        const std::size_t cell_index_;

        // NOT IMPLEMENTED
        ProbeCell &operator=(const ProbeCell &);
    };

    class ProbeElement : public Probe {
    public:
        ProbeElement(const std::size_t element_index,
                     const Eigen::Vector3s starting_cell,
                     const Eigen::Vector3s num_cells,
                     const std::string &name);

        void NumericalValues(const Simulation &simulation,
                             const double point_in_time,
                             std::map<const std::string, double> *dst) const;

        void VectorValues(const Simulation &simulation,
                          const double point_in_time,
                          std::map<const std::string, Eigen::Vector3d> *dst) const;

    private:
        const Eigen::Vector3s starting_cell_;
        const Eigen::Vector3s num_cells_;

        // NOT IMPLEMENTED
        ProbeElement &operator=(const ProbeElement &);
    };
}  // namespace SphaeroSim

#endif  // _SPHAEROSIM_PROBE_H_
