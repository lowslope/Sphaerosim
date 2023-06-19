#include "Probe.h"

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

#ifndef _SPHAEROSIM_SIMULATION_H_

#include "Simulation.h"

#endif

#ifndef _SPHAEROSIM_EXCEPTION_H_
#include "Exception.h"
#endif

#ifndef _SPHAEROSIM_ELEMENT_H_
#include "Element.h"
#endif

namespace SphaeroSim {

    Probe::Probe(const Type type,
                 const std::size_t element_id,
                 const std::string &name) :
            type_(type),
            element_id_(element_id),
            name_(name) {
    }

    ProbeCell::ProbeCell(const std::size_t element_index,
                         const std::size_t cell_index,
                         const std::string &name) :
            Probe(kCell, element_index, name),
            cell_index_(cell_index) {
    }

    void ProbeCell::
    NumericalValues(const Simulation &simulation,
                    const double point_in_time,
                    std::map<const std::string, double> *dst) const {
        std::stringstream ss;
        if (simulation.GetElement(GetElementID()) == NULL) {
            ss << "The element '" << GetElementID() << "' does not exist.";
            throw Exception("Invalid probe index", ss.str());
        }

        simulation.GetElement(GetElementID())->
                GetNumericalFieldsForCell(cell_index_,
                                          point_in_time,
                                          dst);
    }

    void ProbeCell::
    VectorValues(const Simulation &simulation,
                 const double point_in_time,
                 std::map<const std::string, Eigen::Vector3d> *dst) const {
        std::stringstream ss;
        if (simulation.GetElement(GetElementID()) == NULL) {
            ss << "The element '" << GetElementID() << "' does not exist.";
            throw Exception("Invalid probe index", ss.str());
        }

        simulation.GetElement(GetElementID())->
                GetVectorFieldsForCell(cell_index_,
                                       point_in_time,
                                       dst);
    }

    ProbeElement::ProbeElement(const std::size_t element_index,
                               const Eigen::Vector3s starting_cell,
                               const Eigen::Vector3s num_cells,
                               const std::string &name) :
            Probe(kElement, element_index, name),
            starting_cell_(starting_cell),
            num_cells_(num_cells) {
    }

    void ProbeElement::
    NumericalValues(const Simulation &simulation,
                    const double point_in_time,
                    std::map<const std::string, double> *dst) const {
        std::stringstream ss;
        if (simulation.GetElement(GetElementID()) == NULL) {
            ss << "The element '" << GetElementID() << "' does not exist.";
            throw Exception("Invalid probe index", ss.str());
        }
        simulation.GetElement(GetElementID())->
                GetNumericalFieldsForElement(point_in_time,
                                             starting_cell_,
                                             num_cells_,
                                             dst);
    }

    void ProbeElement::
    VectorValues(const Simulation &simulation,
                 const double point_in_time,
                 std::map<const std::string, Eigen::Vector3d> *dst) const {
        std::stringstream ss;
        if (simulation.GetElement(GetElementID()) == NULL) {
            ss << "The element '" << GetElementID() << "' does not exist.";
            throw Exception("Invalid probe index", ss.str());
        }

        simulation.GetElement(GetElementID())->
                GetVectorFieldsForElement(point_in_time,
                                          starting_cell_,
                                          num_cells_,
                                          dst);
    }
}  // namespace SphaeroSim
