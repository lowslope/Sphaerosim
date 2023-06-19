#ifndef _SPHAEROSIM_STATE_EVENT_H_
#define _SPHAEROSIM_STATE_EVENT_H_

#ifndef _INCLUDED_STDINT_H_

#include <stdint.h>

#define _INCLUDED_STDINT_H_
#endif

#ifndef _INCLUDED_VECTOR_H_

#include <vector>

#define _INCLUDED_VECTOR_H_
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

    class Probe;

    class Element;

    typedef std::pair<const Element *, std::size_t> complex_index;

    class StateEvent {
    public:
        enum Type {
            kNucleation,
            kGrowth,
            kWriteResults,
            kNullEvent
        };

        virtual ~StateEvent() {}


        virtual void Execute() = 0;

        virtual bool SameEvent(const StateEvent *compare_to) const = 0;

        // Getter
        inline const double GetPointInTime() const {
            return point_in_time_;
        }

        inline const double GetOffsetNewEvent() const {
            return offset_new_event_;
        }

        inline Element *GetElement() {
            return element_;
        }

        inline const Element *GetElement() const {
            return element_;
        }

        inline const Type GetType() const {
            return type_;
        }

    protected:
        StateEvent(const Type type,
                   const double point_in_time,
                   const double offset_new_event,
                   Element *element);

    private:
        const Type type_;
        const double point_in_time_;
        const double offset_new_event_;
        Element *element_;

        // NOT IMPLEMENTED
        StateEvent &operator=(const StateEvent &);
    };

// Dummy event
    class NullEvent : public StateEvent {
    public:
        NullEvent(const double point_in_time,
                  Element *element);

        void Execute();

        bool SameEvent(const StateEvent *compare_to) const;

    private:
        // NOT IMPLEMENTED
        NullEvent &operator=(const NullEvent &);
    };

// nucleation event for an element
    class NucleationEvent : public StateEvent {
    public:
        NucleationEvent(const double point_in_time,
                        const double offset_new_event,
                        const double last_nucleation,
                        Element *element);

        void Execute();

        bool SameEvent(const StateEvent *compare_to) const;

    private:
        const double last_nucleation_;

        // NOT IMPLEMENTED
        NucleationEvent &operator=(const NucleationEvent &);
    };

// nucleation event for an element or a cell of an element
    class GrowthEvent : public StateEvent {
    public:
        enum Type {
            kMonteCarlo,
            kRayTracing,
            kCCG
        };

        GrowthEvent(const std::size_t cell_index,
                    const complex_index creating_cell,
                    const double point_in_time,
                    const double offset_new_event,
                    const Type growth_type,
                    const Eigen::Vector3i periodic_continuation,
                    Element *element);

        void Execute();

        bool SameEvent(const StateEvent *compare_to) const;

    private:
        void
        IncreaseNeighbourCCG(const complex_index cur_cell,
                             const double value,
                             std::map<complex_index, double> *crystal_state_increase,
                             std::map<complex_index, double> *partial_increase);

        const std::size_t cell_index_;
        const complex_index creating_cell_;
        const Type growth_type_;
        const Eigen::Vector3i periodic_continuation_;

        // NOT IMPLEMENTED
        GrowthEvent &operator=(const GrowthEvent &);
    };

// event for the printout of the current results
    class WriteResultEvent : public StateEvent {
    public:
        WriteResultEvent(const double point_in_time,
                         const double offset_new_event,
                         const Simulation &simulation,
                         Element *element);

        void Execute();

        bool SameEvent(const StateEvent *compare_to) const;

    private:
        void WriteVtkFile(const Simulation &simulation) const;

        void WriteProbeCSVString(const Probe *probe,
                                 const Simulation &simulation,
                                 std::string *header,
                                 std::string *values) const;

        void WriteProbeCSVFile(const std::string &filename,
                               const std::string &header,
                               const std::string &values) const;

        const Simulation *simulation_;

        // NOT IMPLEMENTED
        WriteResultEvent &operator=(const WriteResultEvent &);
    };

// comparator for the comparison of stateevents
    class StateEventCompare {
    public:
        bool operator()(const StateEvent *event_1,
                        const StateEvent *event_2) const {
            return event_1->GetPointInTime() < event_2->GetPointInTime();
        }
    };

}  // namespace SphaeroSim

#endif  // _SPHAEROSIM_STATE_EVENT_H_
