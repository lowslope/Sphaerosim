#ifndef _SPHAEROSIM_SIMULATION_H_
#define _SPHAEROSIM_SIMULATION_H_

#ifndef _INCLUDED_STRING_H_

#include <string>

#define _INCLUDED_STRING_H_
#endif

#ifndef _INCLUDED_VECTOR_H_

#include <vector>

#define _INCLUDED_VECTOR_H_
#endif

#ifndef _INCLUDED_MAP_H_

#include <map>

#define _INCLUDED_MAP_H_
#endif

#ifndef _INCLUDED_MUTEX_H_

#include <mutex>

#define _INCLUDED_MUTEX_H_
#endif

#ifndef _INCLUDED_CTIME_H_

#include <ctime>

#define _INCLUDED_CTIME_H_
#endif

#ifndef _SPHAEROSIM_EIGEN_LIBRARY_H_

#include "EigenLibrary.h"

#endif

#ifndef _SPHAEROSIM_INTERPOLATION_H_

#include "Interpolation.h"

#endif

#ifndef _SPHAEROSIM_ELEMENT_H_

#include "Element.h"

#endif

#ifndef _SPHAEROSIM_MATERIAL_H_

#include "Material.h"

#endif

namespace SphaeroSim {

    class Probe;

    class RandomNumberGenerator;

    class SpheruliteInfo;

    class GrowthModel;

    class NucleationModel;

    class Simulation {
    public:
        class NumericalOptions {
        public:
            NumericalOptions(const int32_t cache_number_entries,
                             const double threshold_pure_flow,
                             const double threshold_incompressibility,
                             const double timeresolution_cache,
                             const double timeresolution_fdm,
                             const double timestep_fdm,
                             const bool consider_pure_shear,
                             const bool consider_pure_elongational) :
                    cache_number_entries_(cache_number_entries),
                    threshold_pure_flow_(threshold_pure_flow),
                    threshold_incompressibility_(threshold_incompressibility),
                    timeresolution_cache_(timeresolution_cache),
                    timeresolution_fdm_(timeresolution_fdm),
                    timestep_fdm_(timestep_fdm),
                    consider_pure_shear_(consider_pure_shear),
                    consider_pure_elongational_(consider_pure_elongational) {
            }

            NumericalOptions(const NumericalOptions &src) :
                    cache_number_entries_(src.cache_number_entries()),
                    threshold_pure_flow_(src.threshold_pure_flow()),
                    threshold_incompressibility_(src.threshold_incompressibility()),
                    timeresolution_cache_(src.timeresolution_cache()),
                    timeresolution_fdm_(src.timeresolution_fdm()),
                    timestep_fdm_(src.timestep_fdm()),
                    consider_pure_shear_(src.consider_pure_shear()),
                    consider_pure_elongational_(src.consider_pure_elongational()) {
            }

            // getter
            inline double threshold_pure_flow() const {
                return threshold_pure_flow_;
            }

            inline double threshold_incompressibility() const {
                return threshold_incompressibility_;
            }

            inline double timeresolution_cache() const {
                return timeresolution_cache_;
            }

            inline double timeresolution_fdm() const {
                return timeresolution_fdm_;
            }

            inline double timestep_fdm() const {
                return timestep_fdm_;
            }

            inline bool consider_pure_shear() const {
                return consider_pure_shear_;
            }

            inline bool consider_pure_elongational() const {
                return consider_pure_elongational_;
            }

            inline int32_t cache_number_entries() const {
                return cache_number_entries_;
            }

        private:
            const int32_t cache_number_entries_;
            const double threshold_pure_flow_;
            const double threshold_incompressibility_;
            const double timeresolution_cache_;
            const double timeresolution_fdm_;
            const double timestep_fdm_;
            const bool consider_pure_shear_;
            const bool consider_pure_elongational_;

            NumericalOptions &operator=(const NumericalOptions &);
        };

        class VTKResultOption {
        public:
            VTKResultOption(const bool save_vtk_files,
                            const std::string &subfolder,
                            const bool all_result_points,
                            const bool absolute_coordinates,
                            const bool temperature,
                            const bool aspectratio,
                            const bool crystallizationtime,
                            const bool phasestate,
                            const bool spheruliteid,
                            const bool spherulitediameter,
                            const bool randomcolor,
                            const bool enthropychange,
                            const bool induction,
                            const bool coolingrate,
                            const bool location,
                            const bool velocity,
                    //new checkboundarylayer option
                            const bool checkboundarylayer) :
                    _subfolder(subfolder),
                    _all_result_points(all_result_points),
                    _save_vtk_files(save_vtk_files),
                    _absolute_coordinates(absolute_coordinates),
                    temperature_(temperature),
                    aspectratio_(aspectratio),
                    crystallizationtime_(crystallizationtime),
                    phasestate_(phasestate),
                    spheruliteid_(spheruliteid),
                    spherulitediameter_(spherulitediameter),
                    randomcolor_(randomcolor),
                    enthropychange_(enthropychange),
                    induction_(induction),
                    coolingrate_(coolingrate),
                    location_(location),
                    velocity_(velocity),
                    //new checkboundarylayer option
                    checkboundarylayer_(checkboundarylayer) {
            }

            inline const bool GetAllResultPoints() const {
                return _all_result_points;
            }

            inline const bool GetSaveVTKFiles() const {
                return _save_vtk_files;
            }

            inline const std::string GetSubfolder() const {
                return _subfolder;
            }

            inline const bool GetAbsoluteCoordinates() const {
                return _absolute_coordinates;
            }

            inline const bool GetTemperature() const {
                return temperature_;
            }

            inline const bool GetAspectratio() const {
                return aspectratio_;
            }

            inline const bool GetCrystallizationtime() const {
                return crystallizationtime_;
            }

            inline const bool GetPhasestate() const {
                return phasestate_;
            }

            inline const bool GetSpheruliteid() const {
                return spheruliteid_;
            }

            inline const bool GetSpherulitediameter() const {
                return spherulitediameter_;
            }

            inline const bool GetRandomcolor() const {
                return randomcolor_;
            }

            inline const bool GetEnthropychange() const {
                return enthropychange_;
            }

            inline const bool GetInduction() const {
                return induction_;
            }

            inline const bool GetCoolingrate() const {
                return coolingrate_;
            }

            inline const bool GetLocation() const {
                return location_;
            }

            inline const bool GetVelocity() const {
                return velocity_;
            }

            //new checkboundarylayer option
            inline const bool Getcheckboundarylayer() const {
                return checkboundarylayer_;
            }


            const VTKResultOption operator=(const VTKResultOption &src) {
                return VTKResultOption(src.GetSaveVTKFiles(),
                                       src.GetSubfolder(),
                                       src.GetAllResultPoints(),
                                       src.GetAbsoluteCoordinates(),
                                       src.GetTemperature(),
                                       src.GetAspectratio(),
                                       src.GetCrystallizationtime(),
                                       src.GetPhasestate(),
                                       src.GetSpheruliteid(),
                                       src.GetSpherulitediameter(),
                                       src.GetRandomcolor(),
                                       src.GetEnthropychange(),
                                       src.GetInduction(),
                                       src.GetCoolingrate(),
                                       src.GetLocation(),
                                       src.GetVelocity(),
                        //new checkboundarylayer option
                                       src.Getcheckboundarylayer());
            }

        private:
            const std::string _subfolder;
            const bool _all_result_points;
            const bool _save_vtk_files;
            const bool _absolute_coordinates;
            const bool temperature_;
            const bool aspectratio_;
            const bool crystallizationtime_;
            const bool phasestate_;
            const bool spheruliteid_;
            const bool spherulitediameter_;
            const bool randomcolor_;
            const bool enthropychange_;
            const bool induction_;
            const bool coolingrate_;
            const bool location_;
            const bool velocity_;
            //new checkboundarylayer
            const bool checkboundarylayer_;
        };

        class ResultStep {
        public:
            enum StepType {
                kTime,
                kCrystallizationDegree
            };

            ResultStep(const StepType type,
                       const double stepsize) :
                    type_(type),
                    stepsize_(stepsize) {
            }

            ResultStep(const ResultStep &src) :
                    type_(src.type_),
                    stepsize_(src.stepsize_) {
            }

            // getter
            inline const StepType GetType() const {
                return type_;
            }

            inline const double GetStepsize() const {
                return stepsize_;
            }

        private:
            const StepType type_;
            const double stepsize_;

            ResultStep &operator=(const ResultStep &);
        };

        Simulation(const std::string &name,
                   const ResultStep convergence_criterium,
                   const double cell_size_m,
                   const InterpolationType temporal_interpolation_type,
                   const InterpolationType spatial_interpolation_type,
                   const std::string &output_folder,
                   const bool create_new_subfolder,
                   const VTKResultOption &vtk_result_option,
                   const char csv_separator,
                   const std::vector<bool> &noslip_neighbourhood,
                   const std::string &calibration_file,
                   const std::vector<uint32_t> &rng_state_registers,
                   const double integral_error,
                   const bool save_memory,
                   const bool calculate_flow_energy,
                   const double cache_convert_temperature,
                   const double cache_convert_velocity,
                   const ResultStep &result_step,
                   const ResultStep &statistics_step,
                   const NumericalOptions &numerical_options,
                   GrowthModel *growth_model,
                   NucleationModel *nucleation_model,
                   Material *material);

        virtual ~Simulation();

        // interface
        virtual void BuildElements(std::vector<Element *> *dst,
                                   RandomNumberGenerator *rng) = 0;

        virtual void
        SetElementNeighbourhood(std::map<const uint64_t, Element *> *dst) const = 0;

        virtual const double Timestep() = 0;

        void Run();

        void ListOfElementIDs(std::vector<uint64_t> *dst) const;

        inline void AddProbeDescription(const std::string &type,
                                        const std::string &name,
                                        const uint64_t element_id,
                                        const std::string &location_type,
                                        const Eigen::Vector3d location,
                                        const Eigen::Vector3d size) {
            ProbeDescription desc;
            desc.type_ = type;
            desc.element_id_ = element_id;
            desc.location_type_ = location_type;
            desc.location_ = location;
            desc.size_ = size;
            desc.name_ = name;
            probe_description_.push_back(desc);
        }

        inline void
        AddPredefinedNucleiList(const uint64_t element_id,
                                const std::string &location_type,
                                const std::vector<Eigen::Vector3d> &locations) {
            predefined_nuclei_.push_back(new PredefinedNucleiList(locations,
                                                                  location_type,
                                                                  element_id));
        }

        inline const std::map<const uint64_t, SpheruliteInfo *> *const getSpherulite_infos() {
            return &spherulite_infos_;
        }

        inline void
        AddPredefinedNucleiDistribution(const uint64_t element_id,
                                        const double num_nuclei_per_volume) {
            predefined_nuclei_.push_back(
                    new PredefinedNucleiDistribution(num_nuclei_per_volume,
                                                     element_id));
        }

        // add a new spherulite to the pool
        const uint64_t NewSpherulite(const std::size_t start_location,
                                     const double creation_time,
                                     const Element *owner);

        // retrieve the spherulite center
        const complex_index SpheruliteCenter(const uint64_t spherulite_id);

        // add a cell to the spherulite
        void AddCellToSpherulite(const uint64_t spherulite_id,
                                 const complex_index &cell_to_add,
                                 const Eigen::Vector3i &periodic_continuation,
                                 const bool spherulitesAdded = true);

        // get the diameter of a spherulite
        const double SpheruliteDiameter(const uint64_t &spherulite_id);
        const double SpheruliteDiameter(const SpheruliteInfo *const spherulite);

        const double SpheruliteDiameterNoDataracePossible(const SpheruliteInfo *const spherulite);

        // retrieve the random color for a spherulite
        const int32_t SpheruliteRandomColor(const uint64_t &spherulite_id);

        const int32_t SpheruliteRandomColor(const SpheruliteInfo *const spherulite);

        // retrieve the aspect ratio of a spherulite
        const double SpheruliteAspectRatio(const uint64_t &spherulite_id);

        const double SpheruliteAspectRatio(const SpheruliteInfo *const spherulite);

        const double SpheruliteAspectRatioNoDatraracePossible(const SpheruliteInfo *const spherulite);

        // get the average diameter and standard deviation of the spherulite diameter
        // for all spherulites which have their origin in a specific element
        // return: <average diameter, standard deviation>
        const std::pair<const double, const double>
        AverageSpheruliteDiameterForElement(const Element *element,
                                            const Eigen::Vector3s start_location,
                                            const Eigen::Vector3s num_cells,
                                            std::size_t *dst_num_spherulites);

        // calculate the relative crystallinity of the entire simulation area
        const double RelativeCrystallinity() const;

        // return the total number of spherulites
        const uint64_t TotalNumberSpherulites() const;

        // Getter
        inline const std::string GetName() const {
            return name_;
        }

        inline const std::size_t GetNumElements() const {
            return elements_.size();
        }

        inline const double GetCellSizeM() const {
            return cell_size_m_;
        }

        inline const std::size_t GetCurrentRunIndex() const {
            return current_run_;
        }

        inline const InterpolationType GetTemporalInterpolationType() const {
            return temporal_interpolation_type_;
        }

        inline const InterpolationType GetSpatialInterpolationType() const {
            return spatial_interpolation_type_;
        }

        inline const std::size_t GetNumRuns() const {
            return num_runs_;
        }

        inline const std::string GetOutputFolder() const {
            return output_folder_;
        }

        inline const VTKResultOption *GetVTKOptions() const {
            return &vtk_result_option_;
        }

        inline const Element *GetElement(const uint64_t id) const {
            if (elements_.find(id) == elements_.end()) {
                return nullptr;
            }
            return elements_.at(id);
        }

        inline Element *GetElement(const uint64_t id) {
            if (elements_.find(id) == elements_.end()) {
                return nullptr;
            }
            return elements_.at(id);
        }

        inline const Probe *GetProbe(const std::size_t index) const {
            return probes_[index];
        }

        inline const std::size_t GetNumProbes() const {
            return probes_.size();
        }

        inline const char GetCSVSeparator() const {
            return csv_separator_;
        }

        // index order: (+/- X_AXIS; +/- Y_AXIS; +/- Z_AXIS)
        inline const bool GetNoslipNeighbourhood(const std::size_t index) const {
            return noslip_neighbourhood_[index];
        }

        inline const GrowthModel *GetGrowthModel() const {
            return growth_model_;
        }

        inline const NucleationModel *GetNucleationModel() const {
            return nucleation_model_;
        }

        inline const double GetIntegralError() const {
            return integral_error_;
        }

        inline const bool GetSaveMemory() const {
            return save_memory_;
        }

        inline const Material *GetMaterial() const {
            return material_;
        }

        inline const bool GetCalculateFlowEnergy() const {
            return calculate_flow_energy_;
        }

        inline const double GetCacheConvertTemperature() const {
            return cache_convert_temperature_;
        }

        inline const double GetCacheConvertVelocity() const {
            return cache_convert_velocity_;
        }

        inline const ResultStep *GetResultStep() const {
            return &result_step_;
        }

        inline const ResultStep *GetStatisticsStep() const {
            return &statistics_step_;
        }

        inline const NumericalOptions *GetNumericalOptions() const {
            return &numerical_options_;
        }

    protected:
        inline void SetNumRuns(const std::size_t num_runs) {
            num_runs_ = num_runs;
        }

    private:
        class PredefinedNuclei {
        public:
            PredefinedNuclei(const uint64_t element_id) :
                    element_id_(element_id) {
            }

            virtual void SpreadInElement(const Simulation &simulation) = 0;

            // getter
            inline const uint64_t GetElementID() const {
                return element_id_;
            }

        private:
            const uint64_t element_id_;

            // not implemente
            PredefinedNuclei &operator=(const PredefinedNuclei &);
        };

        class PredefinedNucleiList : public PredefinedNuclei {
        public:
            PredefinedNucleiList(const std::vector<Eigen::Vector3d> &locations,
                                 const std::string &location_type,
                                 const uint64_t element_id) :
                    PredefinedNuclei(element_id),
                    _locations(locations),
                    location_type_(location_type) {
            }

            void SpreadInElement(const Simulation &simulation) {
                std::vector<uint64_t> element_ids;
                if (GetElementID() == std::numeric_limits<uint64_t>::max()) {
                    simulation.ListOfElementIDs(&element_ids);
                } else {
                    if (simulation.GetElement(GetElementID()) == nullptr) {
                        throw Exception("PredefinedNucleiList",
                                        "Element-ID does not exist.");
                    }
                    element_ids.push_back(GetElementID());
                }

                double start = omp_get_wtime();
                //#pragma omp parallel
                {
                    //#pragma omp for nowait
                    for (auto &element_id : element_ids) {
                        Element *el = simulation.elements_.at(element_id);
                        for (auto &_location : _locations) {
                            Eigen::Vector3s pos = simulation.IndexCoordinate(location_type_,
                                                                             el,
                                                                             _location);
                            el->AddNucleus(pos, 0.0);
                        }
                    }
                    printf("T%i: %fs\n", omp_get_thread_num(), start - omp_get_wtime());
                    //#pragma omp barrier
                };
            }

        private:
            const std::vector<Eigen::Vector3d> _locations;
            const std::string location_type_;

            // not implemented
            PredefinedNucleiList &operator=(const PredefinedNucleiList &);
        };

        class PredefinedNucleiDistribution : public PredefinedNuclei {
        public:
            PredefinedNucleiDistribution(const double number_per_volume,
                                         const uint64_t element_id) :
                    PredefinedNuclei(element_id),
                    _number_per_volume(number_per_volume) {
            }

            void SpreadInElement(const Simulation &simulation) {
                printf("SpreadInElement: 1\n");
                std::vector<uint64_t> element_ids;
                if (GetElementID() == std::numeric_limits<uint64_t>::max()) {
                    printf("SpreadInElement: 2.1\n");
                    simulation.ListOfElementIDs(&element_ids);
                } else {
                    printf("SpreadInElement: 2.2\n");
                    if (simulation.GetElement(GetElementID()) == nullptr) {
                        throw Exception("PredefinedNucleiList",
                                        "Element-ID does not exist.");
                    }
                    element_ids.push_back(GetElementID());
                }

                printf("SpreadInElement: 3\n");
                for (auto &element_id : element_ids) {
                    Element *el = simulation.elements_.at(element_id);
                    const auto num_nuclei = static_cast<uint64_t>(
                            _number_per_volume * el->GetVolume() + 0.5);
                    for (std::size_t cur = 0; cur < num_nuclei; ++cur)
                        el->AddRandomNucleus(0.0);
                }
                printf("SpreadInElement: 4\n");
            }

        private:
            const double _number_per_volume;

            // not implemented
            PredefinedNucleiList &operator=(const PredefinedNucleiList &);
        };

        const double UpdateElement(const double current_time,
                                   Element *element);

        const std::string
        CreateOutputFolder(const bool create_subfolder,
                           const std::size_t current_run_index) const;

        void UpdateStatisticsFile(const int32_t step_index,
                                  const double current_time) const;

        void UpdateSpheruliteFile() const;

        void WriteSpheruliteCenterFile();

        void CopyCalibrationFile() const;

        void CreateInitialEvents();

        void ClearElements();

        void CreateProbes();

        void SaveRNGRegisters(const std::string &filename,
                              const RandomNumberGenerator *rng) const;

        bool ResultStepFullfilled(const ResultStep &result_step,
                                  const double current_time,
                                  double *last_evaluation_point) const;

        inline
        Eigen::Vector3s IndexCoordinate(const std::string &type,
                                        const Element *element,
                                        const Eigen::Vector3d &position) const {
            Eigen::Vector3s result;
            if (type == "Relative") {
                for (std::size_t cur_dim = 0; cur_dim < 3; ++cur_dim)
                    result[cur_dim] =
                            static_cast<std::size_t>(element->GetNumCells()[cur_dim] *
                                                     position[cur_dim]);
            } else if (type == "Absolute") {
                result =
                        element->CellIndicesFromPosition(position, true);
            } else if (type == "Index") {
                result = Eigen::Vector3s(
                        static_cast<std::size_t>(position[0]),
                        static_cast<std::size_t>(position[1]),
                        static_cast<std::size_t>(position[2]));
            } else {
                throw Exception("Invalid coordinate type.",
                                "Allowed values are 'Relative', 'Absolute' and 'Index'");
            }
            return result;
        }

        const std::string name_;
        const ResultStep convergence_criterium_;
        const double cell_size_m_;
        const InterpolationType temporal_interpolation_type_;
        const InterpolationType spatial_interpolation_type_;
        const bool create_subfolder_;
        const std::string root_folder_;
        const std::string calibration_file_;
        const VTKResultOption vtk_result_option_;
        const char csv_separator_;
        const double integral_error_;
        const bool save_memory_;
        const double cache_convert_temperature_;
        const bool calculate_flow_energy_;
        const double cache_convert_velocity_;
        const ResultStep result_step_;
        const ResultStep statistics_step_;

        NumericalOptions numerical_options_;
        Material *material_;
        GrowthModel *growth_model_;
        NucleationModel *nucleation_model_;

        std::map<const uint64_t, Element *> elements_;
        std::vector<Probe *> probes_;
        std::vector<bool> noslip_neighbourhood_;
        std::vector<uint32_t> rng_state_registers_;
        std::map<const uint64_t, SpheruliteInfo *> spherulite_infos_;
        std::mutex spherulite_mutex_;

        struct ProbeDescription {
            std::string type_;
            std::string name_;
            uint64_t element_id_;
            std::string location_type_;
            Eigen::Vector3d location_;
            Eigen::Vector3d size_;
        };
        std::vector<ProbeDescription> probe_description_;
        std::vector<PredefinedNuclei *> predefined_nuclei_;

        std::size_t current_run_;
        std::size_t num_runs_;
        std::string output_folder_;

        clock_t last_clock_value_;

        // NOT IMPLEMENTED
        Simulation &operator=(const Simulation &);
    };

}  // namespace SphaeroSim
#endif  // _SPHAEROSIM_SIMULATION_H_
