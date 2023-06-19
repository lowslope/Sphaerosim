#ifndef _SPHAEROSIM_APPLICATION_H_
#define _SPHAEROSIM_APPLICATION_H_

#ifndef _INCLUDED_STRING_H_

#include <string>

#define _INCLUDED_STRING_H_
#endif

#ifndef _INCLUDED_SSTREAM_H_

#include <sstream>

#define _INCLUDED_SSTREAM_H_
#endif

#ifndef _INCLUDED_VECTOR_H_

#include <vector>

#define _INCLUDED_VECTOR_H_
#endif

#ifndef _SPHAEROSIM_EIGEN_LIBRARY_H_

#include "EigenLibrary.h"

#endif

#ifndef _SPHAEROSIM_EXCEPTION_H_

#include "Exception.h"

#endif

#ifndef _SPHAEROSIM_INTERPOLATION_H_

#include "Interpolation.h"

#endif

#ifndef _SPHAEROSIM_SIMULATION_H_

#include "Simulation.h"

#endif

namespace tinyxml2 {
    class XMLElement;

    class XMLDocument;
}

namespace SphaeroSim {

    class Probe;

    class GrowthModel;

    class Application {
    public:
        explicit Application(const std::string &parameter_file);

        virtual ~Application();

        void Run();

    private:
        const InterpolationType
        ParseInterpolationType(const char *value,
                               const tinyxml2::XMLElement &cur_element) const;

        void ParseParameterFile(const std::string &filename);

        void ParseSimulationRun(const tinyxml2::XMLElement &cur_element);

        Simulation *CreateSimulationObject(const std::string &type,
                                           const tinyxml2::XMLElement &cur_element,
                                           const std::string &name,
                                           const Simulation::ResultStep &convergence,
                                           const double cell_size_m,
                                           const InterpolationType temporal_interpolation,
                                           const InterpolationType spatial_interpolation,
                                           const std::string output_folder,
                                           const bool create_new_subfolder,
                                           const Simulation::VTKResultOption &vtk_options,
                                           const char csv_separator,
                                           const double integral_error,
                                           const bool save_memory,
                                           const bool calculate_flow_energy,
                                           const double cache_convert_temperature,
                                           const double cache_convert_velocity,
                                           const Simulation::ResultStep &result_step,
                                           const Simulation::ResultStep &statistics_step,
                                           const Simulation::NumericalOptions &numerical_options,
                                           Material *material) const;

        void ParseProbeDescription(const tinyxml2::XMLElement &cur_element,
                                   Simulation *simulation) const;

        Simulation::NumericalOptions
        ParseNumericalOptions(const tinyxml2::XMLElement &cur_element) const;

        void ParsePredefinedNuclei(const tinyxml2::XMLElement &cur_element,
                                   Simulation *simulation) const;

        GrowthModel *ParseGrowthModel(const Material *material,
                                      const tinyxml2::XMLElement &cur_element) const;

        Material *ParseMaterial(const tinyxml2::XMLElement &cur_element) const;

        DisengagementTime *
        ParseDisengagementTime(const tinyxml2::XMLElement &cur_element) const;

        Simulation::ResultStep
        ParseResultStep(const tinyxml2::XMLElement &cur_element) const;

        Simulation::VTKResultOption
        ParseVTKOptions(const tinyxml2::XMLElement &cur_element);

        NucleationModel *
        ParseNucleationModel(const Material *material,
                             const tinyxml2::XMLElement &cur_element) const;

        std::string XMLErrorText(const tinyxml2::XMLDocument *xml_file);

        bool XMLValueExists(const char *value,
                            const tinyxml2::XMLElement &element) const;

        std::string XMLValueString(const char *value,
                                   const tinyxml2::XMLElement &element) const;

        template<class T>
        const T
        XMLToValue(const char *value, const tinyxml2::XMLElement &element) const {
            return StringToValue<T>(XMLValueString(value, element));
        }

        template<class T>
        std::vector<T>
        XMLToArray(const char *value, const tinyxml2::XMLElement &element) const {
            if (XMLValueExists(value, element) == false) {
                return std::vector<T>();
            }  // return an empty vector
            return StringToList<T>(XMLValueString(value, element), ';');
        }

        // Convert a string to a value (generic types)
        template<class T>
        const T
        StringToValue(const std::string &str_value) const {
            T result_value;
            std::stringstream stream(str_value);
            stream >> result_value;
            if (stream.fail()) {
                throw Exception("Value conversion failed.",
                                str_value);
            }
            return result_value;
        }

        // Convert a string to an Eigen::Vector<dimension>
        template<class type, std::size_t dimension>
        const Eigen::Matrix<type, dimension, 1, 0, dimension, 1>
        StringToVector(const std::string &str_value) const {
            std::size_t cur_pos = str_value.find('(');
            std::size_t end_pos = str_value.find(')');
            if (cur_pos == std::string::npos ||
                end_pos == std::string::npos ||
                end_pos <= cur_pos + 1) {
                throw Exception("Value conversion to Eigen:Vector failed.",
                                str_value);
            }

            std::string value_str =
                    str_value.substr(cur_pos + 1, end_pos - cur_pos - 1);
            std::vector<type> values = StringToList<type>(value_str, ',');
            if (values.size() != dimension) {
                throw Exception("Value conversion to Eigen:Vector failed.",
                                "Wrong number of values specified for the " \
                      "vector dimension.");
            }

            Eigen::Matrix<type, dimension, 1, 0, dimension, 1> result;
            for (std::size_t cur = 0; cur < dimension; ++cur)
                result[cur] = values[cur];
            return result;
        }

        // Convert a list of value separated by 'separator' to a std::vector
        template<class T>
        std::vector<T>
        StringToList(const std::string &value_list, const char separator) const {
            if (value_list == "") {
                return std::vector<T>();
            }

            std::size_t next_pos = std::string::npos;
            std::size_t cur_pos = 0;
            std::vector<T> result;
            do {
                if (cur_pos != 0) {  // do not choose the next value at the first iteration
                    cur_pos = next_pos + 1;
                }  // next value
                next_pos = value_list.find(separator, cur_pos);

                std::string cur_value_str = value_list.substr(cur_pos, next_pos - cur_pos);
                if (cur_pos == 0) {
                    cur_pos = next_pos + 1;
                }
                const T new_value = StringToValue<T>(cur_value_str);
                result.push_back(new_value);
            }
            while (next_pos != std::string::npos);
            return result;
        }

        std::vector<Simulation *> simulations_;
        std::vector<Material *> materials_;
        std::string calibration_file_;
    };

// Convert a string to an Eigen::Vector2d
    template<>
    inline const Eigen::Vector2d
    Application::StringToValue<Eigen::Vector2d>(const std::string &str_value) const {
        return StringToVector<double, 2>(str_value);
    }

// Convert a string to an Eigen::Vector3d
    template<>
    inline const Eigen::Vector3d
    Application::StringToValue<Eigen::Vector3d>(const std::string &str_value) const {
        return StringToVector<double, 3>(str_value);
    }

// Convert a string to an Eigen::Vector3s
    template<>
    inline const Eigen::Vector3s
    Application::StringToValue<Eigen::Vector3s>(const std::string &str_value) const {
        return StringToVector<std::size_t, 3>(str_value);
    }

}  // namespace SphaeroSim

#endif  // _SPHAEROSIM_APPLICATION_H_