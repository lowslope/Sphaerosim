#include "Application.h"

#ifndef _INCLUDED_STRING_H_
#include <string>
#define _INCLUDED_STRING_H_
#endif

#ifndef _INCLUDED_TINYXML_2_H_

#include "tinyxml2.h"

#define _INCLUDED_TINYXML_2_H_
#endif

#ifndef _SPHAEROSIM_EXCEPTION_H_
#include "Exception.h"
#endif

#ifndef _SPHAEROSIM_SIMULATION_TESTSCENARIO_H_

#include "SimulationTestscenario.h"

#endif

#ifndef _SPHAEROSIM_SIMULATION_VTK_INPUT_H_

#include "SimulationVTKInput.h"

#endif

#ifndef _SPHAEROSIM_INTERPOLATION_H_
#include "Interpolation.h"
#endif

#ifndef _SPHAEROSIM_PROBE_H_

#include "Probe.h"

#endif

#ifndef _SPHAEROSIM_GROWTH_MODEL_H_

#include "GrowthModel.h"

#endif

#ifndef _SPHAEROSIM_NUCLEATION_MODEL_H_

#include "NucleationModel.h"

#endif

#ifndef _SPHAEROSIM_MATERIAL_H_
#include "Material.h"
#endif

#ifndef _SPHAEROSIM_HASH_VALUE_H_

#include "HashValue.h"

#endif

#ifndef _SPHAEROSIM_SIMULATION_COMSOL_INPUT_H_

#include "SimulationComsolInput.h"

#endif

namespace SphaeroSim {

    Application::Application(const std::string &parameter_file) {
        ParseParameterFile(parameter_file);
    }

    Application::~Application() {
        for (std::size_t cur = 0; cur < simulations_.size(); ++cur)
            delete simulations_[cur];
//delete[] &simulations_;

        for (std::size_t cur = 0; cur < materials_.size(); ++cur)
            delete materials_[cur];
    }

    void Application::Run() {
        for (std::size_t cur = 0; cur < simulations_.size(); ++cur)
            simulations_[cur]->Run();
    }

// Convert an error from the tinyxml2 library to a single std::string
    std::string Application::XMLErrorText(const tinyxml2::XMLDocument *xml_file) {
        if (xml_file->Error()) {
            std::string text_1 = "";
            std::string text_2 = "";
            if (xml_file->GetErrorStr1() != NULL) {
                text_1 = xml_file->GetErrorStr1();
            }
            if (xml_file->GetErrorStr2() != NULL) {
                text_2 = xml_file->GetErrorStr2();
            }

            std::string result = "Error in XML document: id=";
            result += xml_file->ErrorName();
            result += " ";
            result += "text_1=";
            result += text_1;
            result += " text_2=";
            result += text_2;
            return result;
        }
        return "";
    }

// return true if a XML value exist in the element/attribute tree of 'element'
    bool Application::XMLValueExists(const char *value,
                                     const tinyxml2::XMLElement &element) const {
        if (element.Attribute(value) != NULL) {
            return true;
        } else {
            const tinyxml2::XMLElement *value_el = element.FirstChildElement(value);
            if (value_el != NULL) {
                return true;
            }
        }
        return false;
    }

// return the value of an xml entry or attribute as std::string
// (the function assumes the element/attribute exists
    std::string Application::XMLValueString(const char *value,
                                            const tinyxml2::XMLElement &element) const {
        std::string str_value = "";
        if (element.Attribute(value) != NULL) {
            str_value = element.Attribute(value);
        } else {
            const tinyxml2::XMLElement *value_el = element.FirstChildElement(value);
            if (value_el != NULL) {
                if (value_el->GetText() == NULL) {
                    str_value = "";
                } else {
                    str_value = value_el->GetText();
                }
            } else {
                throw Exception("Value not found",
                                value);
            }
        }
        return str_value;
    }

// Parse the parameters of a Testscenario simulation run and return a
// new object
    Simulation *Application::
    CreateSimulationObject(const std::string &type,
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
                           Material *material) const {
        // parse the growth model
        const tinyxml2::XMLElement *growth_model =
                cur_element.FirstChildElement("Growthmodel");
        if (growth_model == NULL) {
            throw Exception("Simulation",
                            "The entry 'Growthmodel' is not specified.");
        }

        // parse the nucleation model
        const tinyxml2::XMLElement *nucleation_model =
                cur_element.FirstChildElement("Nucleationmodel");
        if (nucleation_model == NULL) {
            throw Exception("Simulation",
                            "The entry 'Nucleationmodel' is not specified.");
        }

        if (type == "Testscenario") {
            return
                    new SimulationTestscenario(XMLToArray<double>("temperature", cur_element),
                                               XMLToArray<double>("shearrate", cur_element),
                                               XMLToArray<double>("strainrate", cur_element),
                                               XMLToArray<Eigen::Vector2d>("temperature_profile", cur_element),
                                               XMLToArray<Eigen::Vector2d>("velocity_profile", cur_element),
                                               XMLToValue<Eigen::Vector3s>("num_cells", cur_element),
                                               name,
                                               convergence,
                                               XMLToArray<double>("timestep", cur_element),
                                               cell_size_m,
                                               temporal_interpolation,
                                               spatial_interpolation,
                                               output_folder,
                                               create_new_subfolder,
                                               vtk_options,
                                               csv_separator,
                                               XMLToArray<bool>("noslip_neighbourhood", cur_element),
                                               calibration_file_,
                                               XMLToArray<unsigned int>("rng_states", cur_element),
                                               integral_error,
                                               save_memory,
                                               calculate_flow_energy,
                                               cache_convert_temperature,
                                               cache_convert_velocity,
                                               result_step,
                                               statistics_step,
                                               numerical_options,
                                               ParseGrowthModel(material, *growth_model),
                                               ParseNucleationModel(material, *nucleation_model),
                                               material);
        } else if (type == "VTKInput") {
            return new SimulationVTKInput(XMLToValue<double>("timestep", cur_element),
                                          XMLToArray<std::string>("folders", cur_element),
                                          XMLToArray<double>("folder_timesteps", cur_element),
                                          XMLToArray<bool>("folder_start_time_zero", cur_element),
                                          XMLToArray<std::string>("folder_temperature_indicator", cur_element),
                                          XMLToArray<std::string>("folder_velocity_indicator", cur_element),
                                          XMLToArray<uint64_t>("element_ids", cur_element),
                                          name,
                                          convergence,
                                          cell_size_m,
                                          temporal_interpolation,
                                          spatial_interpolation,
                                          output_folder,
                                          create_new_subfolder,
                                          vtk_options,
                                          csv_separator,
                                          XMLToArray<bool>("noslip_neighbourhood", cur_element),
                                          calibration_file_,
                                          XMLToArray<unsigned int>("rng_states", cur_element),
                                          integral_error,
                                          save_memory,
                                          calculate_flow_energy,
                                          cache_convert_temperature,
                                          cache_convert_velocity,
                                          result_step,
                                          statistics_step,
                                          numerical_options,
                                          ParseGrowthModel(material, *growth_model),
                                          ParseNucleationModel(material, *nucleation_model),
                                          material);
        } else if (type == "ComsolInput") {
            return new SimulationComsolInput(XMLToValue<double>("timestep", cur_element),
                                             XMLToArray<std::string>("files", cur_element),
                                             XMLToValue<std::string>("CoordinateXIdentifier", cur_element),
                                             XMLToValue<std::string>("CoordinateYIdentifier", cur_element),
                                             XMLToValue<std::string>("CoordinateZIdentifier", cur_element),
                                             XMLToValue<std::string>("VelocityXIdentifier", cur_element),
                                             XMLToValue<std::string>("VelocityYIdentifier", cur_element),
                                             XMLToValue<std::string>("VelocityZIdentifier", cur_element),
                                             XMLToValue<std::string>("TemperatureIdentifier", cur_element),
                                             XMLToValue<std::string>("LevelsetIdentifier", cur_element),
                                             name,
                                             convergence,
                                             cell_size_m,
                                             temporal_interpolation,
                                             spatial_interpolation,
                                             output_folder,
                                             create_new_subfolder,
                                             vtk_options,
                                             csv_separator,
                                             XMLToArray<bool>("noslip_neighbourhood", cur_element),
                                             calibration_file_,
                                             XMLToArray<unsigned int>("rng_states", cur_element),
                                             integral_error,
                                             save_memory,
                                             calculate_flow_energy,
                                             cache_convert_temperature,
                                             cache_convert_velocity,
                                             result_step,
                                             statistics_step,
                                             numerical_options,
                                             ParseGrowthModel(material, *growth_model),
                                             ParseNucleationModel(material, *nucleation_model),
                                             material);
        } else {
            throw Exception("Siimulation:mode",
                            "Invalid simulation mode. Allowed values are 'Testscenario' and 'VTKInput'");
        }
    }

// parse an interpolation type. If the identifier is not found the function
// returns kLinear
    const InterpolationType Application::
    ParseInterpolationType(const char *value,
                           const tinyxml2::XMLElement &cur_element) const {
        if (XMLValueExists(value, cur_element) == false) {
            return kLinear;
        }

        if (XMLToValue<std::string>(value, cur_element) == "Spline") {
            return kSpline;
        }
        return kLinear;
    }

    Simulation::VTKResultOption Application::ParseVTKOptions(const tinyxml2::XMLElement &cur_element) {
        const bool save_vtk_files = XMLValueExists("VTKResults", cur_element);
        std::string vtk_subfolder = "";
        bool save_all_result_points = true;
        bool absolute_coordinates = true;
        bool temperature = true;
        bool aspectratio = true;
        bool crystallizationtime = true;
        bool phasestate = true;
        bool spheruliteid = true;
        bool spherulitediameter = true;
        bool randomcolor = true;
        bool enthropychange = true;
        bool induction = true;
        bool coolingrate = true;
        bool location = true;
        bool velocity = true;
        //new option checkboundarylayer
        bool checkboundarylayer = true;
        if (save_vtk_files) {
            const tinyxml2::XMLElement *vtk_element =
                    cur_element.FirstChildElement("VTKResults");
            if (XMLValueExists("subfolder", *vtk_element)) {
                vtk_subfolder = XMLToValue<std::string>("subfolder", *vtk_element);
            }
            if (XMLValueExists("all_result_points", *vtk_element)) {
                save_all_result_points =
                        XMLToValue<bool>("all_result_points", *vtk_element);
            }
            if (XMLValueExists("absolute", *vtk_element)) {
                absolute_coordinates = XMLToValue<bool>("absolute", *vtk_element);
            }
            if (XMLValueExists("temperature", *vtk_element)) {
                temperature = XMLToValue<bool>("temperature", *vtk_element);
            }
            if (XMLValueExists("aspectratio", *vtk_element)) {
                aspectratio = XMLToValue<bool>("aspectratio", *vtk_element);
            }
            if (XMLValueExists("crystallizationtime", *vtk_element)) {
                crystallizationtime = XMLToValue<bool>("crystallizationtime", *vtk_element);
            }
            if (XMLValueExists("phasestate", *vtk_element)) {
                phasestate = XMLToValue<bool>("phasestate", *vtk_element);
            }
            if (XMLValueExists("spheruliteid", *vtk_element)) {
                spheruliteid = XMLToValue<bool>("spheruliteid", *vtk_element);
            }
            if (XMLValueExists("spherulitediameter", *vtk_element)) {
                spherulitediameter = XMLToValue<bool>("spherulitediameter", *vtk_element);
            }
            if (XMLValueExists("randomcolor", *vtk_element)) {
                randomcolor = XMLToValue<bool>("randomcolor", *vtk_element);
            }
            if (XMLValueExists("enthropychange", *vtk_element)) {
                enthropychange = XMLToValue<bool>("enthropychange", *vtk_element);
            }
            if (XMLValueExists("induction", *vtk_element)) {
                induction = XMLToValue<bool>("induction", *vtk_element);
            }
            if (XMLValueExists("coolingrate", *vtk_element)) {
                coolingrate = XMLToValue<bool>("coolingrate", *vtk_element);
            }
            if (XMLValueExists("location", *vtk_element)) {
                location = XMLToValue<bool>("location", *vtk_element);
            }
            if (XMLValueExists("velocity", *vtk_element)) {
                velocity = XMLToValue<bool>("velocity", *vtk_element);
            }
            //new option checkboundarylayer
            if (XMLValueExists("checkboundarylayer", *vtk_element)) {
                checkboundarylayer = XMLToValue<bool>("checkboundarylayer", *vtk_element);
            }
        }

        const Simulation::VTKResultOption vtk_options(
                save_vtk_files,
                vtk_subfolder,
                save_all_result_points,
                absolute_coordinates,
                temperature,
                aspectratio,
                crystallizationtime,
                phasestate,
                spheruliteid,
                spherulitediameter,
                randomcolor,
                enthropychange,
                induction,
                coolingrate,
                location,
                velocity,
                //new option checkboundarylayer
                checkboundarylayer);

        return vtk_options;
    }

// Parse the parameters of a simulation run and add it to the 'simulations_'
// array of the application class object
    void Application::ParseSimulationRun(const tinyxml2::XMLElement &cur_element) {
        const std::string output_folder =
                XMLToValue<std::string>("output_folder", cur_element);
        const std::string mode = XMLToValue<std::string>("mode", cur_element);
        const std::string name = XMLToValue<std::string>("name", cur_element);
        const double cell_size_m = XMLToValue<double>("cell_size_m", cur_element);
        if (cell_size_m < 1e-9) {
            throw Exception("Simulation:cell_size_m",
                            "The value should not be lower than 1e-9.");
        }

        bool create_subfolder = true;
        if (XMLValueExists("create_subfolder", cur_element)) {
            create_subfolder = XMLToValue<bool>("create_subfolder", cur_element);
        }

        const InterpolationType temporal_interpolation =
                ParseInterpolationType("temporal_interpolation", cur_element);
        const InterpolationType spatial_interpolation =
                ParseInterpolationType("spatial_interpolation", cur_element);

        char csv_separator = ',';
        if (XMLValueExists("csv_separator", cur_element)) {
            csv_separator = XMLToValue<std::string>("csv_separator", cur_element)[0];
        }

        bool calculate_flow_energy = true;
        if (XMLValueExists("calculate_flow_energy", cur_element)) {
            calculate_flow_energy = XMLToValue<bool>("calculate_flow_energy", cur_element);
        }

        // VTK results
        const Simulation::VTKResultOption vtk_options = ParseVTKOptions(cur_element);

        const double integral_error = XMLToValue<double>("integral_error",
                                                         cur_element);

        bool save_memory = false;
        if (XMLValueExists("save_memory", cur_element)) {
            save_memory = XMLToValue<bool>("save_memory", cur_element);
        }

        const std::string material_id = XMLToValue<std::string>("material_id", cur_element);
        Material *mat = NULL;
        for (std::size_t iX = 0; iX < materials_.size(); ++iX) {
            if (materials_[iX]->GetName() == material_id) {
                mat = materials_[iX];
                break;
            }
        }
        if (mat == NULL) {
            throw Exception("Simulation",
                            "Material not found.");
        }

        double cache_convert_temperature = 100.0;
        if (XMLValueExists("cache_convert_temperature", cur_element)) {
            cache_convert_temperature = XMLToValue<double>("cache_convert_temperature", cur_element);
        }
        double cache_convert_velocity = 10.0;
        if (XMLValueExists("cache_convert_velocity", cur_element)) {
            cache_convert_velocity = XMLToValue<double>("cache_convert_velocity", cur_element);
        }

        // resultsteps
        const tinyxml2::XMLElement *cur_result_step =
                cur_element.FirstChildElement("ResultStep");
        Simulation::ResultStep result_step =
                cur_result_step == NULL ?
                Simulation::ResultStep(Simulation::ResultStep::kTime, -1.0) :
                Simulation::ResultStep(ParseResultStep(*cur_result_step));

        // statistics step
        cur_result_step = cur_element.FirstChildElement("StatisticsStep");
        Simulation::ResultStep statistics_step =
                cur_result_step == NULL ?
                Simulation::ResultStep(Simulation::ResultStep::kTime, -1.0) :
                ParseResultStep(*cur_result_step);

        // convergence criterium
        cur_result_step = cur_element.FirstChildElement("Convergence");
        Simulation::ResultStep convergence =
                cur_result_step == NULL ?
                Simulation::ResultStep(Simulation::ResultStep::kCrystallizationDegree, 1.0) :
                ParseResultStep(*cur_result_step);

        // numerical options
        const tinyxml2::XMLElement *cur_numerical_options =
                cur_element.FirstChildElement("NumericalOptions");
        Simulation::NumericalOptions numerical_options =
                cur_numerical_options == NULL ?
                Simulation::NumericalOptions(0, 1e-6, 0.15, 1000.0, 100000.0, 1e-3, false, false) :
                ParseNumericalOptions(*cur_numerical_options);

        simulations_.push_back(CreateSimulationObject(mode,
                                                      cur_element,
                                                      name,
                                                      convergence,
                                                      cell_size_m,
                                                      temporal_interpolation,
                                                      spatial_interpolation,
                                                      output_folder,
                                                      create_subfolder,
                                                      vtk_options,
                                                      csv_separator,
                                                      integral_error,
                                                      save_memory,
                                                      calculate_flow_energy,
                                                      cache_convert_temperature,
                                                      cache_convert_velocity,
                                                      result_step,
                                                      statistics_step,
                                                      numerical_options,
                                                      mat));

        // parse the probes
        const tinyxml2::XMLElement *cur_probe_element =
                cur_element.FirstChildElement("Probe");
        while (cur_probe_element != NULL) {
            ParseProbeDescription(*cur_probe_element,
                                  simulations_.back());
            cur_probe_element = cur_probe_element->NextSiblingElement("Probe");
        }

        // parse the predefined nuclei
        const tinyxml2::XMLElement *cur_predef_nuclei_element =
                cur_element.FirstChildElement("Predefined_nuclei");
        while (cur_predef_nuclei_element != NULL) {
            ParsePredefinedNuclei(*cur_predef_nuclei_element,
                                  simulations_.back());
            cur_predef_nuclei_element =
                    cur_predef_nuclei_element->NextSiblingElement("Predefined_nuclei");
        }
    }

// parse the numerical options
    Simulation::NumericalOptions Application::
    ParseNumericalOptions(const tinyxml2::XMLElement &cur_element) const {
        const int32_t cache_number_entries =
                XMLToValue<int32_t>("cache_number_entries", cur_element);
        const double timeresolution_cache =
                XMLToValue<double>("cache_timeresolution", cur_element);
        const double timeresolution_fdm =
                XMLToValue<double>("fdm_timeresolution", cur_element);
        const double timestep_fdm =
                XMLToValue<double>("fdm_timestep", cur_element);
        const double threshold_pure_flow =
                XMLToValue<double>("threshold_pure_flow", cur_element);
        const double threshold_incompressibility =
                XMLToValue<double>("threshold_incompressibility", cur_element);
        const int32_t consider_elongational =
                XMLToValue<int32_t>("consider_pure_elongational", cur_element);
        const int32_t consider_shear =
                XMLToValue<int32_t>("consider_pure_shear", cur_element);
        return Simulation::NumericalOptions(cache_number_entries,
                                            threshold_pure_flow,
                                            threshold_incompressibility,
                                            timeresolution_cache,
                                            timeresolution_fdm,
                                            timestep_fdm,
                                            consider_shear == 1,
                                            consider_elongational == 1);
    }

// create a new probe
    void Application::ParseProbeDescription(const tinyxml2::XMLElement &cur_element,
                                            Simulation *simulation) const {
        const std::string type = XMLToValue<std::string>("type", cur_element);
        const uint64_t element_id =
                XMLToValue<uint64_t>("element_id", cur_element);
        const std::string location_type =
                XMLToValue<std::string>("location_type", cur_element);
        const std::string name =
                XMLToValue<std::string>("name", cur_element);
        std::vector<Eigen::Vector3d> locations;
        locations = XMLToArray<Eigen::Vector3d>("location", cur_element);
        if (locations.size() < 2) {
            locations.push_back(Eigen::Vector3d(0.0, 0.0, 0.0));
        }

        simulation->AddProbeDescription(type,
                                        name,
                                        element_id,
                                        location_type,
                                        locations[0],
                                        locations[1]);
    }

    Simulation::ResultStep Application::
    ParseResultStep(const tinyxml2::XMLElement &cur_element) const {
        const std::string type = XMLToValue<std::string>("type", cur_element);
        const double stepsize = XMLToValue<double>("stepsize", cur_element);

        if (type == "Cryst") {
            return Simulation::ResultStep(Simulation::ResultStep::kCrystallizationDegree,
                                          stepsize);
        } else if (type == "Time") {
            return Simulation::ResultStep(Simulation::ResultStep::kTime,
                                          stepsize);
        } else {
            throw Exception(type,
                            "The type of the ResultStep is invalid");
        }
    }

    DisengagementTime *Application::
    ParseDisengagementTime(const tinyxml2::XMLElement &cur_element) const {
        const std::string type = XMLToValue<std::string>("type", cur_element);
        if (type == "Constant") {
            const double value = XMLToValue<double>("value", cur_element);
            return new ConstantDisengagementTime(value);
        } else if (type == "WLF") {
            const double D1 = XMLToValue<double>("D1", cur_element);
            const double D2 = XMLToValue<double>("D2", cur_element);
            const double A1 = XMLToValue<double>("A1", cur_element);
            const double A2 = XMLToValue<double>("A2", cur_element);
            return new WLFDisengagementTime(D1, D2, A1, A2);
        } else {
            throw Exception("Simulation",
                            "Invalid DisengagementTime type.");
        }
    }

// parse a single material definition
    Material *Application::
    ParseMaterial(const tinyxml2::XMLElement &cur_element) const {
        const std::string name = XMLToValue<std::string>("id", cur_element);
        const double density = XMLToValue<double>("density", cur_element);
        const double molecular_weight_entanglement =
                XMLToValue<double>("molecular_weight_entanglement", cur_element);
        const double T_m0 = XMLToValue<double>("T_m0", cur_element);
        const double T_inf = XMLToValue<double>("T_inf", cur_element);
        const double enthalpy = XMLToValue<double>("enthalpy", cur_element);

        const tinyxml2::XMLElement *disengagement_time =
                cur_element.FirstChildElement("DisengagementTime");
        if (disengagement_time == NULL) {
            throw Exception("Simulation",
                            "The entry 'DisengagementTime' is not specified.");
        }

        return new Material(name,
                            density,
                            T_m0,
                            T_inf,
                            molecular_weight_entanglement,
                            enthalpy,
                            ParseDisengagementTime(*disengagement_time));
    }

// parse the definition of the nucleation model
    NucleationModel *Application::
    ParseNucleationModel(const Material *material,
                         const tinyxml2::XMLElement &cur_element) const {
        const std::string type_str =
                XMLToValue<std::string>("type", cur_element);

        if (type_str == "constant") {
            const double rate = XMLToValue<double>("rate", cur_element);
            return new ConstantNucleationModel(rate);
        } else if (type_str == "HDL") {
            const double C = XMLToValue<double>("C", cur_element);
            const double U = XMLToValue<double>("U", cur_element);
            const double K = XMLToValue<double>("K", cur_element);
            const double n = XMLToValue<double>("n", cur_element);
            const int32_t athermal = XMLToValue<bool>("athermal", cur_element);
            const double C_ath = XMLToValue<double>("C_ath", cur_element);
            const int32_t enthalpy_correction = XMLToValue<bool>("enthalpy_correction", cur_element);
            return new HDLNucleationModel(material, C, U, K, n, enthalpy_correction != 0, athermal != 0, C_ath);
        } else {
            throw Exception("Nucleationmodel",
                            "Invalid type. Allowed values are 'constant'.");
        }
    }

// parse the definition of the growth model
    GrowthModel *Application::
    ParseGrowthModel(const Material *material,
                     const tinyxml2::XMLElement &cur_element) const {

        const std::string growth_type_str =
                XMLToValue<std::string>("growth_type", cur_element);
        const std::string type_str =
                XMLToValue<std::string>("type", cur_element);

        GrowthEvent::Type growth_type;
        if (growth_type_str == "Ray-tracing") {
            growth_type = GrowthEvent::kRayTracing;
        } else if (growth_type_str == "Monte-carlo") {
            growth_type = GrowthEvent::kMonteCarlo;
        } else if (growth_type_str == "CCG") {
            growth_type = GrowthEvent::kCCG;
        } else {
            throw Exception("Growthmodel",
                            "Invalid growth-type. Allowed values are 'Ray-tracing', " \
                    "'Monte-carlo' and 'CCG'.");
        }

        if (type_str == "constant") {
            const double speed = XMLToValue<double>("speed", cur_element);
            return new ConstantGrowthModel(growth_type,
                                           speed);
        } else if (type_str == "HM") {
            const double G0 = XMLToValue<double>("G0", cur_element);
            const double U0 = XMLToValue<double>("U0", cur_element);
            const double K_g = XMLToValue<double>("K_g", cur_element);
            return new HDLGrowthModel(growth_type,
                                      material,
                                      G0,
                                      U0,
                                      K_g);
        } else {
            throw Exception("Growthmodel",
                            "Invalid type. Allowed values are 'constant'.");
        }
    }

// parse the defintion of predefined nuclei
    void Application::
    ParsePredefinedNuclei(const tinyxml2::XMLElement &cur_element,
                          Simulation *simulation) const {
        const std::string type = XMLToValue<std::string>("type", cur_element);
        std::vector<uint64_t> element_ids =
                XMLToArray<uint64_t>("element_ids", cur_element);

        if (element_ids.size() == 0) {
            element_ids.push_back(std::numeric_limits<uint64_t>::max());
        }

        for (std::size_t cur_el = 0; cur_el < element_ids.size(); ++cur_el) {
            if (type == "List") {
                simulation->AddPredefinedNucleiList(element_ids[cur_el],
                                                    XMLToValue<std::string>("location_type",
                                                                            cur_element),
                                                    XMLToArray<Eigen::Vector3d>("locations",
                                                                                cur_element));
            } else if (type == "Distribution") {
                simulation->
                        AddPredefinedNucleiDistribution(element_ids[cur_el],
                                                        XMLToValue<double>("number_per_volume",
                                                                           cur_element));
            } else {
                throw Exception("Predefined_nuclei:type",
                                "Invalid type value. Allowed values are 'List' and 'Distribution'");
            }
        }
    }

// Parse the XML-Parameter file
    void Application::ParseParameterFile(const std::string &filename) {
        calibration_file_ = filename;

        tinyxml2::XMLDocument file_handle(true, tinyxml2::PRESERVE_WHITESPACE);
        if (file_handle.LoadFile(filename.c_str()) != tinyxml2::XML_NO_ERROR) {
            throw Exception("Parameter file",
                            XMLErrorText(&file_handle));
        }

        // load the root entry
        tinyxml2::XMLElement *root_element =
                file_handle.FirstChildElement("SphaeroSim");
        if (root_element == NULL) {
            throw Exception("Value missing in parameter file",
                            "The value <SphaeroSim> is missing.");
        }

        // material
        tinyxml2::XMLElement *cur_material =
                root_element->FirstChildElement("Material");
        if (cur_material == NULL) {
            throw Exception("Value missing in parameter file",
                            "At least one <Material> must be defined.");
        }
        do {
            materials_.push_back(ParseMaterial(*cur_material));
            cur_material = cur_material->NextSiblingElement("Material");
        }
        while (cur_material != NULL);

        // simulation runs
        tinyxml2::XMLElement *cur_element =
                root_element->FirstChildElement("Simulation");
        if (cur_element == NULL) {
            throw Exception("Value missing in parameter file",
                            "At least one <Simulation> must be defined.");
        }
        do {
            ParseSimulationRun(*cur_element);
            cur_element = cur_element->NextSiblingElement("Simulation");
        }
        while (cur_element != NULL);
    }

}
