#ifndef _COMSOL_FILE_H_
#define _COMSOL_FILE_H_

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

#ifndef _INCLUDED_SSTREAM_H_

#include <sstream>

#define _INCLUDED_SSTREAM_H_
#endif

#ifndef _INCLUDED_FSTREAM_H_

#include <fstream>

#define _INCLUDED_FSTREAM_H_
#endif

#ifndef _SPHAEROSIM_EIGEN_LIBRARY_H_

#include "EigenLibrary.h"

#endif

namespace SphaeroSim {
    class ComsolFile {
    public:
        enum Fields {
            kLevelSetVariable,
            kVelocityX,
            kVelocityY,
            kVelocityZ,
            kTemperature
        };

        ComsolFile(const std::string &filename,
                   const std::string &coordinate_x_identifier,
                   const std::string &coordinate_y_identifier,
                   const std::string &coordinate_z_identifier,
                   const std::string &velocity_x_identifier,
                   const std::string &velocity_y_identifier,
                   const std::string &velocity_z_identifier,
                   const std::string &temperature_identifier,
                   const std::string &levelset_identifier,
                   const double cell_size_in_m);

        void ReadFromFile();

        void ElementFields(const int64_t index,
                           const Fields field_type,
                           std::vector<double> *points_in_time,
                           std::vector<std::vector<double> > *values) const;

        void ElementIndices(const int64_t index,
                            std::vector<std::size_t> *indices) const;

        // getter
        inline const int64_t GetNumberElements() const {
            return static_cast<int64_t>(geometry_information_.size());
        }

        inline const Eigen::Vector3d GetElementPosition(const int64_t index) const {
            return geometry_information_[index].location;
        }

        inline const Eigen::Vector3d GetElementSize(const int64_t index) const {
            return geometry_information_[index].size;
        }

        inline bool GetHasField(const Fields field_type) const {
            return !(node_fields_.find(field_type) == node_fields_.end());
        }

        inline std::string GetFilename() const {
            return filename_;
        }

    private:
        inline void RemoveBlanks(std::string *dst) {
            dst->erase(remove_if(dst->begin(), dst->end(), ::isspace), dst->end());
        }

        struct GeometryInformation {
            Eigen::Vector3d location;
            Eigen::Vector3d size;
            Eigen::Vector3s indices[8];
        };
        struct NodeField {
            double point_in_time;
            std::vector<std::vector<std::vector<double> > > values;
        };

        bool ParseHeaderLine(const std::string &line);

        bool FindCloseIndex(const int32_t direction,
                            const Eigen::Vector3s index,
                            const std::vector<double> coordinates[3],
                            std::size_t *dst) const;

        Eigen::Vector3d LocationFromIndex(const Eigen::Vector3s &index,
                                          const std::vector<double> coordinates[3]) const;

        double NodefieldValue(const NodeField *field,
                              const Eigen::Vector3s &index) const;

        const std::string filename_;
        const double cell_size_in_m_;

        std::string coordinate_x_identifier_;
        std::string coordinate_y_identifier_;
        std::string coordinate_z_identifier_;
        std::string velocity_x_identifier_;
        std::string velocity_y_identifier_;
        std::string velocity_z_identifier_;
        std::string temperature_identifier_;
        std::string levelset_identifier_;
        int32_t dimension_;
        double to_meter_factor_;
        std::vector<GeometryInformation> geometry_information_;
        std::vector<Fields> field_order_;
        std::map<Fields, std::vector<NodeField> > node_fields_;

        // not implemented
        ComsolFile &operator=(const ComsolFile &);
    };

}  // SphaeroSim
#endif  // _COMSOL_FILE_H_
