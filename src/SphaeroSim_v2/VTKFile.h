#ifndef _SPHAEROSIM_VTKFILE_H_
#define _SPHAEROSIM_VTKFILE_H_

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

#ifndef _SPHAEROSIM_EXCEPTION_H_

#include "Exception.h"

#endif

namespace SphaeroSim {
    class VTKFile {
    public:
        explicit VTKFile(const std::string &filename);

        void SetHeaderLine(const std::string &new_headerline);

        void SetStructuredGrid(const Eigen::Vector3d offset,
                               const Eigen::Vector3d size_m,
                               const std::size_t num_points_x,
                               const std::size_t num_points_y,
                               const std::size_t num_points_z);

        void AppendFloatCellField(const std::vector<double> &cell_field,
                                  const std::string &name);

        void AppendVectorCellField(const std::vector<Eigen::Vector3d> &cell_field,
                                   const std::string &name);

        void WriteToFile(const bool deleteThis = false);

        void ReadFromFile();

        void CellGeometry(std::vector<Eigen::Vector3d> *locations,
                          std::vector<Eigen::Vector3d> *sizes,
                          std::vector<std::vector<std::size_t> > *indices);

        // Getter
        inline void GetPointFieldValues(const std::vector<std::size_t> &indices,
                                        const std::string &indicator,
                                        std::vector<double> *dst) {
            Field<double>(point_data_,
                          num_points_,
                          indices,
                          indicator,
                          dst,
                          &cached_scalar_fields_);
        }

        inline void GetPointFieldValues(const std::vector<std::size_t> &indices,
                                        const std::string &indicator,
                                        std::vector<Eigen::Vector3d> *dst) {
            Field<Eigen::Vector3d>(point_data_,
                                   num_points_,
                                   indices,
                                   indicator,
                                   dst,
                                   &cached_vector_fields_);
        }

    private:
        void ReadGeometryStringUnstructuredGrid(const std::string &firstline,
                                                std::ifstream &file);

        void ReadField(const std::string &firstline,
                       const std::size_t num_entries,
                       std::ifstream &file, std::vector<std::string> *field_storage);

        const std::size_t
        MinimumPointIndex(const std::vector<Eigen::Vector3d> &coordinates,
                          const std::vector<std::size_t> &indices) const;

        const std::size_t
        PointIndex(const Eigen::Vector3d &point,
                   const Eigen::Vector3d &minimum_point) const;

        const std::vector<std::size_t>
        SortIndices(const std::vector<std::size_t> &indices,
                    const std::vector<Eigen::Vector3d> &coordinates,
                    const std::size_t minimum_index) const;

        void
        CellGeometryUnstructuredGrid(std::vector<Eigen::Vector3d> *locations,
                                     std::vector<Eigen::Vector3d> *sizes,
                                     std::vector<std::vector<std::size_t> > *indices);

        template<class dst>
        void AppendScalarField(const std::vector<double> &field,
                               const std::string &name,
                               const std::string &datatype,
                               std::stringstream *output) const {
            std::string final_name = name;
            std::replace(final_name.begin(), final_name.end(), ' ', '_');
            (*output) << "SCALARS " << final_name << " " << datatype << " 1" << std::endl;
            (*output) << "LOOKUP_TABLE default" << std::endl;

            for (auto cur: field) {
                (*output) << static_cast<dst>(cur) << std::endl;
            }
        }

        template<class dst>
        void AppendVectorField(const std::vector<Eigen::Vector3d> &field,
                               const std::string &name,
                               const std::string &datatype,
                               std::stringstream *output) {
            std::string final_name = name;
            std::replace(final_name.begin(), final_name.end(), ' ', '_');
            (*output) << "VECTORS " << final_name << " " << datatype << std::endl;

            for (std::size_t cur = 0; cur < field.size(); ++cur) {
                (*output) << static_cast<dst>(field[cur][0]) << " ";
                (*output) << static_cast<dst>(field[cur][1]) << " ";
                (*output) << static_cast<dst>(field[cur][2]) << std::endl;
            }
        }

        template<class T>
        const T SingleFieldValue(const std::string &text) const {
            std::stringstream ss;
            ss << text;
            T result;
            ss >> result;
            return result;
        }

        template<class T>
        void FieldValueList(const std::string &field_string,
                            const std::size_t num_values,
                            std::vector<T> *dst) const {
            std::stringstream stream(field_string);

            std::string line;
            // header line
            getline(stream, line, '\n');

            for (std::size_t cur = 0; cur < num_values; ++cur) {
                getline(stream, line, '\n');
                dst->push_back(SingleFieldValue<T>(line));
            }
        }

        template<class T>
        void Field(const std::vector<std::string> &field_data,
                   const std::size_t num_values,
                   const std::vector<std::size_t> &indices,
                   const std::string &indicator,
                   std::vector<T> *field,
                   std::map<std::string, std::vector<T> > *cache) {
            if (cache->find(indicator) == cache->end()) {
                // find the correct field
                bool found = false;
                for (std::size_t cur = 0; cur < field_data.size(); ++cur) {
                    std::string line = field_data[cur].substr(0, field_data[cur].find('\n'));
                    if (line.find(indicator) != std::string::npos) {
                        FieldValueList<T>(field_data[cur],
                                          num_values,
                                          &(*cache)[indicator]);
                        found = true;
                        break;
                    }
                }
                if (found == false) {
                    throw Exception("The field does not exist in the vtkfile.",
                                    indicator);
                }
            }

            // copy the correct values
            for (std::size_t cur = 0; cur < indices.size(); ++cur)
                field->push_back((*cache)[indicator][indices[cur]]);
        }

        inline void ClearCache() {
            cached_scalar_fields_.clear();
            cached_vector_fields_.clear();
        }


        enum Dataset {
            kUnstructuredGrid,
            kStructuredGrid,
            kUndefined,
        };
        const std::string filename_;

        Dataset dataset_;
        std::string headerline_;
        std::string geometry_;
        std::vector<std::string> cell_data_;
        std::vector<std::string> point_data_;
        std::size_t num_cells_;
        std::size_t num_points_;

        std::map<std::string, std::vector<double> > cached_scalar_fields_;
        std::map<std::string, std::vector<Eigen::Vector3d> > cached_vector_fields_;

        // not implemented
        VTKFile &operator=(const VTKFile &);
    };

    template<>
    inline const Eigen::Vector3d
    VTKFile::SingleFieldValue<Eigen::Vector3d>(const std::string &text) const {
        std::stringstream ss;
        ss << text;
        double values[3];
        ss >> values[0];
        ss >> values[1];
        ss >> values[2];
        return Eigen::Vector3d(values[0],
                               values[1],
                               values[2]);
    }
}  // namespace SphaeroSim

#endif  // _SPHAEROSIM_VTKFILE_H_