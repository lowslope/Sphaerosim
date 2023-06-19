#include "VTKFile.h"

#ifndef _INCLUDED_STRING_H_
#include <string>
#define _INCLUDED_STRING_H_
#endif

#ifndef _INCLUDED_VECTOR_H_
#include <vector>
#define _INCLUDED_VECTOR_H_
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

    VTKFile::VTKFile(const std::string &filename) :
            filename_(filename) {
        geometry_ = "";
        headerline_ = "";
        num_cells_ = 0;
        num_points_ = 0;
        dataset_ = kUndefined;
    }

    void VTKFile::SetHeaderLine(const std::string &new_headerline) {
        headerline_ = new_headerline;
    }

    void VTKFile::SetStructuredGrid(const Eigen::Vector3d offset,
                                    const Eigen::Vector3d size_m,
                                    const std::size_t num_points_x,
                                    const std::size_t num_points_y,
                                    const std::size_t num_points_z) {
        std::stringstream ss;

        // header
        dataset_ = kStructuredGrid;
        ss << "DATASET STRUCTURED_GRID" << std::endl;
        ss << "DIMENSIONS " << num_points_x;
        ss << " " << num_points_y;
        ss << " " << num_points_z << std::endl;
        ss << "POINTS " << num_points_x * num_points_y * num_points_z;
        ss << " " << "float" << std::endl;

        // coordinates
        const Eigen::Vector3d cell_size_m(
                size_m[0] / static_cast<double>(num_points_x - 1),
                size_m[1] / static_cast<double>(num_points_y - 1),
                size_m[2] / static_cast<double>(num_points_z - 1));
        for (std::size_t cur_z = 0; cur_z < num_points_z; ++cur_z) {
            for (std::size_t cur_y = 0; cur_y < num_points_y; ++cur_y) {
                for (std::size_t cur_x = 0; cur_x < num_points_x; ++cur_x) {
                    const Eigen::Vector3d coordinate(
                            offset[0] + cell_size_m[0] * cur_x,
                            offset[1] + cell_size_m[1] * cur_y,
                            offset[2] + cell_size_m[2] * cur_z);
                    ss << coordinate[0] << " " << coordinate[1] << " " << coordinate[2];
                    ss << std::endl;
                }
            }
        }
        geometry_ = ss.str();
    }

    void VTKFile::AppendFloatCellField(const std::vector<double> &cell_field,
                                       const std::string &name) {
        num_cells_ = cell_field.size();
        std::stringstream cell_field_stream;
        AppendScalarField<float>(cell_field,
                                 name,
                                 "float",
                                 &cell_field_stream);
        cell_data_.push_back(cell_field_stream.str());
    }

    void VTKFile::
    AppendVectorCellField(const std::vector<Eigen::Vector3d> &cell_field,
                          const std::string &name) {
        num_cells_ = cell_field.size();
        std::stringstream cell_field_stream;
        AppendVectorField<float>(cell_field,
                                 name,
                                 "float",
                                 &cell_field_stream);
        cell_data_.push_back(cell_field_stream.str());
    }

    void VTKFile::WriteToFile(const bool deleteThis) {
        FILE *file = fopen(filename_.c_str(), "w");

        fprintf(file, "# vtk DataFile Version 2.0\n");
        fprintf(file, "%s\n", headerline_.c_str());
        fprintf(file, "ASCII\n");
        fprintf(file, "%s", geometry_.c_str());
        if (point_data_.size() > 0) {
            fprintf(file, "POINT_DATA %zu\n", num_points_);
            for (std::size_t cur = 0; cur < point_data_.size(); ++cur)
                fprintf(file, "%s", point_data_[cur].c_str());
        }
        if (cell_data_.size() > 0) {
            fprintf(file, "CELL_DATA %zu\n", num_cells_);
            for (std::size_t cur = 0; cur < cell_data_.size(); ++cur)
                fprintf(file, "%s", cell_data_[cur].c_str());
        }
        fclose(file);
        if (deleteThis) {
            delete this;
        }
    }

    void VTKFile::ReadGeometryStringUnstructuredGrid(const std::string &firstline,
                                                     std::ifstream &file) {
        std::stringstream geometry_stream;
        std::string line = firstline;
        geometry_stream << line << std::endl;

        // number of points
        size_t num_points;
        line = line.substr(7);
        std::stringstream tmpSS(line);
        tmpSS >> num_points;
        num_points_ = num_points;

        // coordinates
        for (int32_t cur_point = 0; cur_point < num_points; ++cur_point) {
            getline(file, line);
            geometry_stream << line << std::endl;
        }

        // cell information
        while (line.find("CELLS") != 0) {
            getline(file, line);
        }
        geometry_stream << line << std::endl;

        // number of elements
        int32_t num_elements;
        line = line.substr(6);
        tmpSS.str(line);
        tmpSS >> num_elements;

        // indices
        for (int32_t cur_element = 0; cur_element < num_elements; ++cur_element) {
            getline(file, line);
            geometry_stream << line << std::endl;
        }

        // cell types
        while (line.find("CELL_TYPES") != 0) {
            getline(file, line);
        }

        for (int32_t cur_element = 0; cur_element < num_elements; ++cur_element) {
            getline(file, line);
            tmpSS.str(line);
            int32_t cell_type;
            tmpSS >> cell_type;
            if (cell_type != 11) {
                throw Exception("VTKFile:ReadFromFile",
                                "Only CELLTYPE=11 is supported for reading operations.");
            }
        }
        geometry_ = geometry_stream.str();
    }

    void VTKFile::ReadField(const std::string &firstline,
                            const std::size_t num_entries,
                            std::ifstream &file,
                            std::vector<std::string> *field_storage) {
        std::string line = firstline;
        std::stringstream field_stream;

        while (line.find("SCALARS") == std::string::npos &&
               line.find("VECTORS") == std::string::npos) {
            getline(file, line);
        }

        field_stream << line << std::endl;

        std::size_t cur_point = 0;
        while (cur_point < num_entries) {
            getline(file, line);

            // skip the LOOKUP_TABLE
            if (line.find("LOOKUP_TABLE") != std::string::npos) {
                continue;
            }
            field_stream << line << std::endl;
            ++cur_point;
        }

        field_storage->push_back(field_stream.str());
    }

    void VTKFile::ReadFromFile() {
        std::ifstream file(filename_.c_str(), std::ifstream::in);
        ClearCache();

        // ignore first two lines (Version, header)
        std::string line;
        for (int iX = 0; iX < 2; iX++)
            getline(file, line);

        // ASCII/BINARY
        getline(file, line);
        if (line != "ASCII") {
            throw Exception("VTKFile:ReadFromFile",
                            "Only ASCII files are supported.");
        }

        // DATASET
        getline(file, line);
        if (line.find("UNSTRUCTURED_GRID") != std::string::npos) {
            dataset_ = kUnstructuredGrid;
        } else {
            throw Exception("VTKFile:ReadFromFile",
                            "Only UNSTRUCTURED_GRID is supported for reading operations.");
        }

        // read the file line by line
        bool point_data = false;
        while (file.good()) {
            getline(file, line);
            if (dataset_ == kUnstructuredGrid) {
                if (line.find("POINTS") == 0) {
                    ReadGeometryStringUnstructuredGrid(line, file);
                }
            }

            if (line.find("POINT_DATA") == 0) {
                point_data = true;
            }
            if (line.find("CELL_DATA") == 0) {
                point_data = false;
            }

            if (line.find("SCALARS") != std::string::npos ||
                line.find("VECTORS") != std::string::npos) {
                if (point_data) {
                    ReadField(line, num_points_, file, &point_data_);
                } else {
                    ReadField(line, num_cells_, file, &cell_data_);
                }
            }
        }
        file.close();
    }

    const std::size_t
    VTKFile::MinimumPointIndex(const std::vector<Eigen::Vector3d> &coordinates,
                               const std::vector<std::size_t> &indices) const {
        std::size_t min_pt = indices[0];
        for (std::size_t cur_index = 1; cur_index < indices.size(); ++cur_index) {
            const std::size_t index = indices[cur_index];
            if (coordinates[index][0] - coordinates[min_pt][0] < -1e-8 ||
                coordinates[index][1] - coordinates[min_pt][1] < -1e-8 ||
                coordinates[index][2] - coordinates[min_pt][2] < -1e-8) {
                min_pt = indices[cur_index];
            }
        }
        return min_pt;
    }

// helper function to sort the corner points in an element box
//                       6--------7
//                      /|       /|                  
//                     2--------3 |                  y
//                     | |      | |                  |  z
//                     | 4------|-5                  | /
//                     |/       |/                   |/
//                     0--------1                    ------x
    const std::size_t
    VTKFile::PointIndex(const Eigen::Vector3d &point,
                        const Eigen::Vector3d &minimum_point) const {
        const double delta = 1e-8;
        if (fabs(minimum_point[2] - point[2]) < delta) { // point 1,2,3
            if (fabs(minimum_point[1] - point[1]) < delta) { // point 1
                return 1;
            }
            if (fabs(minimum_point[0] - point[0]) < delta) { // point 2
                return 2;
            }
            return 3;
        } else { // point 4,5,6,7
            if ((fabs(minimum_point[1] - point[1]) < delta) &&
                (fabs(minimum_point[0] - point[0]) < delta)) { // point 4
                return 4;
            }
            if (fabs(minimum_point[1] - point[1]) < delta) { // point 5
                return 5;
            }
            if (fabs(minimum_point[0] - point[0]) < delta) { // point 6
                return 6;
            }
            return 7;
        }
    }

    const std::vector<std::size_t>
    VTKFile::SortIndices(const std::vector<std::size_t> &indices,
                         const std::vector<Eigen::Vector3d> &coordinates,
                         const std::size_t minimum_index) const {
        std::vector<std::size_t> result;
        result.resize(8);
        result[0] = minimum_index;
        for (std::size_t cur = 0; cur < indices.size(); ++cur) {
            if (indices[cur] == minimum_index) {
                continue;
            }
            const std::size_t cur_index = PointIndex(coordinates[indices[cur]],
                                                     coordinates[minimum_index]);
            result[cur_index] = indices[cur];
        }
        return result;
    }

    void VTKFile::
    CellGeometryUnstructuredGrid(std::vector<Eigen::Vector3d> *locations,
                                 std::vector<Eigen::Vector3d> *sizes,
                                 std::vector<std::vector<std::size_t> > *indices) {
        std::stringstream stream(geometry_);

        std::string line;

        // header line
        getline(stream, line, '\n');

        // number of points
        std::size_t start_pos = line.find(' ');
        if (start_pos == std::string::npos || start_pos == line.size() - 1) {
            throw Exception("VTKFile:CellGeometryUnstructuredGrid",
                            "Invalid headerline for the geometry.");
        }
        std::size_t end_pos = line.find(' ', start_pos + 1);
        std::stringstream point_stream;
        if (end_pos != std::string::npos) {
            point_stream << line.substr(start_pos + 1, end_pos - start_pos - 1);
        } else {
            point_stream << line.substr(start_pos + 1);
        }
        point_stream >> num_points_;

        std::vector<Eigen::Vector3d> coordinates;

        // coordinates
        for (std::size_t cur = 0; cur < num_points_; ++cur) {
            getline(stream, line, '\n');

            std::stringstream coordinate_stream;
            coordinate_stream << line;

            double values[3];
            coordinate_stream >> values[0];
            coordinate_stream >> values[1];
            coordinate_stream >> values[2];

            coordinates.push_back(Eigen::Vector3d(values[0],
                                                  values[1],
                                                  values[2]));
        }

        std::string test = stream.str();

        while (line.find("CELLS") == std::string::npos) {
            getline(stream, line, '\n');
        }

        // number of cells
        start_pos = line.find(' ');
        if (start_pos == std::string::npos || start_pos == line.size() - 1) {
            throw Exception("VTKFile:CellGeometryUnstructuredGrid",
                            "Invalid headerline for the geometry.");
        }
        end_pos = line.find(' ', start_pos + 1);
        std::stringstream cell_stream;
        if (end_pos != std::string::npos) {
            cell_stream << line.substr(start_pos + 1, end_pos - start_pos - 1);
        } else {
            cell_stream << line.substr(start_pos + 1);
        }
        cell_stream >> num_cells_;

        // cell indices
        for (std::size_t cur = 0; cur < num_cells_; ++cur) {
            getline(stream, line, '\n');

            std::stringstream index_stream;
            index_stream << line;

            // number of indices
            std::size_t num_indices;
            index_stream >> num_indices;
            if (num_indices != 8) {
                throw Exception("VTKFile:CellGeometryUnstructuredGrid",
                                "Invalid number of indicesin the cell description");
            }

            // indices
            std::vector<std::size_t> cur_indices;
            cur_indices.resize(8);
            for (std::size_t cur_index = 0; cur_index < 8; ++cur_index)
                index_stream >> cur_indices[cur_index];

            // find the minimum location
            const std::size_t minimum_index = MinimumPointIndex(coordinates,
                                                                cur_indices);
            locations->push_back(coordinates[minimum_index]);

            // sort the indices
            cur_indices = SortIndices(cur_indices,
                                      coordinates,
                                      minimum_index);
            indices->push_back(cur_indices);

            // element sizes
            sizes->push_back(coordinates[cur_indices[7]] - locations->back());
        }
    }

    void VTKFile::CellGeometry(std::vector<Eigen::Vector3d> *locations,
                               std::vector<Eigen::Vector3d> *sizes,
                               std::vector<std::vector<std::size_t> > *indices) {
        locations,
                sizes;
        indices;
        if (dataset_ == kUnstructuredGrid) {
            CellGeometryUnstructuredGrid(locations,
                                         sizes,
                                         indices);
        } else {
            throw Exception("VTKFile:CellGeometry",
                            "Only UNSTRUCTURED_GRID is supported.");
        }
    }

}  // namespace SphaeroSim
