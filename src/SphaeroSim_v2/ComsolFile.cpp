#include "ComsolFile.h"

#ifndef _INCLUDED_STRING_H_
#include <string>
#define _INCLUDED_STRING_H_
#endif

#ifndef _INCLUDED_VECTOR_H_
#include <vector>
#define _INCLUDED_VECTOR_H_
#endif

#ifndef _INCLUDED_SET_H_

#include <set>

#define _INCLUDED_SET_H_
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

    ComsolFile::ComsolFile(const std::string &filename,
                           const std::string &coordinate_x_identifier,
                           const std::string &coordinate_y_identifier,
                           const std::string &coordinate_z_identifier,
                           const std::string &velocity_x_identifier,
                           const std::string &velocity_y_identifier,
                           const std::string &velocity_z_identifier,
                           const std::string &temperature_identifier,
                           const std::string &levelset_identifier,
                           const double cell_size_in_m) :
            filename_(filename),
            cell_size_in_m_(cell_size_in_m) {
        coordinate_x_identifier_ = coordinate_x_identifier;
        coordinate_y_identifier_ = coordinate_y_identifier;
        coordinate_z_identifier_ = coordinate_z_identifier;
        velocity_x_identifier_ = velocity_x_identifier;
        velocity_y_identifier_ = velocity_y_identifier;
        velocity_z_identifier_ = velocity_z_identifier;
        temperature_identifier_ = temperature_identifier;
        levelset_identifier_ = levelset_identifier;

        std::string *ids[] = {&coordinate_x_identifier_,
                              &coordinate_y_identifier_,
                              &coordinate_z_identifier_,
                              &velocity_x_identifier_,
                              &velocity_y_identifier_,
                              &velocity_z_identifier_,
                              &temperature_identifier_,
                              &levelset_identifier_};
        for (std::size_t cur = 0; cur < sizeof(ids) / sizeof(std::string *); ++cur) {
            RemoveBlanks(ids[cur]);
            std::transform(ids[cur]->begin(),
                           ids[cur]->end(),
                           ids[cur]->begin(),
                           ::toupper);
        }
        dimension_ = -1;
        to_meter_factor_ = -1.0;
    }

    void ComsolFile::ElementIndices(const int64_t index,
                                    std::vector<std::size_t> *indices) const {
        for (int32_t cur = 0; cur < 8; ++cur) {
            union IndexUnion {
                int16_t words[4];
                int64_t raw_data;
            } index_union;
            Eigen::Vector3s cur_index = geometry_information_[index].indices[cur];

            for (int32_t dim = 0; dim < 3; ++dim) {
                index_union.words[dim] = static_cast<int16_t>(cur_index[dim]);
            }
            index_union.words[3] = 0;
            indices->push_back(static_cast<std::size_t>(index_union.raw_data));
        }
    }

    void ComsolFile::ReadFromFile() {
        std::ifstream file(filename_.c_str(), std::ifstream::in);
        if (file.bad() || file.is_open() == false) {
            throw Exception("File open error", filename_);
        }

        std::string curline;
        bool is_in_header = true;

        // header
        while (is_in_header == true) {
            getline(file, curline);
            is_in_header = ParseHeaderLine(curline);
            if (file.eof() == true) {
                throw Exception("Comsol File", "Unexpected end of file");
            }
        }

        if (dimension_ == -1 || dimension_ > 3 || dimension_ < 2) {
            throw Exception("Comsol File", "Invalid dimension or dimension no specified");
        }
        if (to_meter_factor_ < 0.0) {
            printf("Lengthunit not specifid, defaulting to 'mm'");
            to_meter_factor_ = 1e-3;
        }

        // got to %Grid
        for (;;) {
            RemoveBlanks(&curline);
            if (curline.find("%Grid") != std::string::npos) {
                break;
            }
            getline(file, curline);
            if (file.eof() == true) {
                throw Exception("Comsol File", "Unexpected end of file");
            }
        }

        // parse the geometry
        std::vector<double> coordinates[3];
        size_t pos;
        getline(file, curline);
        if (file.eof() == true) {
            throw Exception("Comsol File", "Unexpected end of file");
        }
        for (int32_t cur_dim = 0; cur_dim < dimension_;) {
            std::vector<double> *dst = &coordinates[cur_dim];
            pos = curline.find("\t");
            std::stringstream tmpSS(curline.substr(0,pos));
            double new_value;
            tmpSS >> new_value;
            dst->push_back(new_value * to_meter_factor_);
            if(pos != std::string::npos){
                curline.erase(0, pos+1);
            }
            else{
                ++cur_dim;
                if(cur_dim < dimension_){
                    getline(file, curline);
                }
            }
            /*
            if (curline[0] == '%') {
                const std::string identifier = curline.substr(1);
                std::vector<double> *dst;
                if (identifier == coordinate_x_identifier_) {
                    dst = &coordinates[0];
                } else if (identifier == coordinate_y_identifier_) {
                    dst = &coordinates[1];
                } else if (identifier == coordinate_z_identifier_) {
                    dst = &coordinates[2];
                } else {
                    throw Exception("Comsol File - Invalid coordinate identifier", identifier);
                }

                // read in the values
                do {
                    getline(file, curline);
                    RemoveBlanks(&curline);
                    if (file.eof() == true) {
                        throw Exception("Comsol File", "Unexpected end of file");
                    }
                    if (curline[0] != '%') {
                        std::stringstream tmpSS(curline);
                        double new_value;
                        tmpSS >> new_value;
                        dst->push_back(new_value * to_meter_factor_);
                    }
                }
                while (curline[0] != '%');
            } else {
                getline(file, curline);
                RemoveBlanks(&curline);
                if (file.eof() == true) {
                    throw Exception("Comsol File", "Unexpected end of file");
                }
            }
            */
        }
        // ensure there are enough points for a 3D mesh
        if (coordinates[0].size() == 0 && coordinates[1].size() == 0 && coordinates[2].size() == 0) {
            throw Exception("Comsol File", "No coordinates specified");
        }

        // create a mesh from the coordinates
        int32_t number_nodes = 1;
        for (int32_t cur = 0; cur < 3; ++cur) {
            if (coordinates[cur].size() != 0) {
                number_nodes *= static_cast<int32_t>(coordinates[cur].size());
            }
        }

        // read in the node information
        while (!file.eof()) {
            for (std::size_t cur_field = 0; cur_field < field_order_.size(); ++cur_field) {
                if (curline != "%Data") {
                    // % Data
                    getline(file, curline);
                    RemoveBlanks(&curline);
                    if (file.eof() == true) {
                        break;
                    }
                    if (curline != "%Data") {
                        throw Exception("Comsol File", "Expected '% Data' in the file");
                    }
                }

                // % ... @t=XX
                getline(file, curline);
                RemoveBlanks(&curline);
                if (file.eof() == true) {
                    throw Exception("Comsol File", "Unexpected end of file");
                }
                const std::size_t time_pos = curline.find("@t=");
                if (time_pos == std::string::npos) {
                    throw Exception("Comsol File", "Expected '@ t=' for point in time identification");
                }
                const std::string point_in_time_str = curline.substr(time_pos + 3);
                std::stringstream tmpSS(point_in_time_str);
                double point_in_time;
                tmpSS >> point_in_time;

                // reserve memory for all nodes
                std::vector<std::vector<std::vector<double> > > node_values;
                node_values.resize(coordinates[0].size() == 0 ? 1 : coordinates[0].size());
                for (std::size_t cur_x = 0; cur_x < node_values.size(); ++cur_x) {
                    node_values[cur_x].resize(coordinates[1].size() == 0 ? 1 : coordinates[1].size());
                    for (std::size_t cur_y = 0; cur_y < node_values[cur_x].size(); ++cur_y)
                        node_values[cur_x][cur_y].resize(coordinates[2].size() == 0 ? 1 : coordinates[2].size());
                }

                // read in the values
                int32_t node_counter = 0;
                int32_t cur_y_index = 0;
                int32_t cur_z_index = 0;

                while (node_counter != number_nodes) {
                    getline(file, curline, '\n');
                    std::stringstream valSS(curline);
                    for (std::size_t num_x_val = 0; num_x_val < node_values.size(); ++num_x_val) {
                        double cur_val;
                        valSS >> cur_val;
                        node_values[num_x_val][cur_y_index][cur_z_index] = cur_val;
                        ++node_counter;
                    }
                    ++cur_y_index;
                    if (cur_y_index == node_values[0].size()) {
                        cur_y_index = 0;
                        ++cur_z_index;
                    }
                }

                // store the new field
                NodeField new_field;
                new_field.point_in_time = point_in_time;
                new_field.values = node_values;
                node_fields_[field_order_[cur_field]].push_back(new_field);
            }
        }

        // create the mesh
        std::size_t num_points_x = coordinates[0].size() == 0 ? 1 : coordinates[0].size();
        std::size_t num_points_y = coordinates[1].size() == 0 ? 1 : coordinates[1].size();
        std::size_t num_points_z = coordinates[2].size() == 0 ? 1 : coordinates[2].size();

        std::size_t extract_point_x = 3;
        std::size_t extract_point_y = 3;

        // reduce data for x
        if (extract_point_x > 0) {
            for (std::size_t cur = num_points_x - 1; cur >= 1; --cur) {
                if (cur != extract_point_x) {
                    coordinates[0].erase(coordinates[0].begin() + cur);
                    auto field_it = node_fields_.begin();
                    while (field_it != node_fields_.end()) {
                        auto value_x_it = (*field_it).second.begin();
                        while (value_x_it != (*field_it).second.end()) {
                            (*value_x_it).values.erase((*value_x_it).values.begin() + cur);
                            ++value_x_it;
                        }
                        ++field_it;
                    }
                }
            }
            num_points_x = 2;
        }
        // reduce data for y
        if (extract_point_y > 0) {
            for (std::size_t cur = num_points_y - 1; cur >= 1; --cur) {
                if (cur != extract_point_y) {
                    coordinates[1].erase(coordinates[1].begin() + cur);
                    auto field_it = node_fields_.begin();
                    while (field_it != node_fields_.end()) {
                        auto value_x_it = (*field_it).second.begin();
                        while (value_x_it != (*field_it).second.end()) {
                            auto value_y_it = (*value_x_it).values.begin();
                            while (value_y_it != (*value_x_it).values.end()) {
                                (*value_y_it).erase((*value_y_it).begin() + cur);
                                ++value_y_it;
                            }
                            ++value_x_it;
                        }
                        ++field_it;
                    }
                }
            }
            num_points_y = 2;
        }
        // reduce data for z
        std::size_t skip_point_z = 2;
        if (skip_point_z > 0) {
            std::set<std::size_t> remove_points;
            std::size_t cur_point = 1;
            std::size_t skip_counter = skip_point_z;
            while (cur_point < num_points_z / 2) {
                if (skip_counter == 0) {
                    skip_counter = skip_point_z;
                } else {
                    remove_points.insert(cur_point);
                    std::size_t mirror_point = (num_points_z - 1) - cur_point;
                    if (mirror_point > num_points_z / 2) {
                        remove_points.insert(mirror_point);
                    }
                    skip_counter -= 1;
                }
                cur_point += 1;
            }

            if (skip_counter != 0) {
                remove_points.insert(num_points_z / 2);
            }

            auto remove_it = remove_points.rbegin();
            while (remove_it != remove_points.rend()) {
                coordinates[2].erase(coordinates[2].begin() + *remove_it);
                auto field_it = node_fields_.begin();
                while (field_it != node_fields_.end()) {
                    auto value_x_it = (*field_it).second.begin();
                    while (value_x_it != (*field_it).second.end()) {
                        auto value_y_it = (*value_x_it).values.begin();
                        while (value_y_it != (*value_x_it).values.end()) {
                            auto value_z_it = (*value_y_it).begin();
                            while (value_z_it != (*value_y_it).end()) {
                                (*value_z_it).erase((*value_z_it).begin() + *remove_it);
                                ++value_z_it;
                            }
                            ++value_y_it;
                        }
                        ++value_x_it;
                    }
                    ++field_it;
                }
                ++remove_it;
            }

            num_points_z -= remove_points.size();
        }

        //                       6--------7
        //                      /|       /|
        //                     2--------3 |                  y
        //                     | |      | |                  |  z
        //                     | 4------|-5                  | /
        //                     |/       |/                   |/
        //                     0--------1                    ------x
        for (std::size_t cur_x = 0; cur_x < num_points_x; ++cur_x) {
            for (std::size_t cur_y = 0; cur_y < num_points_y; ++cur_y) {
                for (std::size_t cur_z = 0; cur_z < num_points_z; ++cur_z) {
                    GeometryInformation info;
                    info.indices[0] = Eigen::Vector3s(cur_x, cur_y, cur_z);
                    std::size_t close_index;

                    if (FindCloseIndex(0, info.indices[0], coordinates, &close_index) == false) {
                        continue;
                    }
                    info.indices[1] = Eigen::Vector3s(close_index, cur_y, cur_z);

                    if (FindCloseIndex(1, info.indices[0], coordinates, &close_index) == false) {
                        continue;
                    }
                    info.indices[2] = Eigen::Vector3s(cur_x, close_index, cur_z);

                    if (FindCloseIndex(0, info.indices[2], coordinates, &close_index) == false) {
                        continue;
                    }
                    info.indices[3] = Eigen::Vector3s(close_index, info.indices[2][1], cur_z);

                    if (FindCloseIndex(2, info.indices[0], coordinates, &close_index) == false) {
                        continue;
                    }
                    info.indices[4] = Eigen::Vector3s(cur_x, cur_y, close_index);

                    if (FindCloseIndex(0, info.indices[4], coordinates, &close_index) == false) {
                        continue;
                    }
                    info.indices[5] = Eigen::Vector3s(close_index, info.indices[4][1], info.indices[4][2]);

                    if (FindCloseIndex(1, info.indices[4], coordinates, &close_index) == false) {
                        continue;
                    }
                    info.indices[6] = Eigen::Vector3s(info.indices[4][0], close_index, info.indices[4][2]);

                    if (FindCloseIndex(0, info.indices[6], coordinates, &close_index) == false) {
                        continue;
                    }
                    info.indices[7] = Eigen::Vector3s(close_index, info.indices[6][1], info.indices[6][2]);

                    info.location = LocationFromIndex(info.indices[0], coordinates);
                    info.size = LocationFromIndex(info.indices[7], coordinates) -
                                LocationFromIndex(info.indices[0], coordinates);
                    geometry_information_.push_back(info);
                }
            }
        }
    }

    bool ComsolFile::FindCloseIndex(const int32_t direction,
                                    const Eigen::Vector3s index,
                                    const std::vector<double> coordinates[3],
                                    std::size_t *dst) const {
        int32_t start_val = static_cast<int32_t>(index[direction]);
        int32_t number_coordinates = static_cast<int32_t>(coordinates[direction].size());
        if (number_coordinates == 0 || number_coordinates == 1) {
            if (index[direction] == 0) {
                (*dst) = 1;
                return true;
            }
            return false;
        }

        Eigen::Vector3d cur_location = LocationFromIndex(index, coordinates);
        double cur_min_dist = 1e150;
        bool has_one = false;
        // positive direction
        Eigen::Vector3s test_index = index;
        for (int32_t cur = start_val + 1; cur < number_coordinates; ++cur) {
            test_index[direction] = cur;
            Eigen::Vector3d test_location = LocationFromIndex(test_index, coordinates);
            if (test_location[direction] > cur_location[direction]) {
                if (test_location[direction] - cur_location[direction] < cur_min_dist) {
                    cur_min_dist = test_location[direction] - cur_location[direction];
                    *dst = cur;
                    has_one = true;
                }
            }
        }

        // negative direction
        for (int32_t cur = start_val - 1; cur >= 0; --cur) {
            test_index[direction] = cur;
            Eigen::Vector3d test_location = LocationFromIndex(test_index, coordinates);
            if (test_location[direction] > cur_location[direction]) {
                if (test_location[direction] - cur_location[direction] < cur_min_dist) {
                    cur_min_dist = test_location[direction] - cur_location[direction];
                    *dst = cur;
                    has_one = true;
                }
            }
        }
        return has_one;
    }

    Eigen::Vector3d ComsolFile::
    LocationFromIndex(const Eigen::Vector3s &index,
                      const std::vector<double> coordinates[3]) const {
        Eigen::Vector3d location;
        for (int32_t cur = 0; cur < 3; ++cur) {
            if (coordinates[cur].size() == 0 || coordinates[cur].size() == 1) {
                if (index[cur] == 0) {
                    location[cur] = 0.0;
                } else {
                    location[cur] = cell_size_in_m_;
                }
            } else {
                location[cur] = coordinates[cur][index[cur]];
            }
        }
        return location;
    }

    bool ComsolFile::ParseHeaderLine(const std::string &line) {
        std::string to_parse = line;

        RemoveBlanks(&to_parse);
        if (to_parse.length() < 3) {
            return false;
        }

        if (to_parse[0] == '%') {
            // extract the identifier between '%' and ':'
            std::size_t double_point = to_parse.find(':');
            if (double_point == std::string::npos) {
                return false;
            }
            const std::size_t id_length = double_point - 1;
            std::string id = to_parse.substr(1, id_length);
            std::transform(id.begin(), id.end(), id.begin(), ::toupper);

            // extract the value
            const std::size_t value_length = to_parse.size() - double_point - 1;
            const std::string value = to_parse.substr(double_point + 1, value_length);

            if (id == "MODEL") {
                return true;  // model name not used
            } else if (id == "VERSION") {
                return true;  // version not used
            } else if (id == "DATE") {
                return true;  // Date not used
            } else if (id == "DIMENSION") {
                std::stringstream tmpSS(value);
                tmpSS >> dimension_;
                return true;
            } else if (id == "NODES") {
                return true;  // Nodes not used
            } else if (id == "EXPRESSIONS") {
                return true;  // Expressions not used
            } else if (id == "DESCRIPTION") {
                std::size_t cur_pos = 0;
                // read from ',' to ','
                while (cur_pos != std::string::npos) {
                    std::size_t end_pos = value.find(',', cur_pos);
                    std::string current_value =
                            end_pos != std::string::npos ?
                            value.substr(cur_pos, end_pos - cur_pos) :
                            value.substr(cur_pos, end_pos);
                    std::transform(current_value.begin(),
                                   current_value.end(),
                                   current_value.begin(),
                                   ::toupper);
                    if (current_value == levelset_identifier_) {
                        field_order_.push_back(kLevelSetVariable);
                    } else if (current_value == velocity_x_identifier_) {
                        field_order_.push_back(kVelocityX);
                    } else if (current_value == velocity_y_identifier_) {
                        field_order_.push_back(kVelocityY);
                    } else if (current_value == velocity_z_identifier_) {
                        field_order_.push_back(kVelocityZ);
                    } else if (current_value == temperature_identifier_) {
                        field_order_.push_back(kTemperature);
                    } else {
                        throw Exception("Comsol File - Unkown", current_value);
                    }

                    cur_pos = end_pos;
                    if (end_pos != std::string::npos) {
                        cur_pos += 1;
                    }
                }
                return true;
            } else if (id == "LENGTHUNIT") {
                if (value == "mm") {
                    to_meter_factor_ = 1e-3;
                } else if (value == "m") {
                    to_meter_factor_ = 1.0;
                } else {
                    throw Exception("Comsol File", "Unsupported length unit. Should be 'mm' or 'm'");
                }
                return true;
            }
        }
        return false;
    }

    void ComsolFile::ElementFields(const int64_t index,
                                   const Fields field_type,
                                   std::vector<double> *points_in_time,
                                   std::vector<std::vector<double> > *values) const {
        values->resize(node_fields_.at(field_type).size());
        points_in_time->resize(values->size());
        for (std::size_t cur_set = 0; cur_set < values->size(); ++cur_set) {
            const NodeField *cur_field = &node_fields_.at(field_type)[cur_set];
            (*points_in_time)[cur_set] = cur_field->point_in_time;

            (*values)[cur_set].resize(8);
            for (int32_t cur = 0; cur < 8; ++cur) {
                Eigen::Vector3s cur_index = geometry_information_[index].indices[cur];
                (*values)[cur_set][cur] = NodefieldValue(cur_field, cur_index);
            }
        }
    }

    double ComsolFile::NodefieldValue(const NodeField *field,
                                      const Eigen::Vector3s &index) const {
        std::size_t indices[3];

        indices[0] = index[0] >= field->values.size() ?
                     field->values.size() - 1 :
                     index[0];
        indices[1] = index[1] >= field->values[indices[0]].size() ?
                     field->values[indices[0]].size() - 1 :
                     index[1];
        indices[2] = index[2] >= field->values[indices[0]][indices[1]].size() ?
                     field->values[indices[0]][indices[1]].size() - 1 :
                     index[2];
        return field->values[indices[0]][indices[1]][indices[2]];
    }

}  // namespace SphaeroSim