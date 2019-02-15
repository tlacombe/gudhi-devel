/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Vincent Rouvreau
 *
 *    Copyright (C) 2019 Inria
 *
 *    This program is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *
 *    This program is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef BITMAP_CUBICAL_COMPLEX_READ_PERSEUS_FILE_H_
#define BITMAP_CUBICAL_COMPLEX_READ_PERSEUS_FILE_H_

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>

namespace Gudhi {

namespace cubical_complex {

void goto_next_uncommented_line(std::ifstream& input_file, std::string& line) {
  while (line.length() == 0 || line[0] == '#') {
    std::getline(input_file, line);
    if (input_file.eof()) {
      std::string error_str("goto_next_uncommented_line - End of file reached.");
      std::cerr << error_str << std::endl;
      throw std::invalid_argument(error_str);
    }
  }
}

template <typename T>
void read_perseus_style_file(const std::string& perseus_style_file,
                             std::vector<unsigned>& sizes_in_following_directions,
                             std::vector<T>& top_dimensional_cells,
                             std::vector<bool>& directions_in_which_periodic_b_cond_are_to_be_imposed) {
  std::ifstream input_file(perseus_style_file);

  if (!input_file.is_open()) {
    std::string error_str("read_perseus_style_file - Unable to open file ");
    error_str.append(perseus_style_file);
    std::cerr << error_str << std::endl;
    throw std::invalid_argument(error_str);
  }

  std::string line;

  // Read data dimension
  unsigned data_dimension;
  goto_next_uncommented_line(input_file, line);

  int n = sscanf(line.c_str(), "%u", &data_dimension);
  if (n != 1) {
    std::string perseus_error("Bad Perseus file format. This line is incorrect : " + line);
    std::cerr << perseus_error << std::endl;
    throw std::ios_base::failure(perseus_error.c_str());
  }
  // Reset when line is read
  line = "";
#ifdef DEBUG_TRACES
  std::cerr << "data_dimension : " << data_dimension << std::endl;
#endif  // DEBUG_TRACES

  // Read dimensions vector
  sizes_in_following_directions.reserve(data_dimension);
  // all dimensions multiplied
  std::size_t dimensions = 1;
  directions_in_which_periodic_b_cond_are_to_be_imposed.reserve(data_dimension);

  for (unsigned dimension_counter = 0; dimension_counter < data_dimension; dimension_counter++) {
    goto_next_uncommented_line(input_file, line);
    int size_in_this_dimension;
    int n = sscanf(line.c_str(), "%d", &size_in_this_dimension);
    if (n != 1) {
      std::string perseus_error("Bad Perseus file format. This line is incorrect : " + line);
      std::cerr << perseus_error << std::endl;
      throw std::ios_base::failure(perseus_error.c_str());
    }
    // Reset when line is read
    line = "";
#ifdef DEBUG_TRACES
    std::cerr << "size_in_this_dimension[" << dimension_counter << "] : " << size_in_this_dimension << std::endl;
#endif  // DEBUG_TRACES
    // true for periodic
    directions_in_which_periodic_b_cond_are_to_be_imposed.push_back(size_in_this_dimension < 0);
    std::size_t abs_size = std::abs(size_in_this_dimension);
    sizes_in_following_directions.push_back(abs_size);
    dimensions *= abs_size;
  }

  // Read filtrations vector
  top_dimensional_cells.reserve(dimensions);

  for (std::size_t filtration_counter = 0; filtration_counter < dimensions; filtration_counter++) {
    goto_next_uncommented_line(input_file, line);
    T filtration_level;
    int n = sscanf(line.c_str(), "%lf", &filtration_level);
    if (n != 1) {
      std::string perseus_error("Bad Perseus file format. This line is incorrect : " + line);
      std::cerr << perseus_error << std::endl;
      throw std::ios_base::failure(perseus_error.c_str());
    }
    // Reset when line is read
    line = "";
#ifdef DEBUG_TRACES
    std::cerr << "filtration_level[" << filtration_counter << "] : " << filtration_level << std::endl;
#endif  // DEBUG_TRACES
    top_dimensional_cells.push_back(filtration_level);
  }

  // Outputs on standard output all remaining lines
  while (!input_file.eof()) {
    std::getline(input_file, line);
    if (line.length() != 0 && line[0] != '#') {
      std::string perseus_error("Bad Perseus file format. This line will not be read : " + line);
      std::cerr << perseus_error << std::endl;
      throw std::ios_base::failure(perseus_error.c_str());
    }
  }

  input_file.close();
}

}  // namespace cubical_complex

namespace Cubical_complex = cubical_complex;

}  // namespace Gudhi

#endif  // BITMAP_CUBICAL_COMPLEX_READ_PERSEUS_FILE_H_
