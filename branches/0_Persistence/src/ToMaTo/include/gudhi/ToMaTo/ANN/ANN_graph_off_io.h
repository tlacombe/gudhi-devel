/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Vincent Rouvreau
 *
 *    Copyright (C) 2015  INRIA Saclay (France)
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

#ifndef SRC_TOMATO_INCLUDE_GUDHI_TOMATO_ANN_ANN_GRAPH_OFF_IO_H_
#define SRC_TOMATO_INCLUDE_GUDHI_TOMATO_ANN_ANN_GRAPH_OFF_IO_H_

#include <string>
#include <vector>
#include <fstream>

#include "gudhi/Off_reader.h"

namespace Gudhi {

namespace ToMaTo {

/**
 *@brief Off reader visitor that can be passed to Off_reader to read a ANN_graph_complex.
 */
template<typename Complex>
class ANN_graph_off_visitor_reader {
  Complex& complex_;

 public:
  explicit ANN_graph_off_visitor_reader(Complex& complex) :
      complex_(complex) { }

  void init(int dim, int num_vertices, int num_faces, int num_edges) {
#ifdef DEBUG_TRACES
    std::cout << "init - dim=" << dim << " - num_vertices=" << num_vertices << std::endl;
#endif  // DEBUG_TRACES
    complex_.set_dimension(dim);
    complex_.reserve(num_vertices);
  }

  void point(const std::vector<double>& point) {
#ifdef DEBUG_TRACES
    std::cout << "p ";
    for (auto coordinate : point) {
      std::cout << coordinate << " | ";
    }
    std::cout << std::endl;
#endif  // DEBUG_TRACES

    // create new point and corresponding vertex
    ANN_point p(complex_.dimension());
    p.coord = new double[complex_.dimension()];
    for (int i = 0; i < complex_.dimension(); i++)
      p.coord[i] = point[i];
    Vertex< ANN_point > v(p);
    complex_.emplace_back(v);
  }

  void maximal_face(const std::vector<int>& face) {
    // for ANN Graph, only points are read
  }

  void done() {
    complex_.construct_ANN_tree();
#ifdef DEBUG_TRACES
    std::cout << "done" << std::endl;
#endif  // DEBUG_TRACES
  }
};

/**
 *@brief Class that allows to load a ANN_graph_complex from an off file.
 */
template<typename Complex>
class ANN_graph_off_reader {
 public:
  /**
   * name_file : file to read
   * read_complex : complex that will receive the file content
   * read_only_points : specify true if only the points must be read
   */
  ANN_graph_off_reader(const std::string & name_file, Complex& read_complex) : valid_(false) {
    std::ifstream stream(name_file);
    if (stream.is_open()) {
      // For ANN Graph, only points are read
        ANN_graph_off_visitor_reader<Complex> off_visitor(read_complex);
        Off_reader off_reader(stream);
        valid_ = off_reader.read(off_visitor);
    }
  }

  /**
   * return true if reading did not meet problems.
   */
  bool is_valid() const {
    return valid_;
  }

 private:
  bool valid_;
};

template<typename Complex>
class ANN_graph_off_writer {
 public:
  /**
   * name_file : file where the off will be written
   * save_complex : complex that be outputted in the file
   * for now only save triangles.
   */
  ANN_graph_off_writer(const std::string & name_file, const Complex& save_complex) {
    std::ofstream stream(name_file);
    if (stream.is_open()) {
      // OFF header
      stream << "OFF" << std::endl;
      stream.close();
    } else {
      std::cerr << "could not open file " << name_file << std::endl;
    }
  }
};

}  // namespace ToMaTo
}  // namespace Gudhi

#endif  // SRC_TOMATO_INCLUDE_GUDHI_TOMATO_ANN_ANN_GRAPH_OFF_IO_H_
