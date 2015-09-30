/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Primoz Skraba
 *
 *    Copyright (C) 2009 Primoz Skraba.  All Rights Reserved.
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

#ifndef _ANN_POINT_H_
#define _ANN_POINT_H_


#include <vector>
#include <string>
#include <limits>  // for numeric_limits<>
#include <algorithm>  // for min

namespace Gudhi {

namespace ANN_graph {

/** \brief ANN_point container class for ANN points of an ANN_graph. */
class ANN_point {
 public:
  /** \brief ANN_point Iterator definition.
   * \warning This is so that if you want to use set you must change the above type definition to take the set (i.e.
   * include a partial ordering).
   */
  typedef typename std::vector<ANN_point>::iterator Iterator;

 private:
  /** \brief Store point coordinates through a vector of double.
   * \note Points for ANN are double* types. They are retrieved from <CODE>coord_.data()<\CODE>. */
  std::vector<double> coord_;
  /** \brief Store point dimension. */
  int dim_;
  /** \brief Store point function.
   * \note Function can be set manually on each point or computed from Density functions.*/
  double func_;
  /** \brief Store point sink. */
  Iterator sink_;

 private:
  // private default constructor for not to be used.
  ANN_point() { }

 public:
  template<class CoordinateRange>
  ANN_point(CoordinateRange& coordinates, int dim, double func = std::numeric_limits<double>::quiet_NaN())
      : dim_(dim),
      func_(func) {
    assert(dim_ > 0);
    auto first = std::begin(coordinates);
    auto last = std::end(coordinates);
    if (first >= last)
      std::cerr << "ANN_point::ANN_point - Error on constructor: empty coordinates" << std::endl;
    assert(first < last);
    if ((last - first) < dim_)
      std::cerr << "ANN_point::ANN_point - Error on constructor: dimension is less than coordinates size" << std::endl;
    assert((last - first) >= dim_);

    coord_.assign(first, first + dim_);
  }

  /** \brief Sets the point sink.*/
  void set_sink(Iterator sink) {
    sink_ = sink;
  }

  /** \brief Returns the point sink.*/
  Iterator get_sink() const {
    return sink_;
  }

  /** \brief Sets the point function value.*/
  void set_func(double func) {
    func_ = func;
  }

  /** \brief Returns the point function value.*/
  double func() const {
    return func_;
  }

  /** \brief Returns the point coordinates.*/
  double* get_coord() {
    return coord_.data();
  }

  /** \brief Returns the first 3 point coordinates, completed with 0 if dimension is less than 0.*/
  std::string first_3_cordinates() const {
    std::string first_3_cordinates;
    // truncate to first 3 coordinates
    for (int i = 0; i < std::min(3, dim_); i++) {
      first_3_cordinates.append(std::to_string(coord_[i]));
      first_3_cordinates.append(" ");
    }
    // add additional 0 coordinates if necessary
    if (dim_ <= 0)
      first_3_cordinates.append("0 ");
    if (dim_ <= 1)
      first_3_cordinates.append("0 ");
    if (dim_ <= 2)
      first_3_cordinates.append("0 ");
    return first_3_cordinates;
  }

  /** \brief Less_Than structure definition in order to sort a vector of ANN_point.*/
  typedef struct {
    bool operator()(const ANN_point& a, const ANN_point& b) const {
      assert(a.dim_ == b.dim_);
      if (a.func() > b.func()) {
        return true;
      } else if (a.func() < b.func()) {
        return false;
      } else {
        for (int i = 0; i < a.dim_; i++) {
          if (a.coord_[i] < b.coord_[i]) return true;
          else if (a.coord_[i] > b.coord_[i]) return false;
        }
        return false;
      }
    }
  } Less_Than;
};

}  // namespace ANN_graph

}  // namespace Gudhi

#endif  // _ANN_POINT_H_
