/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Primoz Skraba
 *
 *    Copyright (C) 2009-2015
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

#ifndef SRC_0_PERSISTENCE_INCLUDE_GUDHI_0_PERSISTENCE_DISTANCE_H_
#define SRC_0_PERSISTENCE_INCLUDE_GUDHI_0_PERSISTENCE_DISTANCE_H_

#include <vector>
#include <cassert>

/** \brief Distance Oracle Concept
 * provide a list of neighbors
 * and return Distance
 */
template <class Iterator>
class Distance {
 public:
  /**
   *****************************************************************************************
   *  @brief      Get neighbors
   *
   *  @usage      Sets a vector of neighbors points from a query point
   * 
   *  @param[in]     queryPoint Iterator on the point from which the function finds the neighbors.
   *  @param[in,out] out        Vector of Iterator on neighbors points.
   *
   *  @return     None
   ****************************************************************************************/
  virtual void get_neighbors(Iterator queryPoint, std::vector<Iterator>& out) = 0;
};

#endif  // SRC_0_PERSISTENCE_INCLUDE_GUDHI_0_PERSISTENCE_DISTANCE_H_
