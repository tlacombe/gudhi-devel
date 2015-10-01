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

#ifndef _CLUSTER_NEIGHBORHOOD_GRAPH_H
#define _CLUSTER_NEIGHBORHOOD_GRAPH_H

#include <vector>
#include <string>

/** \brief Neighborhood graph concept for cluster computation function.*/
template <class Iterator>
class Cluster_neigborhood_graph {
 public:
  /** \brief Type for the Iterator on points of the neighborhood graph. */
  typedef Iterator Neighborhood_iterator;

  /** \brief get_neighbors is returning Iterator's neighbors in the graph excluding the given Iterator.
   *  @param[in] Iterator Iterator on the point of the graph we want to find the neighbors.
   *  @param[out] std::vector<Iterator>& Vector is fed with neighbors Iterator. */
  void get_neighbors(Iterator, std::vector<Iterator>&) { }

  /** \brief Returns the Iterator on the first point of the neighborhood graph. */
  Iterator get_start() const;
  /** \brief Returns the Iterator on the last point of the neighborhood graph. */
  Iterator get_end() const;

  /** \brief Returns a string in the format "X Y Z " composed from the Iterator coordinates.
   If dimension is less than 3, a 0 shall replace the coordinate value
   If dimension is more than 3, coordinates shall be truncated*/
  std::string get_xyz(Iterator in) const;

  /** \brief Returns an Iterator on the sink of the Iterator*/
  Iterator get_sink(Iterator in) const;

  /** \brief Sets the Iterator on the sink of the Iterator*/
  void set_sink(Iterator in, Iterator sink);

  /** \brief Returns the function value of an Iterator*/
  double get_func(Iterator in) const;
};

#endif  // _CLUSTER_NEIGHBORHOOD_GRAPH_H
