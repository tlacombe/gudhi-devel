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

#ifndef _DENSITY_NEIGHBORHOOD_GRAPH_H_
#define _DENSITY_NEIGHBORHOOD_GRAPH_H_

#include <vector>
#include <string>

/** \brief Neighborhood graph concept for density computation function.*/
template <class Iterator>
class Density_neigborhood_graph {
 public:
  /** \brief Type for the Iterator on points of the neighborhood graph. */
  typedef Iterator Neighborhood_iterator;

  /** \brief get_num_neighbors is returning Iterator's neighbors in the graph including the given Iterator.
   *  @param[in] Iterator Iterator on the point of the graph we want to find the neighbors.
   *  @param[in] radius Radius in which the number of neighbors are found.
   *  @return num_neighbors Number of neighbors. */
  int get_num_neighbors(Iterator q, double radius) const { }

  /** \brief get_neighbors_dist instantiates an array of k-closest distance from an Iterator in the graph including the
   * given Iterator.
   *  @param[in] Iterator Iterator on the point of the graph we want to find the neighbors.
   *  @param[in] k Number of the closest neighbors to find the distance.
   *  @param[out] ndist Array of k-closest distance. */
  void get_neighbors_dist(Iterator q, int k, double *ndist) const { }

  /** \brief get_neighbors_dist_r instantiates an array of k-distances within a given radius from an Iterator in the
   * graph including the given Iterator.
   *  @param[in] Iterator Iterator on the point of the graph we want to find the neighbors.
   *  @param[in] radius Radius of the closest neighbors to find the distance.
   *  @param[out] ndist Array of k-closest distance. */
  int get_neighbors_dist_r(Iterator q, double radius, double *ndist) const { }

  /** \brief Returns the number of points in the neighborhood graph. */
  int get_num_points() const;
  /** \brief Returns the Iterator on the first point of the neighborhood graph. */
  Iterator get_start() const;
  /** \brief Returns the Iterator on the last point of the neighborhood graph. */
  Iterator get_end() const;

  /** \brief Sets the function value of an Iterator*/
  set_func(Iterator in, double func);
};

#endif  // _DENSITY_NEIGHBORHOOD_GRAPH_H_
