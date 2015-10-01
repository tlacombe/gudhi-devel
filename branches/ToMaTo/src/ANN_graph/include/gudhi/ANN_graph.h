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

#ifndef _ANN_GRAPH_H_
#define _ANN_GRAPH_H_

#include "ANN/ANN.h"

#include <vector>
#include <string>
#include <cassert>
#include <algorithm>  // for copy

#include <gudhi/ANN_point.h>

namespace Gudhi {

namespace ANN_graph {

/** \brief ANN_graph based on Approximate Nearest Neighbor kd_tree
 * provides a list of neighbors. */
template <class Iterator>
class ANN_graph {
 private:
  //----------------------
  // persistent data structure
  // for kd_tree_ tree
  //----------------------
  ANNkd_tree *kd_tree_;
  int num_points_;
  int dim_;
  Iterator start_;
  Iterator end_;

  //----------------------------------
  // rips parameter
  //----------------------------------
  double mu_;


 public:
  /** \brief Type for the Iterator on points of the neighborhood graph. */
  typedef Iterator Neighborhood_iterator;

  /** \brief Returns the number of points in the neighborhood graph. */
  int get_num_points() const {
    return num_points_;
  }

  /** \brief Returns the Iterator on the first point of the neighborhood graph. */
  Iterator get_start() const {
    return start_;
  }

  /** \brief Returns the Iterator on the last point of the neighborhood graph. */
  Iterator get_end() const {
    return end_;
  }

  /** \brief Returns a string in the format "X Y Z " composed from the Iterator coordinates.
   If dimension is less than 3, a 0 shall replace the coordinate value
   If dimension is more than 3, coordinates shall be truncated*/
  std::string get_xyz(Iterator it) const {
    return it->first_3_cordinates();
  }

  /** \brief Returns an Iterator on the sink of the Iterator*/
  Iterator get_sink(Iterator in) const {
    return in->get_sink();
  }

  /** \brief Sets the Iterator on the sink of the Iterator*/
  void set_sink(Iterator in, Iterator sink) {
    in->set_sink(sink);
  }

  /** \brief Returns the function value of an Iterator*/
  double get_func(Iterator in) const {
    return in->func();
  }

  /** \brief Returns the function value of an Iterator*/
  void set_func(Iterator in, double func) {
    return in->set_func(func);
  }

 private:
  //----------------------------------
  // initialize data structure
  //----------------------------------

  void initialize() {
    int i = 0;
    //------------------
    // memory allocation
    //------------------
    double** data = new double*[num_points_];
    for (Iterator st = start_; st != end_; st++) {
      data[i] = new double[dim_];
      std::copy(st->get_coord(), st->get_coord() + dim_, data[i]);
      i++;
    }

    kd_tree_ = new ANNkd_tree(data, num_points_, dim_);
  }

  // Private default constructor for not to be used
  ANN_graph() { }

 public:
  ANN_graph(Iterator start, Iterator end, int _dim, double mu = 0)
      : num_points_(end - start),
      dim_(_dim),
      start_(start),
      end_(end),
      mu_(mu) {
    initialize();
  }

  ~ANN_graph() {
    // data used by ANN kd tree are also destructed.
    delete kd_tree_;
    annClose();
  }

  //----------------------------------
  // get rips neighbors, this is what
  // Cluster algorithm calls
  //----------------------------------

  /** \brief get_neighbors is returning Iterator's neighbors in the graph excluding the given Iterator.
   *  @param[in] Iterator Iterator on the point of the graph we want to find the neighbors.
   *  @param[out] std::vector<Iterator>& Vector is fed with neighbors Iterator. */
  void get_neighbors(Iterator q, std::vector<Iterator> &out) const {
    out.clear();
    int nb_neighb = kd_tree_->annkFRSearch(q->get_coord(), mu_, 0);

    //---------------------------
    // memory efficient version
    //---------------------------
    // allocate mem each time
    //---------------------------
    int* neighb = new int[nb_neighb];
    double* ndist = new double[nb_neighb];

    int test = kd_tree_->annkFRSearch(q->get_coord(), mu_, nb_neighb, neighb, ndist);
    if (test != nb_neighb)
      std::cerr << "ANN_graph::get_neighbors - Error on number of neighbors found" << std::endl;
    assert(test == nb_neighb);
    for (int i = 0; i < nb_neighb; i++) {
      Iterator nit = start_ + neighb[i];
      if (nit != q) {
        out.push_back(nit);
      }
    }

    delete[] ndist;
    delete[] neighb;
  }

  /** \brief get_neighbors is returning Iterator's neighbors in the graph including the given Iterator.
   *  @param[in] Iterator Iterator on the point of the graph we want to find the neighbors.
   *  @param[in] radius Radius in which the number of neighbors are found.
   *  @return num_neighbors Number of neighbors. */
  int get_num_neighbors(Iterator q, double radius) const {
    int nb_neighb = kd_tree_->annkFRSearch(q->get_coord(), radius*radius, 0, nullptr, nullptr, 0);
    return nb_neighb;
  }

  /** \brief get_neighbors_dist instantiates an array of k-closest distance from an Iterator in the graph including the
   * given Iterator.
   *  @param[in] Iterator Iterator on the point of the graph we want to find the neighbors.
   *  @param[in] k Number of the closest neighbors to find the distance.
   *  @param[out] ndist Array of k-closest distance. */
  void get_neighbors_dist(Iterator q, int k, double *ndist) const {
    assert(ndist != nullptr);
    int *neighb = new int[k];
    kd_tree_->annkSearch(q->get_coord(), k, neighb, ndist, 0);
    delete[] neighb;
  }

  /** \brief get_neighbors_dist_r instantiates an array of k-distances within a given radius from an Iterator in the
   * graph including the given Iterator.
   *  @param[in] Iterator Iterator on the point of the graph we want to find the neighbors.
   *  @param[in] radius Radius of the closest neighbors to find the distance.
   *  @param[out] ndist Array of the closest distance within the radius. */
  int get_neighbors_dist_r(Iterator q, double radius, double *ndist) const {
    assert(ndist != nullptr);
    int nb_neighb = get_num_neighbors(q, radius);

    int* neighb = new int[nb_neighb];
    int test = kd_tree_->annkFRSearch(q->get_coord(), radius*radius, nb_neighb, neighb, ndist, 0);

    delete[] neighb;
    return test;
  }
};

}  // namespace ANN_graph

}  // namespace Gudhi

#endif  // _ANN_GRAPH_H_
