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

#ifndef _DENSITY_H_
#define _DENSITY_H_

#include <cmath>
#include <cassert>
#include <vector>

namespace Gudhi {

namespace density {

/**
 * \brief Density function class.
 *
 * \details
 * The density function can be applied on a Density_neigborhood_graph data structure, please refer to its concept for
 * more details.
 *
 */
template<class Density_neigborhood_graph>
class Density {
 public:
  typedef typename Density_neigborhood_graph::Neighborhood_iterator Iterator;

 private:
  /** \brief Pointer copy on the Density_neigborhood_graph data structure.*/
  Density_neigborhood_graph* ngbh_graph_;

  // Private default constructor for not to be used.

  Density() { }

 public:
  /** \brief Density constructor from a Density_neigborhood_graph data structure.
   *
   * @param[in] n_graph Density_neigborhood_graph data structure.
   */
  Density(Density_neigborhood_graph& n_graph) : ngbh_graph_(&n_graph) { }

  /** \brief Simple ball density estimator.
   * Iterates over all points of the graph and sets the function to -(Number of neighbors in a ball of a given radius
   * + 1) / Number of points in the graph.
   * @param[in] radius Radius of the ball to find the neighbors in.
   */
  void ball_density(double radius) {
    for (Iterator it = ngbh_graph_->get_start(); it != ngbh_graph_->get_end(); it++) {
      int nb_neighb = ngbh_graph_->get_num_neighbors(it, radius);
      ngbh_graph_->set_func(it, -(static_cast<double> (nb_neighb + 1)) /
                            static_cast<double> (ngbh_graph_->get_num_points()));
    }
  }

  /** \brief Distance to measure density estimator.
   * Iterates over all points of the graph and sets the function with the square root of (k / sum of the squared
   * distance).
   * @param[in] k Number of the closest neighbors to find the neighbors distance with.
   */
  void distance_to_density(int k) {
    double sum;
    double *ndist = new double[k];

    for (Iterator it = ngbh_graph_->get_start(); it != ngbh_graph_->get_end(); it++) {
      ngbh_graph_->get_neighbors_dist(it, k, ndist);
      sum = 0;
      for (int i = 0; i < k; ++i) {
        sum += ndist[i] * ndist[i];
      }
      ngbh_graph_->set_func(it, sqrt(static_cast<double> (k) / static_cast<double> (sum)));
    }

    delete[] ndist;
  }

  /** \brief Gaussian kernel bounded by the number of nearest neighbors density estimator.
   * Iterates over all points of the graph and sets the function to the sum of the Gaussian on the k-neighbors
   * distances.
   * @param[in] k Number of the closest neighbors to find the neighbors distance with.
   * @param[in] h Height of the Gaussian.
   */
  void gaussian_NN(int k, double h) {
    double *ndist = new double[k];

    for (Iterator it = ngbh_graph_->get_start(); it != ngbh_graph_->get_end(); it++) {
      ngbh_graph_->get_neighbors_dist(it, k, ndist);
      double sum = 0;
      for (int i = 0; i < k; ++i) {
        sum += exp(ndist[i]*(-0.5) / h);
      }
      ngbh_graph_->set_func(it, -sum);
    }

    delete[] ndist;
  }

  /** \brief Gaussian kernel bounded by cuttoff point density estimator.
   * Iterates over all points of the graph and sets the function to the sum of the Gaussian on the neighbors distances
   * in a ball of a given radius.
   * @param[in] radius Radius of the ball to find the neighbors in.
   * @param[in] h Height of the Gaussian.
   */
  void gaussian_mu(double radius, double h) {
    double *ndist = new double[ngbh_graph_->get_num_points()];

    for (Iterator it = ngbh_graph_->get_start(); it != ngbh_graph_->get_end(); it++) {
      // ndist is instantiate in get_neighbors_dist_r - needs the number of neighbors
      int nb_neighb = ngbh_graph_->get_neighbors_dist_r(it, radius, ndist);
      double sum = 0;
      for (int i = 0; i < nb_neighb; ++i) {
        sum += exp(ndist[i]*(-0.5) / h);
      }
      ngbh_graph_->set_func(it, -sum);
    }
    delete[] ndist;
  }
};

}  // namespace density

}  // namespace Gudhi

#endif  // _DENSITY_H_
