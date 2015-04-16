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

#ifndef SRC_0_PERSISTENCE_INCLUDE_GUDHI_0_PERSISTENCE_DISTANCE_ANN_H_
#define SRC_0_PERSISTENCE_INCLUDE_GUDHI_0_PERSISTENCE_DISTANCE_ANN_H_

#include <ANN/ANN.h>

#include <vector>

#include "gudhi/0_persistence/Distance.h"

//--------------------------
// when the vertex class
// has a *double coordinate
// so we can use ANN
// assumes a coord member
//--------------------------

 /** \brief Distance_ANN is an implementation of the Distance concept using ANN points.
  * Uses annkFRSearch on a ANNkd_tree to get rips neighbors.
  *  @param  ANNkd_tree* <I>kd:</I>         Persistent data structure for ANN kd tree.
  *  @param  int         <I>num_points:</I> Number of points in the data structure.
  *  @param  int         <I>dim:</I>        Dimension of points in the data structure.
  *  @param  Iterator    <I>start:</I>      Iterator on the first point of the data structure.
  *  @param  double      <I>mu:</I>         Rips parameter of the data structure.
  */
template <class Iterator>
class Distance_ANN : public Distance<Iterator> {
 private:
  //----------------------
  // persistent data structure
  // for kd tree
  //----------------------
  ANNkd_tree *kd;
  int num_points;
  int dim;
  Iterator start;
  //----------------------------------
  // rips parameter
  //----------------------------------
  double mu;

  /** \name Constructor/Destructor
   * @{ */

 private:
  /** \brief Default copy constructor is disabled.*/
  Distance_ANN(const Distance_ANN& ref) = delete;
  /** \brief Default operator= is disabled.*/
  Distance_ANN& operator=(const Distance_ANN& ref) = delete;

 public:
  /**
 *****************************************************************************************
 *  @brief      Constructs the tree structure.
 *
 *  @usage      Call initialize ANN kd tree structure with start_ and end_ iterators on a list of points.
 *              Sets the dimension and te rips parameter value of the class
 *  @note       num_points is computed again. 
 * 
 *  @param[in] start_ Iterator on the first point.
 *  @param[in] end_   Iterator on the last point.
 *  @param[in] dim_   Dimension of points in the data structure.
 *  @param[in] mu_    Rips parameter of the data structure.
 *
 *  @return     None
 ****************************************************************************************/
Distance_ANN(Iterator start_, Iterator end_, int dim_, double mu_)
      : kd(nullptr),
      num_points(end_ - start_),
      dim(dim_),
      start(start_),
      mu(mu_) {
    initialize(start_, end_);
  }

  /** \brief Destructor; deallocates the kd tree structure from ANN. */
  ~Distance_ANN() {
    delete kd;
  }

  /** @} */ // end constructor/destructor

/**
 *****************************************************************************************
 *  @brief      Initializes ANN kd tree structure
 *
 *  @usage      Constructs an ANN kd tree structure from start and end iterators on a list of points.
 *  @note       num_points is computed again. 
 * 
 *  @param start_ Iterator on the first point.
 *  @param end_   Iterator on the last point.
 *  @warning    Former ANN kd tree structure is deleted
 *
 *  @return     None
 ****************************************************************************************/
  void initialize(Iterator start_, Iterator end_) {
    num_points = end_ - start_;
    start = start_;
    int i = 0;
    //------------------
    // memory allocation
    //------------------
    double** data = new double*[num_points];
    for (; start_ != end_; start_++) {
      data[i++] = start_->geometry.coord;
    }
    if (kd != nullptr) {
      // former kd tree removal
      delete kd;
    }
    kd = new ANNkd_tree(data, num_points, dim);
    // TODO(VR): kd deletion does not delete the array of double* - memory leak to be fixed
  }


  //----------------------------------
  // get rips neighbors, this is what
  // algorithm calls
  //----------------------------------
/**
 *****************************************************************************************
 *  @brief         Get rips neighbors
 *
 *  @usage         Sets a vector of neighbors points from a query point.
 *  @note          Neighbors are found from rips. Rips parameter is set on Distance_ANN constructor. 
 * 
 *  @param[in]     queryPoint Iterator on the point from which the function finds the neighbors.
 *  @param[in,out] out        Vector of Iterator on neighbors points.
 *  @warning       Out vector is cleared on call before to be filled.
 *
 *  @return     None
 ****************************************************************************************/
  virtual void get_neighbors(Iterator queryPoint, std::vector<Iterator>& out) {
    out.clear();
    // get number of neighbors to allocate memory
    int nb_neighb = kd->annkFRSearch(queryPoint->geometry.coord, mu, 0);

    //---------------------------
    // memory efficient version
    //---------------------------
    // allocate mem each time
    //---------------------------
    int *neighb = new int[nb_neighb];
    memset(neighb, 0, sizeof (int)*nb_neighb);
    double *ndist = new double[nb_neighb];
    memset(ndist, 0, sizeof (double)*nb_neighb);

    int test = kd->annkFRSearch(queryPoint->geometry.coord, mu, nb_neighb, neighb, ndist);

    assert(test == nb_neighb);
    for (int i = 0; i < nb_neighb; i++) {
      Iterator nit = start + neighb[i];
      if (nit != queryPoint) {
        out.push_back(nit);
      }
    }
    delete[] ndist;
    delete[] neighb;
  }

  //----------------------------------
  // get number of neighbors
  // for density estimation
  //----------------------------------
/**
 *****************************************************************************************
 *  @brief      Get number of neighbors in a given radius.
 *
 *  @note       Rips parameter set on Distance_ANN constructor is not used. 
 * 
 *  @param[in]  queryPoint Iterator on the point from which the function finds the neighbors.
 *  @param[in]  radius     radius of rips parameter.
 *
 *  @return     Number of neighbors.
 ****************************************************************************************/
  int get_num_neighbors(Iterator queryPoint, double radius) {
    int nb_neighb = kd->annkFRSearch(queryPoint->geometry.coord, radius*radius, 0);
    return nb_neighb;
  }

  //----------------------------------
  // get distances to k nearest neighbors
  // for density estimation
  //----------------------------------

  void get_neighbors_dist(Iterator queryPoint, int k, double *ndist) {
    int *neighb = new int[k];
    memset(neighb, 0, sizeof (int)*k);

    kd->annkSearch(queryPoint->geometry.coord, k, neighb, ndist, 0);
    delete[] neighb;
  }

  //----------------------------------
  // get distances neighbors within r
  // for density estimation
  //----------------------------------  

  int get_neighbors_dist_r(Iterator queryPoint, double radius, double *ndist) {
    int nb_neighb = kd->annkFRSearch(queryPoint->geometry.coord, radius*radius, 0);

    int *neighb = new int[nb_neighb];
    memset(neighb, 0, sizeof (int)*nb_neighb);
    int test = kd->annkFRSearch(queryPoint->geometry.coord, radius*radius, nb_neighb, neighb, ndist, 0);

    delete[] neighb;
    return test;
  }
 
};

#endif  // SRC_0_PERSISTENCE_INCLUDE_GUDHI_0_PERSISTENCE_DISTANCE_ANN_H_
