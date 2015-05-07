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

#ifndef SRC_TOMATO_INCLUDE_GUDHI_TOMATO_ANN_ANN_GRAPH__H_
#define SRC_TOMATO_INCLUDE_GUDHI_TOMATO_ANN_ANN_GRAPH__H_

#include <ANN/ANN.h>

#include <vector>

#include "gudhi/ToMaTo/Graph.h"

//--------------------------
// when the vertex class
// has a *double coordinate
// so we can use ANN
// assumes a coord member
//--------------------------

 /** \brief Graph_ANN is an implementation of the Graph concept using ANN points.
  * Uses annkFRSearch on a ANNkd_tree to get rips neighbors.
  *  @param  ANNkd_tree* <I>kd:</I>         Persistent data structure for ANN kd tree.
  *  @param  int         <I>num_points:</I> Number of points in the data structure.
  *  @param  int         <I>dim:</I>        Dimension of points in the data structure.
  *  @param  Iterator    <I>start:</I>      Iterator on the first point of the data structure.
  *  @param  double      <I>mu:</I>         Rips parameter of the data structure.
  */
template <class Iterator>
class Graph_ANN : public Graph<Iterator> {
 protected:
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
  Graph_ANN(const Graph_ANN& ref) = delete;
  /** \brief Default operator= is disabled.*/
  Graph_ANN& operator=(const Graph_ANN& ref) = delete;

 public:
    /**
   *****************************************************************************************
   *  @brief      Constructs the tree structure.
   *
   *  @usage      Call initialize ANN kd tree structure with start_ and end_ iterators on a list of points.
   *              Sets the dimension and te rips parameter value of the class
   *  @note       num_points is computed again. 
   * 
   *  @param[in] off_file_name OFF file name and path.
   *  @param[in] mu_    Rips parameter of the data structure.
   *
   *  @return     None
   ****************************************************************************************/
Graph_ANN(const std::string& off_file_name, double mu_)
      : kd(nullptr),
      num_points(-1),
      dim(-1),
      start(nullptr),
      mu(mu_) {
  }
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
Graph_ANN(Iterator start_, Iterator end_, int dim_, double mu_)
      : kd(nullptr),
      num_points(end_ - start_),
      dim(dim_),
      start(start_),
      mu(mu_) {
    initialize(start_, end_);
  }

  /** \brief Destructor; deallocates the kd tree structure from ANN. */
  ~Graph_ANN() {
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
 *  @note          Neighbors are found from rips. Rips parameter is set on Graph_ANN constructor. 
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

    assert(kd->annkFRSearch(queryPoint->geometry.coord, mu, nb_neighb, neighb, ndist) == nb_neighb);
    for (int i = 0; i < nb_neighb; i++) {
      Iterator nit = start + neighb[i];
      if (nit != queryPoint) {
        out.push_back(nit);
      }
    }
    delete[] ndist;
    delete[] neighb;
  }
 
};

#endif  // SRC_TOMATO_INCLUDE_GUDHI_TOMATO_GRAPH_ANN__H_
