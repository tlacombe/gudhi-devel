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

#ifndef SRC_TOMATO_INCLUDE_GUDHI_TOMATO_ANN_ANN_GRAPH_WITH_DENSITY_ANN__H_
#define SRC_TOMATO_INCLUDE_GUDHI_TOMATO_ANN_ANN_GRAPH_WITH_DENSITY__H_

#include <vector>

#include "gudhi/ToMaTo.h"
#include "gudhi/ToMaTo/ANN/ANN_Graph.h"

//--------------------------
// when the vertex class
// has a *double coordinate
// so we can use ANN
// assumes a coord member
//--------------------------

/** \brief ANN_Graph_with_density is an implementation of the ANN_Graph.
 */
template <class Iterator>
class ANN_Graph_with_density : public Graph_ANN<Iterator> {
 private:
  /** \name Constructor/Destructor
   * @{ */

 private:
  /** \brief Default copy constructor is disabled.*/
  ANN_Graph_with_density(const ANN_Graph_with_density& ref) = delete;
  /** \brief Default operator= is disabled.*/
  ANN_Graph_with_density& operator=(const ANN_Graph_with_density& ref) = delete;

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
  ANN_Graph_with_density(const std::string& off_file_name, double mu_)
      : Graph_ANN<Iterator>(off_file_name, mu_) { }

  /** \brief Destructor. */
  ~ANN_Graph_with_density() { }

  /** @} */ // end constructor/destructor

  //----------------------------------
  // get number of neighbors
  // for density estimation
  //----------------------------------

  /**
   *****************************************************************************************
   *  @brief      Get number of neighbors in a given radius.
   *
   *  @note       Rips parameter set on ANN_Graph_with_density constructor is not used. 
   * 
   *  @param[in]  queryPoint Iterator on the point from which the function finds the neighbors.
   *  @param[in]  radius     radius of rips parameter.
   *
   *  @return     Number of neighbors.
   ****************************************************************************************/
/*    int get_num_neighbors(Iterator queryPoint, double radius) {
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
   */
};

#endif  // SRC_TOMATO_INCLUDE_GUDHI_TOMATO_GRAPH_ANN__H_
