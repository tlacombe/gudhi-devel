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
#include "gudhi/ToMaTo/ANN/ANN_graph.h"

//--------------------------
// when the vertex class
// has a *double coordinate
// so we can use ANN
// assumes a coord member
//--------------------------

/** \brief ANN_graph_with_density is an implementation of the ANN_Graph.
 */
template <class Vertex>
class ANN_graph_with_density : public ANN_graph<Vertex> {
  typedef typename std::vector< Vertex >::iterator Iterator;
  
  /** \name Constructor/Destructor
   * @{ */

 private:
  /** \brief Default operator= is disabled.*/
  ANN_graph_with_density& operator=(const ANN_graph_with_density& ref) = delete;

 public:
/**
   *****************************************************************************************
   *  @brief      Constructs the tree structure.
   *
   *  @usage      Default constructor
   *  @note       num_points is computed again. 
   * 
   *  @return     None
   ****************************************************************************************/
  ANN_graph_with_density()
      : ANN_graph<Vertex>() { }

  /** \brief Destructor. */
  ~ANN_graph_with_density() { }

  /** @} */ // end constructor/destructor
  
//------------------------------
// distance to measure
//------------------------------
void distance_to_density(int k)
{
  Iterator it;
  double sum;
  double *ndist = new double[k];

  for(it=Graph< Vertex >::point_cloud.begin();it!=Graph< Vertex >::point_cloud.end();it++){
    get_neighbors_dist(it,k,ndist);  
    sum=0;
    for (int i=0; i<k; ++i) {
      sum += ndist[i]*ndist[i];
    }
    it->set_func(sqrt((double)k/(double)sum));
  }  

  delete ndist;
  Graph<Vertex>::sort_and_permute();
}


 private:
  //----------------------------------
  // get number of neighbors
  // for density estimation
  //----------------------------------
  /**
   *****************************************************************************************
   *  @brief      Get number of neighbors in a given radius.
   *
   *  @note       Rips parameter set on ANN_graph_with_density constructor is not used. 
   * 
   *  @param[in]  queryPoint Iterator on the point from which the function finds the neighbors.
   *  @param[in]  radius     radius of rips parameter.
   *
   *  @return     Number of neighbors.
   ****************************************************************************************/
    int get_num_neighbors(Iterator queryPoint, double radius) {
      int nb_neighb = ANN_graph<Vertex>::kd->annkFRSearch(queryPoint->geometry.coord, radius*radius, 0);
      return nb_neighb;
    }

    //----------------------------------
    // get distances to k nearest neighbors
    // for density estimation
    //----------------------------------

    void get_neighbors_dist(Iterator queryPoint, int k, double *ndist) {
      int *neighb = new int[k];
      memset(neighb, 0, sizeof (int)*k);

      ANN_graph<Vertex>::kd->annkSearch(queryPoint->geometry.coord, k, neighb, ndist, 0);
      delete[] neighb;
    }

    //----------------------------------
    // get distances neighbors within r
    // for density estimation
    //----------------------------------  

    int get_neighbors_dist_r(Iterator queryPoint, double radius, double *ndist) {
      int nb_neighb = ANN_graph<Vertex>::kd->annkFRSearch(queryPoint->geometry.coord, radius*radius, 0);

      int *neighb = new int[nb_neighb];
      memset(neighb, 0, sizeof (int)*nb_neighb);
      int test = ANN_graph<Vertex>::kd->annkFRSearch(queryPoint->geometry.coord, radius*radius, nb_neighb, neighb, ndist, 0);

      delete[] neighb;
      return test;
    }
};

#endif  // SRC_TOMATO_INCLUDE_GUDHI_TOMATO_GRAPH_ANN__H_
