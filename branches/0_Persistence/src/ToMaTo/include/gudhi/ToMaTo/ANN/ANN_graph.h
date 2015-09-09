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

#ifndef SRC_TOMATO_INCLUDE_GUDHI_TOMATO_ANN_ANN_GRAPH_H_
#define SRC_TOMATO_INCLUDE_GUDHI_TOMATO_ANN_ANN_GRAPH_H_

#include <ANN/ANN.h>

#include <vector>

#include "gudhi/ToMaTo/Graph.h"
#include "gudhi/ToMaTo/Cluster.h"
#include "gudhi/ToMaTo/ANN/ANN_point.h"

//--------------------------
// when the vertex class
// has a *double coordinate
// so we can use ANN
// assumes a coord member
//--------------------------

/** \brief ANN_graph is an implementation of the Graph concept using ANN points.
 * Uses annkFRSearch on a ANNkd_tree to get rips neighbors.
 *  @param  ANNkd_tree*            <I>kd:</I>                   Persistent data structure for ANN kd tree.
 *  @param  int                    <I>dim:</I>                  Dimension of points in the data structure.
 *  @param  double                 <I>squared_radius:</I>       Rips parameter of the data structure.
 */
template <class Vertex>
class ANN_graph : public Graph< Vertex > {
 public:
  typedef typename std::vector< Vertex >::iterator Iterator;

 protected:
  //----------------------
  // persistent data structure
  // for kd tree
  //----------------------
  ANNkd_tree *kd;
  int dim;
  double squared_radius;

 private:
  // comparison function object for vector indices


  /** \name Constructor/Destructor
   * @{ */

 private:
  /** \brief Default copy constructor is disabled.*/
  ANN_graph(const ANN_graph& ref) = delete;
  /** \brief Default operator= is disabled.*/
  ANN_graph& operator=(const ANN_graph& ref) = delete;

 public:
  /**
   *****************************************************************************************
   *  @brief      Constructs the tree structure.
   *
   *  @details    Default constructor
   * 
   *  @return     None
   ****************************************************************************************/
  ANN_graph()
      : kd(nullptr),
      dim(-1),
      squared_radius(-1.0) { }

  /** \brief Destructor; deallocates the kd tree structure from ANN. */
  virtual ~ANN_graph() {
    delete kd;
  }

  /** @} */  // end constructor/destructor

  /**
   *****************************************************************************************
   *  @brief      Sets the graph structure dimension.
   *
   *  @param[in] dim_   Dimension of the graph structure.
   *
   *  @return     None
   ****************************************************************************************/
  void set_dimension(const int dim_) {
    dim = dim_;
  }

  /**
   *****************************************************************************************
   *  @brief      Returns the graph structure dimension
   *
   *  @return     Dimension of points.
   ****************************************************************************************/
  const int dimension() const {
    return dim;
  }

  /**
   *****************************************************************************************
   *  @brief      Sets the squared radius for ANN neighbor search function.
   *
   *  @param[in]  squared_radius_ The persistence threshold.
   *
   *  @return     None
   ****************************************************************************************/
  void set_sqrad(const double squared_radius_) {
    squared_radius = squared_radius_;
  }

  /**
   *****************************************************************************************
   *  @brief      Constructs the ANN kd tree structure
   *
   *  @details    Constructs an ANN kd tree structure from the vector of vertex.
   * 
   *  @warning    The graph structure dimension must be set before
   *  @warning    Former ANN kd tree structure is deleted
   *
   *  @return     None.
   ****************************************************************************************/
  void construct_ANN_tree() {
    // Assert on default constructed value
    assert(dim != -1);
    int i = 0;
    int num_points = Graph< Vertex >::point_cloud.size();
    std::cout << "initialize - dim=" << dim << " - num_points=" << num_points << std::endl;
    //------------------
    // memory allocation
    //------------------
    double** data = new double*[num_points];
    for (Iterator it = Graph< Vertex >::point_cloud.begin(); it != Graph< Vertex >::point_cloud.end(); it++) {
      data[i++] = it->geometry.coord;
    }
    if (kd != nullptr) {
      // former kd tree removal
      delete kd;
    }
    kd = new ANNkd_tree(data, num_points, dim);
    // TODO(VR): kd deletion does not delete the array of double* - memory leak to be fixed
  }

  void reserve(const unsigned int size) {
    Graph< Vertex >::point_cloud.reserve(size);
  }

  void emplace_back(Vertex point) {
    Graph< Vertex >::point_cloud.emplace_back(point);
  }

  void push_back(Vertex point) {
    Graph< Vertex >::point_cloud.push_back(point);
  }

  //----------------------------------
  // get rips neighbors, this is what
  // algorithm calls
  //----------------------------------

  /**
   *****************************************************************************************
   *  @brief         Get rips neighbors
   *
   *  @details       Sets a vector of neighbors points from a query point.
   *  @note          Neighbors are found from rips. Rips parameter is set on ANN_graph constructor. 
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
    int nb_neighb = kd->annkFRSearch(queryPoint->geometry.coord, squared_radius, 0);

    //---------------------------
    // memory efficient version
    //---------------------------
    // allocate mem each time
    //---------------------------
    int *neighb = new int[nb_neighb];
    memset(neighb, 0, sizeof (int)*nb_neighb);
    double *ndist = new double[nb_neighb];
    memset(ndist, 0, sizeof (double)*nb_neighb);

    assert(kd->annkFRSearch(queryPoint->geometry.coord, squared_radius, nb_neighb, neighb, ndist) == nb_neighb);
    for (int i = 0; i < nb_neighb; i++) {
      Iterator nit = Graph< Vertex >::point_cloud.begin() + neighb[i];
      if (nit != queryPoint) {
        out.push_back(nit);
      }
    }

    delete[] ndist;
    delete[] neighb;
  }

  virtual void compute_persistence() {
    std::cout << "ANN_graph compute_persistence" << std::endl;
    Graph< Vertex >::sort_and_permute();
    construct_ANN_tree();
    Graph< Vertex >::compute_persistence();
  }
};

#endif  // SRC_TOMATO_INCLUDE_GUDHI_TOMATO_ANN_ANN_GRAPH_H_
