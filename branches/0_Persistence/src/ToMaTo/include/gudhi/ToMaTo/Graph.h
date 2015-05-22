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

#ifndef SRC_TOMATO_INCLUDE_GUDHI_TOMATO_GRAPH_H_
#define SRC_TOMATO_INCLUDE_GUDHI_TOMATO_GRAPH_H_

#include <vector>
#include <cassert>
#include <string>
#include <algorithm>  // for sort

#include "gudhi/ToMaTo/Cluster.h"
#include "gudhi/ToMaTo/Vertex.h"

template<class V> class less_than {
 protected:
  std::vector<V>& v;

 public:
  less_than(std::vector<V>& v_) : v(v_) { }

  bool operator()(const int a, const int b) const {
    return typename V::less_than()(v[a], v[b]);
  }
};

/** \brief Graph Oracle Concept
 * provides a list of neighbors and return Graph
 */
template <class Vertex>
class Graph {
 public:
  typedef typename std::vector< Vertex >::iterator Iterator;

  /** \brief Graph concept.
   * Graph can compute_persistence
   * Graph must be  get rips neighbors.
   *  @param  std::vector<Vertex>    <I>point_cloud:</I>             Points determining the graph.
   *  @param  std::vector<Iterator>  <I>inverse_permutation:</I>     Inverse permutation as array of iterators on initial point cloud.
   *  @param  Cluster<Iterator>      <I>cluster_data_structure:</I>  Cluster data structure.
   */

 protected:
  std::vector<Vertex> point_cloud;
  std::vector< Iterator > inverse_permutation;
  Cluster< Iterator > cluster_data_structure;


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

  /**
   *****************************************************************************************
   *  @brief     Sets the persistence threshold of the cluster data structure.
   *
   *  @param[in] persistence_threshold_   Dimension of the graph structure.
   *
   *  @return     None
   ****************************************************************************************/
  void set_persistence_threshold(const double persistence_threshold_) {
    cluster_data_structure.tau = persistence_threshold_;
  }

  /**
   *****************************************************************************************
   *  @brief      Returns Sets the persistence threshold of the cluster data structure.
   *
   *  @param      None.
   *
   *  @return     Persistence threshold of the cluster data structure.
   ****************************************************************************************/
  const double persistence_threshold() const {
    return cluster_data_structure.tau;
  }

  virtual void compute_persistence() {
    Iterator gradient;
    Iterator sink;
    std::vector<Iterator> adjacent_nodes;
    typename std::vector<Iterator>::iterator neighb;


    //=================================
    // assume vertex container is
    // pre-sorted by function value
    //=================================
    for (Iterator vit = point_cloud.begin(); vit != point_cloud.end(); vit++) {
      //--------------------------------
      // clear adjacency list
      //--------------------------------
      adjacent_nodes.clear();

      //-----------------------------------
      // get neighbors
      //-----------------------------------
      get_neighbors(vit, adjacent_nodes);

      //-----------------------------------
      // find gradient, if it exists
      //-----------------------------------
      gradient = vit;
      for (neighb = adjacent_nodes.begin(); neighb != adjacent_nodes.end(); neighb++) {
        assert(*neighb != vit);
        if (*neighb < gradient)
          gradient = *neighb;
      }


      //------------------------------
      // if no gradient, then declare
      // vit a peak
      //------------------------------
      if (gradient == vit) {
        //-----------------------
        // set the sink to itself
        //-----------------------
        vit->set_sink(vit);

        //-----------------------------
        // check to make sure it worked
        //-----------------------------
        assert(find_sink(vit) == vit);

        //-----------------------
        // create new cluster
        //-----------------------
        cluster_data_structure.new_cluster(vit);
      } else {
        //----------------------------------
        // gradient has been found
        // attach vit to gradient's cluster
        //----------------------------------
        vit->set_sink(find_sink(gradient));

        //-----------------------------------------------
        // check that vit is right below the cluster root
        // (important invariant in the following)
        //-----------------------------------------------
        assert(vit->get_sink() == find_sink(vit));

        //------------------------------
        // Go through the neighbors again
        // to see if their clusters
        // can be merged with vit's
        //------------------------------
        for (neighb = adjacent_nodes.begin();
             neighb != adjacent_nodes.end();
             neighb++) {
          //----------------------
          // only consider older
          // neighbors
          //----------------------
          if (vit < *neighb)
            continue;

          //----------------------
          // find sink of neighbor
          //----------------------
          sink = find_sink(*neighb);

          //----------------------------------
          // check that vit is still right
          // below its cluster's root
          // (cf invariant)
          //----------------------------------
          assert(vit->get_sink() == find_sink(vit));

          //----------------------------
          // no need to do anything
          // if neighb's sink is the
          // same as vit's
          //----------------------------
          if (sink == vit->get_sink())
            continue;

          //----------------------------
          // check if merge conditions
          // are met
          //----------------------------
          if (cluster_data_structure.merge(vit, sink)) {
            //--------------------------
            // If so, then merge cluster
            // with lower peak into
            // cluster with higher peak
            //--------------------------
            if (vit->get_sink()->func() >= sink->func()) {
              sink->set_sink(vit->get_sink());
            } else {
              vit->get_sink()->set_sink(sink);

              //-------------------------
              // compress path from vit
              // to root on the fly
              // (cf invariant)
              //-------------------------
              vit->set_sink(sink);
            }
          }
        }
      }
    }
    attach_to_clusterheads();
  }

  void sort_and_permute() {
    std::cout << "sort_and_permute" << std::endl;
    int num_points = point_cloud.size();
    // sort point cloud and retrieve permutation (for pretty output)
    std::vector<int> perm;
    perm.reserve(num_points);

    for (int i = 0; i < num_points; i++) {
      perm.push_back(i);
    }

    std::sort(perm.begin(), perm.end(), less_than<Vertex>(point_cloud));
    // store inverse permutation as array of iterators on initial point cloud
    inverse_permutation.reserve(num_points);
    for (int i = 0; i < num_points; i++)
      inverse_permutation.push_back(point_cloud.begin());
    for (int i = 0; i < num_points; i++)
      inverse_permutation[perm[i]] = (point_cloud.begin() + i);
    // operate permutation on initial point cloud
    std::vector<Vertex> pc;
    pc.reserve(num_points);
    for (int i = 0; i < num_points; i++)
      pc.push_back(point_cloud[i]);
    for (int i = 0; i < num_points; i++)
      point_cloud[i] = pc[perm[i]];
  }

  bool output_intervals(const std::string& file_name) {
    std::ofstream out;
    // output barcode
    out.open(file_name);
    if (out.is_open()) {
      cluster_data_structure.output_intervals(out);
      out.close();
      return true;
    }
    return false;
  }

  bool output_clusters_off(const std::string& file_name) {
    std::ofstream out;
    // output clusters in a OFF file
    out.open(file_name);
    if (out.is_open()) {
      cluster_data_structure.output_clusters_coff(out, point_cloud.begin(), point_cloud.end());
      out.close();
      return true;
    }
    return false;
  }

  bool output_clusters(const std::string& file_name) {
    std::ofstream out;
    // output clusters (use permutation to preserve original point order)
    out.open(file_name);
    if (out.is_open()) {
      cluster_data_structure.output_clusters(out, inverse_permutation.begin(), inverse_permutation.end());
      out.close();
      return true;
    }
    return false;
  }

 private:
  Iterator find_sink(Iterator in) {
    assert(in != Iterator());
    while (in->get_sink() != in) {
      in = in->get_sink();
      assert(in != Iterator());
    }
    return in;
  }

  void attach_to_clusterheads() {
    for (Iterator vit = point_cloud.begin(); vit != point_cloud.end(); vit++) {
      vit->set_sink(find_sink(vit));
    }
  }
};

#endif  // SRC_TOMATO_INCLUDE_GUDHI_TOMATO_GRAPH_H_
