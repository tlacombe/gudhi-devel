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


#ifndef _CLUSTER_H
#define _CLUSTER_H

#include <iostream>
#include <fstream>
#include <cmath>
#include <ctime>
#include <algorithm>
#include <map>
#include <set>
#include <vector>  // for vector

#include <gudhi/Cluster_interval.h>

namespace Gudhi {

namespace cluster {

/**
 * \brief Cluster function class.
 *
 * \details
 * The cluster function can be applied on a Cluster_neigborhood_graph data structure, please refer to its concept for
 * more details.
 *
 */
template<class Cluster_neigborhood_graph>
class Cluster {
 public:
  typedef typename Cluster_neigborhood_graph::Neighborhood_iterator Iterator;

 private:
  /** \brief Pointer copy on the Cluster_neigborhood_graph data structure.*/
  Cluster_neigborhood_graph* ngbh_graph_;
  /** \brief Vector of Cluster_interval for 0-persistence computation and output (refer to output_intervals).*/
  std::vector<Cluster_interval> interval_vector_;
  //--------------------------------
  // since we are storing it as a
  // vector, it is equivalent to
  // storing an iterator
  // for faster searching
  //--------------------------------
  /** \brief Map of Cluster_neigborhood_graph iterator to its cluster.*/
  std::map<Iterator, int> generator_;
  /** \brief Persistence threshold.*/
  double tau_;

  // Private default constructor for not to be used.
  Cluster() { }

 public:
  /** \brief Cluster constructor from a Cluster_neigborhood_graph data structure and persistence threshold value.
   * 0-persistence is computed and each element of the Cluster_neigborhood_graph is clustered.
   *
   * @param[in] n_graph Cluster_neigborhood_graph data structure.
   * @param[in] tau Persistence threshold.
   */
  Cluster(Cluster_neigborhood_graph& n_graph, double tau) : ngbh_graph_(&n_graph),
      tau_(tau) {
    compute_persistence();
    attach_to_clusterheads();
  }

  /** \brief Returns the number of clusters. */
  size_t get_nb_clusters() const {
    std::set<Iterator> peaks;

    for (Iterator it = ngbh_graph_->get_start(); it != ngbh_graph_->get_end(); it++) {
      Iterator sink = find_sink(it);
      if (ngbh_graph_->get_func(sink) >= tau_)
        peaks.insert(sink);
    }
    return peaks.size();
  }

  /** \brief Returns the cluster number of an Iterator. */
  int get_cluster(Iterator x) {
    return generator_[find_sink(x)];
  }

  /** \brief Output 0-persistence in a given standard output stream.
   * @param[out] out Standard output stream.
   */
  void output_intervals(std::ostream &out) const {
    for (auto& it : interval_vector_) {
      out << it.get_birth() << " ";
      (it.is_infinite()) ? out << "-inf" : out << it.get_death();
      out << std::endl;
    }
  }

  /** \brief Output graph elements cluster value in a given standard output stream.
   * @param[out] out Standard output stream.
   * \note Input array is an array of iterators on the real point cloud array (to cope with permutation from initial 
   * sort)
   */
  template <class IIterator> void output_clusters(std::ostream &out, IIterator start, IIterator finish) const {
    // run through and create map from prominent peaks
    std::map<Iterator, int> cluster_ids;

    cluster_ids.clear();

    int i = 1;
    //-------------
    // number nodes
    //-------------
    for (IIterator it = start; it != finish; it++) {
      // assert(find_sink(*it) == ngbh_graph_->get_sink(*it));
      typename std::map<Iterator, int>::iterator mit =
          cluster_ids.find(find_sink(*it));
      if (mit == cluster_ids.end()) {
        //--------------
        // check peak height
        //--------------
        if (ngbh_graph_->get_func(find_sink(*it)) >= tau_)
          cluster_ids[find_sink(*it)] = i++;
      }
    }

    for (IIterator it = start; it != finish; it++) {
      //--------------
      // check peak height
      // if it is ok, then
      // output cluster number
      // otherwise output NaN
      //--------------
      if (ngbh_graph_->get_func(find_sink(*it)) >= tau_)
        out << cluster_ids[find_sink(*it)];
      else
        out << "NaN";
      out << std::endl;
    }
  }

  /** \brief Output graph elements and its cluster value in the COFF format (graph elements are the coordinates of the
   * COFF points, and cluster value is represented by its color value) in a given standard output stream.
   * @param[out] out Standard output stream.
   * \warning works only in 2 and 3 dim and outputs height as function value in 2 dim
   */
  void output_clusters_coff(std::ostream &out) const {
    // run through and create map from prominent peaks
    std::map<Iterator, double*> colors;
    std::set<Iterator> peaks;

    typename std::set<Iterator>::iterator lit;

    colors.clear();

    int num_nodes = 0;
    //------------
    // number nodes
    //------------
    for (Iterator it = ngbh_graph_->get_start(); it != ngbh_graph_->get_end(); it++) {
      peaks.insert(find_sink(it));
      num_nodes++;
    }

    //--------------------------------
    // create color scheme
    //--------------------------------
    std::cout << "Number of clusters: "
        << assign_colors(colors, peaks) << std::endl;
    out << "COFF" << std::endl << num_nodes << " " << num_nodes << " 0" << std::endl;

    for (Iterator it = ngbh_graph_->get_start(); it != ngbh_graph_->get_end(); it++) {
      out << ngbh_graph_->get_xyz(it);

      Iterator sink = find_sink(it);
      out << colors[sink][0] << " "
          << colors[sink][1] << " "
          << colors[sink][2] << " "
          << "0" << std::endl;
    }
    for (int i = 0; i < num_nodes; i++) {
      out << "1 " << i << std::endl;
    }
    for (auto& color : colors) {
      delete[] color.second;
    }
  }

 private:
  //-----------------------------
  // Follow linked list to the end
  //-----------------------------

  Iterator find_sink(Iterator in) const {
    assert(in != Iterator());
    while (ngbh_graph_->get_sink(in) != in) {
      in = ngbh_graph_->get_sink(in);
      assert(in != Iterator());
    }
    return in;
  }

  //-----------------------------
  // collapse union-find data structure
  //-----------------------------

  void attach_to_clusterheads() {
    for (Iterator it = ngbh_graph_->get_start(); it != ngbh_graph_->get_end(); it++) {
      ngbh_graph_->set_sink(it, find_sink(it));
    }
  }

  //-----------------------------
  // main algorithm
  //-----------------------------

  void compute_persistence() {
    //=================================
    // assume vertex container is
    // presorted by function value
    //=================================
    for (Iterator vit = ngbh_graph_->get_start(); vit != ngbh_graph_->get_end(); vit++) {
      std::vector<Iterator> adjacent_nodes;
      //-----------------------------------
      // get adjacent_nodes neighbors for each iterator
      //-----------------------------------
      ngbh_graph_->get_neighbors(vit, adjacent_nodes);

      //-----------------------------------
      // find gradient, if it exists
      //-----------------------------------
      Iterator gradient = vit;
      for (auto& neighb : adjacent_nodes) {
        assert(neighb != vit);
        if (neighb < gradient)
          gradient = neighb;
      }


      //------------------------------
      // if no gradient, then declare
      // vit a peak
      //------------------------------
      if (gradient == vit) {
        //-----------------------
        // set the sink to itself
        //-----------------------
        ngbh_graph_->set_sink(vit, vit);

        //-----------------------
        // check to make sure it
        // worked
        //-----------------------
        assert(find_sink(vit) == vit);

        //-----------------------
        // create new cluster
        //-----------------------
        new_cluster(vit);
      } else {
        //----------------------------------
        // gradient has been found
        // attach vit to gradient's cluster
        //----------------------------------
        ngbh_graph_->set_sink(vit, find_sink(gradient));

        //----------------------------------
        // check that vit is right
        // below the cluster root
        // (important invariant in the
        // following)
        //----------------------------------
        assert(ngbh_graph_->get_sink(vit) == find_sink(vit));

        //------------------------------
        // Go through the neighbors again
        // to see if their clusters
        // can be merged with vit's
        //------------------------------
        for (auto& neighb : adjacent_nodes) {
          //----------------------
          // only consider older
          // neighbors
          //----------------------
          if (vit < neighb)
            continue;

          //----------------------
          // find sink of neighbor
          //----------------------
          Iterator sink = find_sink(neighb);

          //----------------------------------
          // check that vit is still right
          // below its cluster's root
          // (cf invariant)
          //----------------------------------
          assert(ngbh_graph_->get_sink(vit) == find_sink(vit));

          //----------------------------
          // no need to do anything
          // if neighb's sink is the
          // same as vit's
          //----------------------------
          if (sink == ngbh_graph_->get_sink(vit))
            continue;

          //----------------------------
          // check if merge conditions
          // are met
          //----------------------------
          if (merge(vit, sink)) {
            //-------------------------
            // If so, then merge cluster
            // with lower peak into
            // cluster with higher peak
            //-------------------------
            if (ngbh_graph_->get_func(ngbh_graph_->get_sink(vit)) >= ngbh_graph_->get_func(sink)) {
              ngbh_graph_->set_sink(sink, ngbh_graph_->get_sink(vit));
            } else {
              ngbh_graph_->set_sink(ngbh_graph_->get_sink(vit), sink);

              //-------------------------
              // compress path from vit
              // to root on the fly
              // (cf invariant)
              //-------------------------
              ngbh_graph_->set_sink(vit, sink);
            }
          }
        }
      }
    }
  }

  //----------------------
  // create a new interval
  //----------------------

  void new_cluster(Iterator x) {
    generator_[x] = interval_vector_.size();
    interval_vector_.push_back(Cluster_interval(ngbh_graph_->get_func(x)));
  }

  //----------------------
  // Merge two intervals
  // note this will only output
  // correctly if persistence
  // threshold is set to infinity
  //----------------------

  bool merge(Iterator x, Iterator y) {
    // note: by hypothesis, func(y) >= func(x)
    assert(ngbh_graph_->get_func(y) >= ngbh_graph_->get_func(x));

    //---------------------------------
    // test prominences of both clusters
    // assumptions:
    //   - y is its cluster's root
    //   - x is attached to its root directly
    //---------------------------------
    if (std::min(ngbh_graph_->get_func(ngbh_graph_->get_sink(x)),
                 ngbh_graph_->get_func(y)) < ngbh_graph_->get_func(x) + tau_) {
      //---------------------------------
      // kill younger interval
      //---------------------------------
      int i = generator_[ngbh_graph_->get_sink(x)];
      int j = generator_[y];
      if (ngbh_graph_->get_func(y) <= ngbh_graph_->get_func(ngbh_graph_->get_sink(x))) {
        assert(interval_vector_[j].is_infinite());
        interval_vector_[j].close(ngbh_graph_->get_func(x));
      } else {
        assert(interval_vector_[i].is_infinite());
        interval_vector_[i].close(ngbh_graph_->get_func(x));
      }
      return true;
    }
    return false;
  }

  //---------------------------------------------
  // This is heuristic to assign colors
  // for the clusters for outputting to an OFF
  // file. Also filters out those clusters whose peak
  // is lower than the threshold tau_
  //---------------------------------------------

  int assign_colors(std::map<Iterator, double*>& cluster_colors,
                    std::set<Iterator> &clusters) const {
    int nb_clusters = clusters.size();

    // reduce nb_clusters to the actual number of clusters
    for (auto& cit : clusters)
      // check if cluster does not appear below tau_
      if (ngbh_graph_->get_func(find_sink(cit)) < tau_)
        nb_clusters--;

    //---------------------------------
    // assign distinct colors to clusters (in H,S,V)
    // first, create array of colors;
    //---------------------------------
    std::vector<double*> colors(nb_clusters + 1, nullptr);
    //---------------------------------
    // just choose circularly
    //---------------------------------
    for (int i = 0; i < nb_clusters; ++i) {
      colors[i] = new double[3];
      double theta = static_cast<double>(2 * i * M_PI) / static_cast<double>(nb_clusters);
      if (theta >= 0 && theta < M_PI / 3) {
        colors[i][0] = 1;
        colors[i][1] = 3 * theta / M_PI;
        colors[i][2] = 0;
      } else if (theta >= M_PI / 3 && theta < 2 * M_PI / 3) {
        colors[i][0] = 1 - 3 * (theta - M_PI / 3) / M_PI;
        colors[i][1] = 1;
        colors[i][2] = 0;
      } else if (theta >= 2 * M_PI / 3 && theta < M_PI) {
        colors[i][0] = 0;
        colors[i][1] = 1;
        colors[i][2] = 3 * (theta - 2 * M_PI / 3) / M_PI;
      } else if (theta >= M_PI && theta < 4 * M_PI / 3) {
        colors[i][0] = 0;
        colors[i][1] = 1 - 3 * (theta - M_PI) / M_PI;
        colors[i][2] = 1;
      } else if (theta >= 4 * M_PI / 3 && theta < 5 * M_PI / 3) {
        colors[i][0] = 3 * (theta - 4 * M_PI / 3) / M_PI;
        colors[i][1] = 0;
        colors[i][2] = 1;
      } else {
        // theta >= 5*M_PI/3 && theta < 2*M_PI
        colors[i][0] = 1;
        colors[i][1] = 0;
        colors[i][2] = 1 - 3 * (theta - 5 * M_PI / 3) / M_PI;
      }
    }
    // force black at index nb_clusters
    colors[nb_clusters] = new double[3];
    colors[nb_clusters][0] = 0;
    colors[nb_clusters][1] = 0;
    colors[nb_clusters][2] = 0;

    //---------------------------------
    // compute random permutation
    //---------------------------------
    int *perm = new int[nb_clusters];
    //---------------------------------
    // first, create identity
    //---------------------------------
    for (int i = 0; i < nb_clusters; ++i)
      perm[i] = i;
    //---------------------------------
    // then, permute it
    //---------------------------------
    srand((unsigned) time(NULL));
    for (int i = 0; i < nb_clusters - 1; ++i) {
      int tmp = perm[i];
      int j = i + static_cast<int>(static_cast<double>(rand()) * (nb_clusters - i) / RAND_MAX);
      perm[i] = perm[j];
      perm[j] = tmp;
    }

    //---------------------------------
    // now assign color entries to the cluster centers
    //---------------------------------
    int k = 0;
    int res = 0;
    for (auto& cit : clusters)
      // check if cluster does not appear below tau_
      if (ngbh_graph_->get_func(find_sink(cit)) >= tau_) {
        cluster_colors[cit] = new double[3];
        std::copy(colors[perm[k]], colors[perm[k]] + 3, cluster_colors[cit]);
        k++;
        res++;
      } else {
        // Black color for non clustered element of graph
        cluster_colors[cit] = new double[3];
        std::copy(colors[nb_clusters], colors[nb_clusters] + 3, cluster_colors[cit]);
      }

    delete[] perm;
    // Colors have been copied, we can delete them.
    for (auto& color : colors) {
      delete[] color;
    }
    return res;
  }
};

}  // namespace cluster

}  // namespace Gudhi

#endif  // _CLUSTER_H
