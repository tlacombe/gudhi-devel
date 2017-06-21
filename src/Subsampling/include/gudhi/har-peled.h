/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Siargey Kachanovich
 *
 *    Copyright (C) 2017  INRIA Sophia Antipolis-Méditerranée (France)
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

#ifndef LANDMARK_CHOICE_BY_FARTHEST_POINT_H_
#define LANDMARK_CHOICE_BY_FARTHEST_POINT_H_

#include <gudhi/Spatial_tree_data_structure.h>

#include <iterator>
#include <algorithm>  // for sort
#include <vector>
#include <list>
#include <random>
#include <boost/heap/fibonacci_heap.hpp>

namespace Gudhi {

namespace subsampling {
  

  template < typename Point_d,
             typename Heap,
             typename Tree,
             typename Presence_table >
  void update_heap( Point_d &l,
                    unsigned nbL,
                    Heap &heap,
                    Tree &tree,
                    Presence_table &table)
  {
    auto search = tree.query_incremental_ANN(l);
    for (auto w: search) {
      if (table[w.first].first)
        if (w.second < table[w.first].second->second) {
          heap.update(table[w.first].second, w);
        }
    }
  }
  
  /** 
   *  \ingroup witness_complex
   *  \brief Landmark choice strategy by iteratively adding the farthest witness from the
   *  current landmark set as the new landmark. 
   *  \details It chooses nbL landmarks from a random access range `points` and
   *  writes {witness}*{closest landmarks} matrix in `knn`.
   *
   *  The type KNearestNeighbors can be seen as 
   *  Witness_range<Closest_landmark_range<Vertex_handle>>, where
   *  Witness_range and Closest_landmark_range are random access ranges 
   *  
   */

  template < typename Kernel,
             typename Point_range,
             typename OutputIterator>
  void harpeled_mendel(Kernel const &k,
                       Point_range const &input_pts,
                       std::size_t final_size,
                       std::size_t starting_point,
                       OutputIterator output_it)
  {

    typedef typename Kernel::FT FT;
    typedef std::size_t Point_id;
    typedef std::pair<Point_id, FT> Heap_node;
    struct R_max_compare
    {
      bool operator()(const Heap_node &rmh1, const Heap_node &rmh2) const
      {
        return rmh1.second < rmh2.second;
      }
    };
    typedef boost::heap::fibonacci_heap<Heap_node, boost::heap::compare<R_max_compare>> Cluster_heap;
    typedef std::vector<Cluster_heap> Clusters;
    typedef typename Clusters::iterator Max_heap_node;
    struct M_max_compare
    {
      bool operator()(const Max_heap_node &mhn1, const Max_heap_node &mhn2) const
      {
        return mhn1->top().second < mhn2->top().second;
      }
    };
    typedef boost::heap::fibonacci_heap<Max_heap_node, boost::heap::compare<M_max_compare>> Max_heap;
    typedef std::list<Point_id> Friend_list;
    
    typename Kernel::Squared_distance_d sqdist = k.squared_distance_d_object();

    std::size_t nb_points = boost::size(input_pts);
    //assert(nb_points >= final_size);

    std::vector<Point_id> centers, prev_centers, prev_prev_centers;
    centers.reserve(final_size);

    Clusters clusters;
    clusters.reserve(final_size);    

    Max_heap max_heap;

    std::vector<Friend_list> friend_lists;
    friend_lists.reserve(final_size);

    centers.push_back(starting_point);
    *output_it++ = input_pts[starting_point];
    for (std::size_t i = 0; i < nb_points; ++i) {
      prev_centers.push_back(starting_point);
      prev_prev_centers.push_back(starting_point);
    }
    clusters.push_back(Cluster_heap());
    for (std::size_t i = 0; i < nb_points; ++i)
      clusters[0].push(Heap_node(i, sqdist(input_pts[centers[0]], input_pts[i])));
    max_heap.push(clusters.begin());
    friend_lists.push_back(Friend_list());
    FT r = std::numeric_limits<double>::infinity();
    
    // Assumption: Lazy max_heap update. The friends' lists don't pop to max unless they are max.
    for (std::size_t i = 1; i < final_size; ++i) {
      centers.push_back(max_heap.top()->top().first);
      *output_it++ = input_pts[max_heap.top()->top().first];
      r = max_heap.top()->top().second;
      clusters.push_back(Cluster_heap());
      for (std::size_t j = 0; j < nb_points; ++j)
        prev_prev_centers[j] = prev_centers[j];
      std::vector<typename Cluster_heap::iterator> to_erase;
      std::size_t c_k = std::distance(clusters.begin(), max_heap.top());
      for (auto c_it = max_heap.top()->begin(); c_it != max_heap.top()->end(); ++c_it) {
        FT new_dist = sqdist(input_pts[centers[i]], input_pts[c_it->first]);
        if (new_dist < c_it->second) {
          clusters[i].push(Heap_node(c_it->first, new_dist));
          prev_centers[c_it->first] = c_k;
          to_erase.push_back(c_it);
        }
      }
      for (auto c_it = to_erase.begin(); c_it != to_erase.end(); ++c_it) 
        max_heap.top()->erase(Cluster_heap::s_handle_from_iterator(*c_it));
      auto fr_it = friend_lists[c_k].begin();
      while (fr_it != friend_lists[c_k].end()) {
        to_erase.clear();
        for (auto c_it = clusters[*fr_it].begin(); c_it != clusters[*fr_it].end(); ++c_it) {
          FT new_dist = sqdist(input_pts[centers[i]], input_pts[c_it->first]);
          if (new_dist < c_it->second) {
            clusters[i].push(Heap_node(c_it->first, new_dist));
            prev_centers[c_it->first] = *fr_it;
            to_erase.push_back(c_it);
          }
        }
        for (auto c_it = to_erase.begin(); c_it != to_erase.end(); ++c_it) 
          clusters[*fr_it].erase(Cluster_heap::s_handle_from_iterator(*c_it));
        if (sqdist(input_pts[centers[*fr_it]], input_pts[centers[c_k]]) > 64*r)
          friend_lists[c_k].erase(fr_it++);
        else
          fr_it++;
      }
      max_heap.top()->pop();
      friend_lists.push_back(Friend_list());
      fr_it = friend_lists[prev_prev_centers[centers[i]]].begin();
      while (fr_it != friend_lists[prev_prev_centers[centers[i]]].end()) {
        if (sqdist(input_pts[centers[*fr_it]], input_pts[centers[i]]) <= 16*r)
          friend_lists[i].push_back(*fr_it);
        if (sqdist(input_pts[centers[*fr_it]], input_pts[prev_prev_centers[centers[i]]]) > 64*r)
          friend_lists[prev_prev_centers[centers[i]]].erase(fr_it++);
        else
          fr_it++;
      }
    }
    
    // typedef boost::heap::fibonacci_heap<Heap_node, boost::heap::compare<R_max_compare>> Heap;
    // typedef Spatial_tree_data_structure<Kernel, Point_container> Tree;
    // typedef std::vector< std::pair<bool, Heap_node*> > Presence_table;
    
  //   Tree tree(points);
  //   Heap heap;
  //   Presence_table table(points.size());
  //   for (auto p: table)
  //     std::cout << p.first << "\n";
  //   int number_landmarks = 0; // number of treated landmarks

  //   double curr_max_dist = 0;                                      // used for defining the furhest point from L
  //   const double infty = std::numeric_limits<double>::infinity();  // infinity (see next entry)
  //   std::vector< double > dist_to_L(points.size(), infty);         // vector of current distances to L from points
    
  //   // Choose randomly the first landmark 
  //   std::random_device rd;
  //   std::mt19937 gen(rd());
  //   std::uniform_int_distribution<> dis(1, 6);
  //   int curr_landmark = dis(gen);
    
  //   do {
  //     *output_landmarks++ = points[curr_landmark];
  //     std::cout << curr_landmark << "\n";
  //     number_landmarks++;
  //   }
  //   while (number_landmarks < nbL);
  // }

    // int nb_points = boost::size(points);
    // assert(nb_points >= nbL);

    // int current_number_of_landmarks = 0;  // counter for landmarks
    // double curr_max_dist = 0;  // used for defining the furhest point from L
    // const double infty = std::numeric_limits<double>::infinity();  // infinity (see next entry)
    // std::vector< double > dist_to_L(nb_points, infty);  // vector of current distances to L from points

    // // Choose randomly the first landmark 
    // std::random_device rd;
    // std::mt19937 gen(rd());
    // std::uniform_int_distribution<> dis(1, 6);
    // int curr_max_w = dis(gen);

    
    // for (current_number_of_landmarks = 0; current_number_of_landmarks != nbL; current_number_of_landmarks++) {
    //   // curr_max_w at this point is the next landmark
    //   *output_it++ = points[curr_max_w];
    //   std::cout << curr_max_w << "\n";
    //   unsigned i = 0;
    //   for (auto& p : points) {
    //     double curr_dist = sqdist(p, *(std::begin(points) + curr_max_w));
    //     if (curr_dist < dist_to_L[i])
    //       dist_to_L[i] = curr_dist;
    //     ++i;
    //   }
    //   // choose the next curr_max_w
    //   curr_max_dist = 0;
    //   for (i = 0; i < dist_to_L.size(); i++)
    //     if (dist_to_L[i] > curr_max_dist) {
    //       curr_max_dist = dist_to_L[i];
    //       curr_max_w = i;
    //     }
    // }
  }

}  // namespace subsampling  
  
}  // namespace Gudhi

#endif  // LANDMARK_CHOICE_BY_FARTHEST_POINT_H_
