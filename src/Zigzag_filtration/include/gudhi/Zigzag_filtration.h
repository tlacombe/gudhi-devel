
// #include "tbb/tbb.h"
#include <iostream>
#include <fstream>

/** Given a set of points p_1, ... , p_n ordered by their insertion order in the 
  * filtration, 
  * computes the edge-filtration corresponding to the oscillating Rips zigzag 
  * filtration of the set of points, i.e.,
  * ... <- R({p_0, ... , p_{i}}, nu * eps_i) ->
  *                   R({p_0, ... , p_i, p_{i+1}}, mu * eps_i) <- 
  *                              R({p_0, ... , p_i, p_{i+1}}, nu * eps_{i+1}) -> ...
  * where 0 < nu <= mu, and eps_i is defined as the sparsity of the point cloud 
  * {p_1, ... , p_i}, i.e., the shortest distance between two points in the set. 
  * This is a decreasing sequence of numbers.
  *
  * The function computes the eps_i in filtration_value[i] = eps_i, with eps_0 = 
  * infinity. A simplex appearing in the inclusion  
  * R({p_0, ... , p_{i}}, nu * eps_i) -> R({p_1, ... , p_i, p_{i+1}}, mu * eps_i)
  * is given filtration value eps_i.
  *
  * filtration_values must be empty, receives the filtration values.
  * edge_filtration must be empty, receives the edge filtration.
  */
template<typename Point_container,
         typename Distance, //furnish()
         typename FiltrationValue,
         typename Edge_t >
void points_to_edge_filtration(Point_container const        &points,
                               Distance                      distance,
                               double                        nu,
                               double                        mu,
                               std::vector<FiltrationValue> &filtration_values,
                               std::vector<Edge_t>          &edge_filtration )
{
  //computes the eps_i = min_{p,q \in P_i} d(p,q) naively, in parallel
  size_t n = points.size();
  filtration_values.resize(n);
  filtration_values[0] = std::numeric_limits<double>::infinity();//eps_0
  for(size_t i = 1; i < n; ++i) {//truly parallelisable
    double dist = std::numeric_limits<double>::infinity();
    for(size_t j = 0; j < i; ++j) { //d(P_{i-1}, p_i)
      auto curr_dist = distance(points[i],points[j]); 
      if(dist > curr_dist) { dist = curr_dist; } //maintain shortest distance
    }
    filtration_values[i] = dist; //d(P_{i-1}, p_i)
  }
  //turn filtration_value[i] into sparsity of {p_0, ... , p_i}
  for(size_t i = 1; i < n; ++i) {
    if(filtration_values[i] > filtration_values[i-1]) //make decreasing
    {  filtration_values[i] = filtration_values[i-1];  }
  }
  //initialise R({p_0}, \nu * eps_0)
  edge_filtration.emplace_back(0, 0, filtration_values[0], true);//add p_0
  for(size_t i = 0; i < n-1; ++i) {//all ascending arrows eps_i
    //R({p_0, ... , p_i}, nu * eps_i) -> R({p_1, ... , p_i}, mu * eps_i)   radius   
    for(size_t j = 1; j <= i; ++j) {//nu eps_i < length(p_j, p_k) <= mu eps_i
      for(size_t k = 0; k < j; ++k) {
        if(distance(points[j],points[k]) <= mu * filtration_values[i] && 
           distance(points[j],points[k]) > nu * filtration_values[i]) {
          edge_filtration.emplace_back(k, j, filtration_values[i], true);//edge k,j 
        }
      }
    }
    //R({p_0, ... , p_i}, mu * eps_i) -> R({p_1, ... , p_i, p_i+1}, mu * eps_i)
    edge_filtration.emplace_back(i+1, i+1, filtration_values[i], true);//add p_{i+1}
    for(size_t j = 0; j < i+1; ++j) {//edges (p_{i+1}, p_j) of length <= mu * eps_i 
      if(distance(points[j],points[i+1]) <= mu * filtration_values[i]) {
        edge_filtration.emplace_back(j, i+1, filtration_values[i], true);//edge 
      }
    }
    //R({p_0, ... , p_{i+1}}, mu * eps_i) <- R({p_0, ... , p_{i+1}}, nu * eps_{i+1})
    // std::cout << "Remove edges going down eps_" << i << "to eps_" << i+1 <<" \n";
    for(size_t j = 1; j <= i+1; ++j) {//nu eps_i+1 < length(p_j, p_k) <= mu eps_i
      for(size_t k = 0; k < j; ++k) {
        auto dist = distance(points[j],points[k]);
        if(dist <= mu *filtration_values[i] && dist > nu *filtration_values[i+1]) {
          // std::cout << "  " << j << " " << k << "   " << nu * filtration_values[i+1] << " < *" << dist << "* <= " << mu * filtration_values[i] << "\n";
          edge_filtration.emplace_back(k, j, filtration_values[i+1], false);
        }
      }
    }
    // std::cout << "done.\n";
  }
}

