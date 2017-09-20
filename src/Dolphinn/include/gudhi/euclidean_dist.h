#ifndef DOLPHINN_EUCLIDEAN_DIST_H
#define DOLPHINN_EUCLIDEAN_DIST_H

#include <vector>
#include <random>
#include <math.h>
#include <utility> //for pair



namespace Gudhi {
namespace dolphinn {

/** \brief Euclidean distance squared.
 *
 * @param p1       	first point
 * @param p2   			second point
 * @return          the Euclidean distance between p1 and p2
 */
template<typename iterator>
double squared_Eucl_distance(const iterator p1, const iterator p2)
{
	double res=0;
	double tmp;
	for(size_t i=0;i<(*p1).size();++i){
		tmp = (*p1)[i] - (*p2)[i];
		res += tmp*tmp;
	}
  return res;
}

/** \brief Report a point's index (if any) that has Euclidean distance
 * less or equal than a given radius.
 *
 * @param pointset        vector of all points
 * @param points_idxs     indices of candidate points
 * @param D               dimension of points
 * @param query_point     vector containing only the coordinates of the query point
 * @param squared_radius  square value of given radius
 * @param threshold       max number of points to check
 * @return                the index of the point. -1 if not found.
 */
template <typename iterator>
int Euclidean_distance_within_radius(iterator pointset, const std::vector<int>& points_idxs,
 const int D, iterator query_point, const double squared_radius, const int threshold)
{
  const int size = points_idxs.size();
  for(int i = 0; i < threshold && i < size; ++i)
  {
    if(squared_Eucl_distance(query_point, pointset + points_idxs[i]) <= squared_radius)
      return points_idxs[i];
  }
  return -1;
}

/** \brief Collects all the indices of the points that have Euclidean distance
 * less or equal than a given radius.
 *
 * @param pointset         vector of all points
 * @param points_idxs      indices of candidate points
 * @param D                dimension of points
 * @param query_point      vector containing only the coordinates of the query point
 * @param squared_radius   square value of given radius
 * @param threshold        max number of points to check
 * @param answer_point_idx indices of the neighbors already found
 */
template <typename iterator>
void Euclidean_distance_within_radius_all(iterator pointset, const std::vector<int>& points_idxs,
 const int D, iterator query_point, const double squared_radius, const int threshold, std::vector<int>& answer_point_idx)
{
  const int size = points_idxs.size();
  for(int i = 0; i < threshold && i < size; ++i)
  {
    if(squared_Eucl_distance(query_point, pointset + points_idxs[i]) <= squared_radius)
      answer_point_idx.push_back(points_idxs[i]);
  }
}

/** \brief Report M Nearest Neighbors' indices, if something better than the current NN is found.
 *
 * @param pointset              vector of all points
 * @param points_idxs           indices of candidate points
 * @param D                     dimension of points
 * @param M
 * @param query_point           vector containing only the coordinates of the query point
 * @param answer_point_idx_dist current best NN points. Will be updated if a point closer to the query is found.
 * @param threshold             max number of points to check
 */
template <typename iterator>
void find_M_Nearest_Neighbor_indices(iterator pointset, const std::vector<int>& points_idxs,
 const int D, const int M, iterator query_point, std::vector<std::pair<int, double>>& answer_point_idx_dist, const int threshold)
{
  const int size = points_idxs.size();
  for(int i = 0; i < threshold && i < size; ++i)
  {
  	double current_dist = squared_Eucl_distance(query_point, pointset + points_idxs[i]);
  	if(current_dist<answer_point_idx_dist[M-1].second){
  		int j=M-2;
  		while(j>=0 && current_dist<answer_point_idx_dist[j].second){
  			answer_point_idx_dist[j+1].second = answer_point_idx_dist[j].second;
  			answer_point_idx_dist[j+1].first = answer_point_idx_dist[j].first;
  			--j;
  		}
  		answer_point_idx_dist[j+1].second = current_dist;
  		answer_point_idx_dist[j+1].first = points_idxs[i];
  	}
  }
}
}
}

#endif /*DOLPHINN_EUCLIDEAN_DIST_H*/
