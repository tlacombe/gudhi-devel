#ifndef DOLPHINN_H
#define DOLPHINN_H

#include <gudhi/Hypercube.h>
#include <thread>

//#define Point std::vector<T>

namespace Gudhi {
namespace dolphinn {

/**
 * \class Dolphinn Dolphinn.h gudhi/Dolphinn.h
 * \brief Method for approximate neighbor search.
 * 
 * \ingroup dolphinn
 *
 * \details
 *
 * Dolphinn projects all the data on the vertices of an hypercube aiming to send close points to 
 * close vertices (w.r.t. the hamming distance). It offers two types of query: the k-nearest neigbor 
 * search and the range search.
 *
 * \remark When the class is built, the hypercube is immediately built.
 * 
 * 
 *
 */


	template <typename T>
  class Dolphinn
  {
  	
  	public:
  	/** \brief A point in Euclidean space.*/
  	typedef typename std::vector<T> Point;
  	/** \brief type of the coordinates of the points*/
		typedef T data_type;
		private:
  	// See constructor
  	int N;
  	const int D,K;
		double hashing_method;
		std::vector<Point>& pointset;
		Hypercube<Point, T> hypercube;

		public:
		
		/** An accessor to the hypercube for tests purposes */
		Hypercube<Point, T> get_hypercube(){
			return hypercube;
		}
		
		

  	/** \brief Constructor of the class, fills the hypercube.
      *
      * @param pointset    		Set of points
      * @param N           		number of points
      * @param D           		dimension of the points
      * @param K           		dimension of Hypercube (and of the mapped points)
      * @param hashing_method if positive, the parameter of Stable Distribution, if nul, the LSH method used is the hyperplanes.
      *                      Neighbor Search, to adapt to the average distance of the NN, 'r' is the hashing window.
   */
  	
  	Dolphinn(std::vector<Point>& pointset, int N, const int D, const int K, const double hashing_method) : N(pointset.size()), D(D), K(K), hashing_method(hashing_method), pointset(pointset), hypercube(pointset, N, D, K, 1, hashing_method) 
  		{}
  	
  	/** \brief Radius query the Hamming cube, returns the index of only ONE point per query.
      *
      * @param query               vector of queries
      * @param Q                   number of queries
      * @param radius              find a point within r with query
      * @param max_pnts_to_search  number of candidates checked before giving up
      * @param results_idxs        indices of Q points, where Eucl(point[i], query[i]) <= r
    */
  	void radius_query(const std::vector<Point>& query, const int Q, const double radius, const int max_pnts_to_search, std::vector<int>& results_idxs) {
  		if(max_pnts_to_search>N){
  			hypercube.radius_query(query, Q, radius, N, results_idxs);
  		} else {
  			hypercube.radius_query(query, Q, radius, max_pnts_to_search, results_idxs);
  		}
  	}
  	
  	/** \brief Nearest Neighbor query in the Hamming cube.
      *
      * @param query               vector of queries
      * @param Q                   number of queries
      * @param k                   number of neighbors to search 
      * @param max_pnts_to_search  number of candidates checked before giving up
      * @param results_idxs_dists  indices and distances of Q points, where the approximate Nearest Neighbors are stored.
    	*/
  	void k_nearest_neighbors_query(const std::vector<Point>& query, const int Q, const int k, const int max_pnts_to_search, std::vector<std::vector<std::pair<int, double>>>& results_idxs_dists) {
  		if(max_pnts_to_search>N){
  			hypercube.m_nearest_neighbors_query(query, Q, k, N, results_idxs_dists);
  		} else {
  			hypercube.m_nearest_neighbors_query(query, Q, k, max_pnts_to_search, results_idxs_dists);
  		}
  		
  	}
  	
  	
  	
  };
}
}

#endif /* DOLPHINN_H*/
