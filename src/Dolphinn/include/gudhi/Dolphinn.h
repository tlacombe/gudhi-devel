#include "Hypercube.h"
#include <thread>



namespace Gudhi {
namespace dolphinn {

	template <typename T, typename bitT>
  class Dolphinn
  {
  public:
  	int N,D,K;
		double hashing_method;
		std::vector<T>& pointset;
		Hypercube<T, bitT> hypercube;

  	/** \brief Constructor that only fills the members of the class apart from the hypercube.
      *
      * @param pointset    		- 1D vector of points, emulating a 2D, with N rows and D columns per row.
      * @param N           		- number of points
      * @param D           		- dimension of the points
      * @param K           		- dimension of Hypercube (and of the mapped points)
      * @param hashing_method - if positive, the parameter of Stable Distribution, if nul, the LSH method used is the hyperplanes.
      *                      Neighbor Search, to adapt to the average distance of the NN, 'r' is the hashing window.
   */
  	
  	Dolphinn(std::vector<T>& pointset, const int N, const int D, const int K, const double hashing_method) : N(N), D(D), K(K), hashing_method(hashing_method), pointset(pointset), hypercube(pointset, N, D, K, std::thread::hardware_concurrency(), hashing_method) 
  		{}
  	
  	/*TODO*/
  	void add_data_point() {
  		
  	}
  	
  	
  	
  	void radius_query(const std::vector<T>& query, const int Q, const int radius, const int max_pnts_to_search, std::vector<int>& results_idxs, const int threads_no = std::thread::hardware_concurrency()) {
  		hypercube.radius_query(query, Q, radius, max_pnts_to_search, results_idxs, threads_no);
  	}
  	
  	void m_nearest_neighbors_query(const std::vector<T>& query, const int Q, const int m, const int max_pnts_to_search, std::vector<std::vector<std::pair<int, float>>>& results_idxs_dists, const int threads_no = std::thread::hardware_concurrency()) {
  		hypercube.m_nearest_neighbors_query(query, Q, m, max_pnts_to_search, results_idxs_dists, threads_no);
  	}
  	
  	
  	
  };
}
}
