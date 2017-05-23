#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE "dolphinn"
#include <boost/test/unit_test.hpp>

#include <gudhi/Dolphinn.h>
#include <gudhi/Points_off_io.h>
#include <vector>
#include <random>
#include <cmath>
#include <utility> //for pair

#define T float
#define bitT char
#define Point std::vector<T>
#define N 10000
#define K 4
#define D 10
#define Q 5

using Dolphinn = Gudhi::dolphinn::Dolphinn<T, bitT>;

BOOST_AUTO_TEST_CASE(hypercube_building_lines) {
  std::vector<Point> pointset;
  std::default_random_engine generator(std::chrono::system_clock::now().time_since_epoch().count());
  std::normal_distribution<float> distribution(0.0,1.0/std::sqrt((float)D));
  for(size_t i=0;i<N;++i){
  	Point p;
  	for(size_t j=0;j<D;++j){
  		p.push_back(distribution(generator));
  	}
  	pointset.push_back(p);
  }
    	
  Dolphinn dolphi(pointset, N, D, K, 0);
  const Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic> m = dolphi.get_hypercube().get_H()[0].get_m();
  const float* ptr = &pointset[0][0];
  Eigen::Matrix<float, Eigen::Dynamic, 1> v =  m * Eigen::Map<const Eigen::Matrix<float ,Eigen::Dynamic, 1>>(ptr,D,1);
  float s=0;
  for(int i=0;i<D;i++){s+=v[i]*v[i];}
  int num=0;
  for(auto x:dolphi.get_hypercube().get_H()[0].get_hashtable())
  	num += x.second.size();
  
  //no point lost
  BOOST_CHECK(num==10000);
  //no absurd hash
  BOOST_CHECK(s>0);
}


BOOST_AUTO_TEST_CASE(hypercube_building_hyperplanes) {
	std::vector<Point> pointset;
  std::default_random_engine generator(std::chrono::system_clock::now().time_since_epoch().count());
  std::normal_distribution<float> distribution(0.0,1.0/std::sqrt((float)D));
  for(size_t i=0;i<N;++i){
  	Point p;
  	for(size_t j=0;j<D;++j){
  		p.push_back(distribution(generator));
  	}
  	pointset.push_back(p);
  }
  Dolphinn dolphi(pointset, N, D, K, 0.0001);
  int num=0;
  for(auto x:dolphi.get_hypercube().get_H()[K-1].get_hashtable())
  	num += x.second.size();
  
  
  //no point lost
  BOOST_CHECK(num==N);
  //no absurd hash
  for(auto& h:dolphi.get_hypercube().get_H()){
  	float norm = 0;
		for(size_t i=0;i<D;++i){
			norm += h.get_a()[i] * h.get_a()[i];
		}
		BOOST_CHECK(norm != 0);
		BOOST_CHECK(h.get_b()!=0);
		BOOST_CHECK(h.get_r()!=0);
  }
}


BOOST_AUTO_TEST_CASE(knn_query) {
	std::vector<Point> pointset;
  std::default_random_engine generator(std::chrono::system_clock::now().time_since_epoch().count());
  std::normal_distribution<float> distribution(0.0,1.0/std::sqrt((float)D));
  for(size_t i=0;i<N;++i){
  	Point p;
  	for(size_t j=0;j<D;++j){
  		p.push_back(distribution(generator));
  	}
  	pointset.push_back(p);
  }
  Dolphinn dolphi(pointset, N, D, K, 0);
  
  std::vector<Point> queries;
  std::vector<std::vector<std::pair<int, float>>> result;
	for(size_t i=0;i<Q;++i){
  	Point p;
  	for(size_t j=0;j<D;++j){
  		p.push_back(distribution(generator));
  	}
  	std::vector<std::pair<int, float>> dummy;
  	result.push_back(dummy);
  	queries.push_back(p);
  }
	
	dolphi.m_nearest_neighbors_query(queries, Q, 10, N/100, result, 1);
	float mean=0;
	int no_nn =0;
	for(size_t j=0;j<Q;++j){
		for(auto& x:result[j]){
			if(x.first !=-1) {
				float d=0;
				for(size_t i=0;i<D;++i){
					d+=(queries[j][i]-pointset[x.first][i])*(queries[j][i]-pointset[x.first][i]);
				}
				//is the output distance correct
				BOOST_CHECK(d==x.second);
				mean += std::sqrt(x.second);
			} else {++no_nn;}
		}
	}
	//are the ouput points some kind of neighbours (will fail with D too large)
	BOOST_CHECK(mean<std::sqrt(2)*(Q*10-no_nn)-5);
	
	std::vector<std::vector<std::pair<int, float>>> result2;
	std::vector<std::pair<int, float>> dummy;
  result2.push_back(dummy);
	dolphi.m_nearest_neighbors_query(queries, 1, 1, N*2, result2, 1);
	int min_idx = 0;
	float min_dist = 500;
	
	for(size_t i=0;i<N;++i){
		float d=0;
		for(size_t j=0;j<D;++j){
			d+=(queries[0][j]-pointset[i][j])*(queries[0][j]-pointset[i][j]);
		}
		if(d<min_dist) {min_dist=d; min_idx=i;}
	}
		
	//if all the points are checked, is the actual nearest neighbour output 
	BOOST_CHECK(result2[0][0].first==min_idx && result2[0][0].second==min_dist);
}



BOOST_AUTO_TEST_CASE(radius_query) {
	std::vector<Point> pointset;
  std::default_random_engine generator(std::chrono::system_clock::now().time_since_epoch().count());
  std::normal_distribution<float> distribution(0.0,1.0/std::sqrt((float)D));
  for(size_t i=0;i<N;++i){
  	Point p;
  	for(size_t j=0;j<D;++j){
  		p.push_back(distribution(generator));
  	}
  	pointset.push_back(p);
  }
  Dolphinn dolphi(pointset, N, D, K, 0.5);
  
  std::vector<Point> queries;
  std::vector<int> result(Q);
	for(size_t i=0;i<Q;++i){
  	Point p;
  	for(size_t j=0;j<D;++j){
  		p.push_back(distribution(generator));
  	}
  	queries.push_back(p);
  }
	
	dolphi.radius_query(queries, Q, 0.8, N/100, result, 1);
	for(int i=0; i<Q;++i){
		if(result[i]!=-1){
			float d=0;
			for(size_t j=0;j<D;++j){
				d+=(queries[i][j]-pointset[result[i]][j])*(queries[i][j]-pointset[result[i]][j]);
			}
			//is the answer inside the radius?
			BOOST_CHECK(d<=0.8*0.8);
		}
	}
	
	// for find_strings_with_fixed_Hamming_dist_for_radius_query (Stable_hash_function.h)
	std::vector<int> result2(Q);
	Dolphinn dolphi2(pointset, N, D, K+5, 0);
	dolphi2.radius_query(queries, Q, 0.000001, N, result2, 1);
	
}
