#include <iostream>
#include "gudhi/Dolphinn.h"
#include <gudhi/Points_off_io.h>
#include <vector>
#include <random>
#include <cmath> // for sqrt
#include <gudhi/IO.h>

#include <ctime>
#include <ratio>
#include <chrono>

#define T float
#define bitT char
#define Point std::vector<T>


using namespace Gudhi;
using namespace std::chrono;

float compare (size_t q, size_t m, const std::vector<std::vector<std::pair<int, float>>> result, const std::vector<std::vector<std::pair<int, float>>> resultn){
	float correct =0.0;
  for(size_t i=0;i<q;++i){
  	int c=0;
  	for(size_t j=0;j<m;j++){
  		if(result[i][c].first == resultn[i][j].first){++correct;++c;}
  	}
  }
  return (correct*100)/(q*m);
}

double query_test(std::vector<std::vector<std::pair<int, float>>>& result, size_t n, size_t d, size_t m, size_t q, std::vector<Point> pointset, std::vector<Point> queries, size_t k, size_t exp){
	dolphinn::Dolphinn<T, bitT> dolphi(pointset, n, d, k, 0);
	
	high_resolution_clock::time_point t1 = high_resolution_clock::now();
	dolphi.m_nearest_neighbors_query(queries, q, m, exp, result, 1);
	high_resolution_clock::time_point t2 = high_resolution_clock::now();
	duration<double> time_span = duration_cast<duration<double>>(t2 - t1);
	
	for (auto& x:result){
  	std::sort(x.begin(),x.end(),[](std::pair<int, float> a, std::pair<int, float> b) { return a.second < b.second; });
  }

	return time_span.count();
}

int main(int argc, char **argv) {
	
	if (argc != 5 ) {
		std::cerr << "Error: Wrong number of arguments\n";
		std::cerr << "First argument is the number of points\n";
		std::cerr << "Second argument is the dimension of the pointset\n";
		std::cerr << "Third argument is the number of nearest neighbour\n";
		std::cerr << "Fourth argument is a parameter for the LSH function\n";
		exit(-1);
	}
  
  size_t n = 10000;
  //size_t k = floor(log2(n)/2);
  size_t d = 128;
  size_t q = 100;
  size_t m = 10;
  
  std::vector<Point> pointset(n, Point(d));
  
  std::vector<Point> queries(q, Point(d));
  
  
  
  
  readfvecs(pointset, n, d, "./siftsmall/siftsmall_base.fvecs");
  readfvecs(queries, q, d, "./siftsmall/siftsmall_query.fvecs");
  
  Point mean;
  for(size_t i=0;i<d;++i){
  	mean.push_back(0);
  }
  for(size_t j=0;j<n;++j){
  	for(size_t i=0;i<d;++i){
  		mean[i] += pointset[j][i] ;
  	}
  }
  for(size_t i=0;i<d;++i){
  	mean[i]=mean[i]/n;
  }
  for(size_t j=0;j<n;++j){
  	for(size_t i=0;i<d;++i){
  		pointset[j][i] = pointset[j][i] - mean[i];
  	}
  }
  for(size_t j=0;j<q;++j){
  	for(size_t i=0;i<d;++i){
  		queries[j][i] = queries[j][i] - mean[i];
  	}
  }
  
	
	std::cout << "Data generated\n";
	
	
  /*high_resolution_clock::time_point t1 = high_resolution_clock::now();
  
  dolphinn::Dolphinn<T, bitT> dolphi(pointset, n, d, k, atof(argv[4])); 
  std::cout << "Hypercube built\n";
  
	high_resolution_clock::time_point t2 = high_resolution_clock::now();

	duration<double> time_span = duration_cast<duration<double>>(t2 - t1);

	std::cout << "Build: " << time_span.count() << " seconds.\n";
  
  t1 = high_resolution_clock::now();
  dolphi.m_nearest_neighbors_query(queries, q, m, 1000, result, 1);
  t2 = high_resolution_clock::now();
  time_span = duration_cast<duration<double>>(t2 - t1);
  std::cout << "Query: " << time_span.count() << " seconds.\n";
  std::cout << "radius queries done\n";*/
  std::vector<std::vector<std::pair<int, float>>> resultn(q);
  
  //t1 = high_resolution_clock::now();
  for(size_t i=0;i<q;++i){
  	if(i%100==0){std::cout <<i<<"\n";}
  	for(size_t j=0;j<n;j++){
  		if(j<m){
				resultn[i].push_back(std::pair<int, float> (j,squared_Eucl_distance(&pointset[j],&queries[i])));
				std::sort(resultn[i].begin(),resultn[i].end(),[](std::pair<int, float> a, std::pair<int, float> b) { return a.second < b.second; });
			} else {
				float tmpd=squared_Eucl_distance(&pointset[j],&queries[i]);
				if(tmpd<resultn[i][m-1].second){
					resultn[i][m-1].first=j;
					resultn[i][m-1].second=tmpd;
					std::sort(resultn[i].begin(),resultn[i].end(),[](std::pair<int, float> a, std::pair<int, float> b) { return a.second < b.second; });
				}
			}
  	}
  }
  //t2 = high_resolution_clock::now();
  //time_span = duration_cast<duration<double>>(t2 - t1);
  //std::cout << "Naive: " << time_span.count() << " seconds.\n";
  
	std::cout << "naive done\n";
  std::vector<std::vector<std::pair<int, float>>> result;
	for(size_t i=0;i<q;++i){
		std::vector<std::pair<int, float>> dummy(m);
		result.push_back(dummy);
	}
  
  std::map<float, std::pair<double , std::pair<size_t,size_t>>> times;
  
  for(size_t K = 1; K<16; ++K) {
  	for(size_t ex = 10; ex<n; ex=floor(1.3*ex)) {
  		double t = query_test(result, n, d, m, q, pointset, queries, K, ex);
  		float comp = compare (q, m, result, resultn);
			std::cout << "K=" << K << "  ex=" << ex << "  t=";
  		std::cout << t;
  		std::cout << "  acc=";
  		std::cout << comp;
  		std::cout << "\n";
  		if(times.find(comp) != times.end()){
  			if(times[comp].first>t) {times[comp] = std::pair<double , std::pair<size_t,size_t>> (t,std::pair<size_t,size_t> (K,ex));}
  		} else {
  			auto it = times.upper_bound(comp);
  			if(it == times.end()) {
  				times[comp] = std::pair<double , std::pair<size_t,size_t>> (t,std::pair<size_t,size_t> (K,ex));
  			} else {
  				if(it->second.first>t){
  					times[comp] = std::pair<double , std::pair<size_t,size_t>> (t,std::pair<size_t,size_t> (K,ex));
  				}
  				auto it2 = times.find(comp);
  				if(it2 != times.begin()){
  					--it2;
  					if(it2->second.first>t){
  						times.erase(it2);
  					}
  				}
  			}
  		}
  	}
  }
  std::cout << "\n";std::cout << "\n";std::cout << "\n";
  for(auto x: times){
  	std::cout << "acc= " << x.first << "  t= " << x.second.first << "  K= " << x.second.second.first << "  ex= " << x.second.second.second << "\n";
  }
  
  
  
  
  /*float correct =0;
  std::vector<int>vect(11);
  for(size_t i=0;i<q;++i){
  	int c=0;
  	float tmp = correct;
  	for(size_t j=0;j<m;j++){
  		if(result[i][c].first == resultn[i][j].first){++correct;++c;}
  	}
  	std::cout << correct-tmp<<"\n";
  	vect[correct-tmp] += 1;
  	if(correct-tmp==0){
			for (auto& x:resultn[i]){
				std::cout <<"("<< x.first<<","<< x.second<<")  ";
			}
			std::cout << "\n";
			for (auto& x:result[i]){
				std::cout <<"("<< x.first<<","<< x.second<<")  ";
			}std::cout << "\n";
  	}
  }
	for(auto x:vect)
		std::cout << x <<" ";std::cout << "\n";
  std::cout << "correct:" << correct*100/(q*m);std::cout << "\n";*/
  return 0;
}
