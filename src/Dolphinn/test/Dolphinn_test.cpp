#include <iostream>
#include <assert.h> 
#include "gudhi/Dolphinn.h"
//#include <CGAL/Epick_d.h>
#include <gudhi/Points_off_io.h>

#include <ctime>
#include <ratio>
#include <chrono>
#include <vector>

#define N 10000
#define D 128
#define K floor(log2(N)/2)
#define M 1// #NN to find 
#define Q 100
#define T float
#define bitT char
#define MAX_PNTS_TO_SEARCH N/100
#define THREADS_NO 1
#define RADIUS 1
//#define Point CGAL::Cartesian_d<T>::Point_d
#define Point std::vector<float>


using namespace Gudhi;

int main()
{
	
  Gudhi::Points_off_reader<Point> off_reader("./../../../../data/points/bunny_5000.off");
  if (!off_reader.is_valid()) {
    std::cerr << "Unable to read file " << "\n";
    exit(-1);  // ----- >>
  }
  std::vector<Point> pointset = off_reader.get_point_cloud();
  size_t n = pointset.size();
  size_t k = floor(log2(N)/2);
  size_t d = pointset[0].size();
  
  Point query = pointset[n-1];
  pointset.pop_back();
  
  std::vector<Point> queries;
  queries.push_back(query);
  
  std::vector<int> result;
  
  result.push_back(21); 
  
  dolphinn::Dolphinn<T, bitT> dolphi(pointset, n, d, k, 0.001); 
  dolphi.radius_query(queries, 1, 10, 500, result, 1);
  
  float res=0;
	float tmp;
	std::cout << "\n";
	for(auto x:result) std::cout << x << " ";
	std::cout << "\n";
	for(auto x:pointset[result[0]]) std::cout << x << " ";
	std::cout << "\n";
	for(auto x:query) std::cout << x << " ";
	std::cout << "\n";
	
	for(size_t i=0;i<d;++i){
		tmp = pointset[result[0]][i] - query[i];
		res += tmp*tmp;
	}
  std::cout << res << " res\n";
  return 0;
}
