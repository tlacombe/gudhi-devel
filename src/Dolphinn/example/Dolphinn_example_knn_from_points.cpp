#include <iostream>
#include <gudhi/Dolphinn.h>
#include <gudhi/Points_off_io.h>
#include <vector>
#include <random>
#include <cmath> // for sqrt
#include <utility> // for pair

#include <ctime>
#include <ratio>
#include <chrono>

#define T double
#define Point std::vector<T>


using Dolphinn = Gudhi::dolphinn::Dolphinn<T>;

int main(int argc, char **argv) {
	
	if (argc != 5 ) {
		std::cerr << "Error: Wrong number of arguments\n";
		std::cerr << "First argument is the number of points\n";
		std::cerr << "Second argument is the dimension of the pointset\n";
		std::cerr << "Third argument is the number of nearest neighbors\n";
		std::cerr << "Fourth argument is a parameter for the LSH function\n";
		exit(-1);
	}
  
  size_t n = atoi(argv[1]);
  size_t k = floor(log2(n)/2);
  size_t d = atoi(argv[2]);
  
  
  
  std::vector<Point> pointset;
  
  std::vector<Point> queries;
  
  std::default_random_engine generator(std::chrono::system_clock::now().time_since_epoch().count());
  
  std::normal_distribution<double> distribution(0.0,1.0/std::sqrt((double)d));
  
  for(size_t i=0;i<n;++i){
  	Point p;
  	for(size_t j=0;j<d;++j){
  		p.push_back(distribution(generator));
  	}
  	pointset.push_back(p);
  }
  
  Point query;
  for(size_t j=0;j<d;++j){
		query.push_back(distribution(generator));
	}
	queries.push_back(query);
	
	
	std::vector<std::vector<std::pair<int, double>>> result;
	std::vector<std::pair<int, double>> dummy;
	result.push_back(dummy);
	
	
	std::cout << "Data generated\n";
	
	using namespace std::chrono;
  high_resolution_clock::time_point t1 = high_resolution_clock::now();
  
  Dolphinn dolphi(pointset, n, d, k, atof(argv[4])); 
  std::cout << "Hypercube built\n";
  
	high_resolution_clock::time_point t2 = high_resolution_clock::now();

	duration<double> time_span = duration_cast<duration<double>>(t2 - t1);

	std::cout << "Build: " << time_span.count() << " seconds.\n";
  
  dolphi.k_nearest_neighbors_query(queries, 1, atoi(argv[3]), n/100+atoi(argv[3]), result);
  std::cout << "radius queries done\n";
  
  std::cout << "The query is:";
  for(auto& x:query) std::cout << x << " ";
  std::cout << "\n";
  
  
  std::cout << "The nearest neighbors are: \n";
  for(auto& x:result[0]){
  	for(auto y:pointset[x.first])
  		std::cout << y << " ";
  	//std::cout << x.first << " ";
  	std::cout << " for a distance of " << std::sqrt(x.second) << "\n";
  }
  
  return 0;
}
