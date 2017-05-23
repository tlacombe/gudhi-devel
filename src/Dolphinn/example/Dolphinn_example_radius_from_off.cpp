#include <iostream>
#include <gudhi/Dolphinn.h>
#include <gudhi/Points_off_io.h>
#include <vector>

#define T float
#define bitT char
#define Point std::vector<T>


using Dolphinn = Gudhi::dolphinn::Dolphinn<T, bitT>;

int main(int argc, char **argv) {
	
	if (argc != 5 && argc != 6) {
		std::cerr << "Error: Wrong number of arguments\n";
		std::cerr << "First argument is an .off data file\n";
		std::cerr << "Second argument is an .off query file\n";
		std::cerr << "Third argument is the radius\n";
		std::cerr << "Fourth argument is the LSH's windows size\n";
		std::cerr << "(optionnal) Fifth argument is the number of queries (the n firsts queries of the query file)\n";
		exit(-1);
	}
  
  Gudhi::Points_off_reader<Point> off_readerd(argv[1]);
  if (!off_readerd.is_valid()) {
    std::cerr << "Unable to read file " << "\n";
    exit(-1);  
  }
  std::vector<Point> pointset = off_readerd.get_point_cloud();
  size_t n = pointset.size();
  size_t k = floor(log2(n)/2)+1;
  size_t d = pointset[0].size();
  
  Gudhi::Points_off_reader<Point> off_readerq(argv[2]);
  if (!off_readerq.is_valid()) {
    std::cerr << "Unable to read file " << "\n";
    exit(-1);  
  }
	std::vector<Point> queries = off_readerq.get_point_cloud();  
  
  int nq;
  if(argc == 6){nq=atoi(argv[5]);} else{nq=queries.size();};
  
  std::vector<int> result(nq);
  
  Dolphinn dolphi(pointset, n, d, k, atof(argv[4])); 
  std::cout << "Hypercube built \n";
  
  dolphi.radius_query(queries, nq, atof(argv[3]), n/100 + 10, result, 1);
  std::cout << "radius queries done\n";
  
  for(int i=0; i<nq; ++i){
  	if(result[i]!=-1){
  		std::cout << "Neighbour of ";
  		for(auto x: queries[i]) std::cout << x << " ";
  		std::cout << "found: ";
  		for(auto x: pointset[result[i]]) std::cout << x << " ";
  		std::cout << "d=";
  		float res=0;
  		for(size_t j=0;j<d;++j){
  			float tmp = pointset[result[i]][j] - queries[i][j];
				res += tmp*tmp;
			}
  		std::cout << res <<"\n";
  	} else {
  		std::cout << "No neighbour found for ";
  		for(auto x: queries[i]) std::cout << x << " ";
  		std::cout << "\n";
  	}
  }
  
  return 0;
}
