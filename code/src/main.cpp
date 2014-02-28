#include <iostream>
#include <ctime>

#include "utils/iofile.h"

#include "boost/pending/disjoint_sets.hpp"
#include "boost/iterator/counting_iterator.hpp"
#include "boost/range/counting_range.hpp"

#include "boost/container/flat_map.hpp"

#include "geometry/Euclidean_geometry.h"
#include "geometry/Rips_graph_naive.h"
#include "combinatorics/Simplex_tree/Simplex_tree.h"

#include "topology/Persistent_cohomology.h"

using namespace std;

typedef std::vector<double> Point;           // in order not to copy many times the Points.
typedef std::vector<Point>  Point_range;



Simplex_tree< Euclidean_geometry< Point > > * test_simplex_tree ()
{
  // Extract data points from file file_name.
  // Turn them into a Point_range object:    <--- has to be templated; 
                                                //best, use a istream_iterator
  Point_range points;
  string file_name = "/Users/cmaria/Desktop/Points.txt";
  read_points(file_name,points);
  // Create a metric space from the points, with euclidean metric:
  Euclidean_geometry< Point > ms;
  ms.init(points);                 
  // Create a NeighborGraph with the space:
  double threshold = 100;
  Rips_graph_naive< Euclidean_geometry< Point > > ng(ms,threshold);  
  // Create a simplex_tree
  Simplex_tree< Euclidean_geometry< Point > > * st_ptr 
  = new Simplex_tree< Euclidean_geometry< Point > >(ms); //constructor
  st_ptr->insert_graph(ng); //insert the graph

  int max_dim = 10;
  std::cout << "Expand the flag complex \n";
  st_ptr->expansion(max_dim);
//  std::cout << st << std::endl;
  return st_ptr;
}


int main (int argc, char * const argv[]) 
{
  auto st_ptr = test_simplex_tree (); //++st_ptr;
  Persistent_cohomology < Simplex_tree< Euclidean_geometry< Point > >  > pcoh(st_ptr);
  pcoh.compute_persistent_cohomology();

  return 0;

  /* TO DO
   PcoH<Flag_simplex_tree> co;
   co.init(rho_max);  // CAM->does nothing
                      // PH initialize the matrix
   intervals_begin(); // <- does the work
  //end
  */

}
