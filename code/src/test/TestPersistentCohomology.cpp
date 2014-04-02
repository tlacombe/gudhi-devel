/*
 *  TestPersistentCohomology.h
 *  Gudhi
 *
 *  Created by Clément Maria on 02/28/14.
 *  Copyright 2013 INRIA. All rights reserved.
 *
 */

#include <ctime>
#include "utils/iofile.h"
#include "geometry/Euclidean_geometry.h"
#include "geometry/Rips_graph_naive.h"
#include "combinatorics/Simplex_tree/Simplex_tree.h"
#include "topology/Persistent_cohomology.h"

//#include "topology/compare_diagrams.h"

using namespace std;



typedef vector< double >                          Point;
typedef vector< Point >                           Point_range;
typedef Euclidean_geometry< Point >               Metric_space;
typedef Rips_graph_naive< Metric_space >          Neighbor_graph;
typedef Simplex_tree<>                            Simplicial_cpx;
typedef Persistent_cohomology< Simplicial_cpx >   PcoH;




void compute_rips_graph( std::string filepoints
                       , std::string filegraph 
                       , double threshold )
{
  Point_range points;              //read the points from the file
  read_points( filepoints, points ); //and turn them into a Point_range
  // Create a metric space from the points, with Euclidean metric:
  Euclidean_geometry< Point > ms;
  ms.init(points);
  // Create a NeighborGraph with the space:
  Neighbor_graph ng( ms, threshold );
  
  ofstream outgraph;
  outgraph.open ( filegraph );
  if( !outgraph.is_open() ) {
    std::cerr << "Unable to open file " << filegraph << std::endl;
    return;}

  for(auto u : ms.space_vertex_range())
  {
    outgraph << 0 << " " << u << " " << 0. << std::endl;
    for(auto v : ng.adjacent_vertices(u)) 
    { outgraph << 1 << " " << u << " " << v << " " << ms.distance(u,v) << std::endl; }
  }
  outgraph.close();
}


int main (int argc, char * const argv[]) {


  clock_t start, end;

// Construct a filtered Rips complex:
  // string filename = "../data/Cy8.txt";
  string filepoints = "/Users/cmaria/Downloads/rings.txt";
  string filegraph  = "/Users/cmaria/Desktop/skel_graph.txt";
  double threshold = 1.5; //max distance between neighbors
  int max_dim = 3;

  compute_rips_graph( filepoints, filegraph, threshold );

  //Construct the Simplex Tree
  start = clock();
  Simplex_tree<> st;
  end = clock();
  cout << "Initialize the simplex tree in: "
       << (double)(end-start)/CLOCKS_PER_SEC << " s. \n";


  auto g =  read_graph(filegraph);

  start = clock();
  st.insert_graph (g);   //insert the graph in the simplex tree as 1-skeleton
  end = clock();
  cout << "Insert the 1-skeleton in the simplex tree in: "
       << (double)(end-start)/CLOCKS_PER_SEC << " s. \n";

  start = clock();
  st.expansion ( max_dim ); //expand the 1-skeleton until dimension max_dim
  end = clock();
  cout << "Expand the simplex tree in: "
       << (double)(end-start)/CLOCKS_PER_SEC << " s. \n";

  cout << "Num simplices = " << st.num_simplices() << endl;

  cout << "                  Construction time per simplex = " << 
  ((double)(end-start)/CLOCKS_PER_SEC) / (double)st.num_simplices() << endl;


  start = clock();
  PcoH pcoh (st);
  end = clock();
  cout << "Initialize the PcoH structures in: "
       << (double)(end-start)/CLOCKS_PER_SEC << " s. \n";

  start = clock();
  pcoh.compute_persistent_cohomology ();
  end = clock();

  cout << "Compute persistent cohomology in: " 
       << (double)(end-start)/CLOCKS_PER_SEC << " s. \n";

  cout << "                  PcoH time per simplex = " << 
  ((double)(end-start)/CLOCKS_PER_SEC) / (double)st.num_simplices() << endl;

  ofstream outdiag;
  outdiag.open ("/Users/cmaria/Desktop/gudhi_diag.txt");
  pcoh.output_diagram(outdiag);
  outdiag.close();


  return 0;
}
