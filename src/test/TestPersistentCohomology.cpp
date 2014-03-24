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

#include "topology/compare_diagrams.h"

using namespace std;

typedef vector< double >                          Point;
typedef vector< Point >                           Point_range;
typedef Euclidean_geometry< Point >               Metric_space;
typedef Rips_graph_naive< Metric_space >          Neighbor_graph;
typedef Simplex_tree< Metric_space >              Simplicial_cpx;
typedef Persistent_cohomology< Simplicial_cpx >   PcoH;

int main (int argc, char * const argv[]) {




std::cout << "\n \n \n \n \n \n";
compare_diagrams ( "/Users/cmaria/Desktop/gudhi_diag.txt"
               //  , "/Users/cmaria/Desktop/dio_diag.txt" );
                  , "/Users/cmaria/Desktop/diag_steve.txt" );
return 0;


  clock_t start, end;

// Construct a filtered Rips complex:
  string filename = "../data/Bro_fp.dat";
  double threshold = 0.75; //max distance between neighbors
  int max_dim = 5;

  Point_range points;              //read the points from the file
  read_points( filename, points ); //and turn them into a Point_range
  // Create a metric space from the points, with Euclidean metric:
  Euclidean_geometry< Point > ms;
  ms.init(points);
  // Create a NeighborGraph with the space:
  Rips_graph_naive< Euclidean_geometry< Point > > ng( ms, threshold );
  //Construct the Simplex Tree
  start = clock();
  Simplex_tree< Euclidean_geometry< Point > > st ( ms );
  end = clock();
  cout << "Initialize the simplex tree in: "
       << (double)(end-start)/CLOCKS_PER_SEC << " s. \n";

  start = clock();
  st.insert_graph ( ng );   //insert the graph in the simplex tree as 1-skeleton
  end = clock();
  cout << "Insert the 1-skeleton in the simplex tree in: "
       << (double)(end-start)/CLOCKS_PER_SEC << " s. \n";

  start = clock();
  st.expansion ( max_dim ); //expand the 1-skeleton until dimension max_dim
  end = clock();
  cout << "Expand the simplex tree in: "
       << (double)(end-start)/CLOCKS_PER_SEC << " s. \n";

  cout << "Num simplices = " << st.num_simplices() << endl;


  ofstream outfil;
  outfil.open ("/Users/cmaria/Desktop/st_fil.txt");
  cout << "Iterator on Simplices in the filtration, with [filtration value]:" << endl;
  for( auto f_simplex : st.filtration_simplex_range() )
  { 
    outfil << st.dimension(f_simplex) << " ";
    for( auto vertex : st.simplex_vertex_range(f_simplex) ) { outfil << vertex << " "; } 
    outfil << st.filtration(f_simplex) << " \n"; 
  }
  outfil.close();


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


  ofstream outdiag;
  outdiag.open ("/Users/cmaria/Desktop/gudhi_diag.txt");
  pcoh.output_diagram(outdiag);
  outdiag.close();




// std::cout << "\n \n \n \n \n \n";
// compare_diagrams ( "/Users/cmaria/Desktop/gudhi_diag.txt"
//                //  , "/Users/cmaria/Desktop/dio_diag.txt" );
//                   , "/Users/cmaria/Desktop/diag_steve.txt" );
// return 0;





  return 0;
}
