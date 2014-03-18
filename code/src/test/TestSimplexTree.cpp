/*
 *  TestSimplexTree.h
 *  Gudhi
 *
 *  Created by Clément Maria on 02/23/14.
 *  Copyright 2013 INRIA. All rights reserved.
 *
 */

#include <iostream>
#include <ctime>
#include "utils/iofile.h"
#include "geometry/Euclidean_geometry.h"
#include "geometry/Rips_graph_naive.h"
#include "combinatorics/Simplex_tree/Simplex_tree.h"

using namespace std;

typedef vector< double > Point;
typedef vector< Point >  Point_range;

bool test_points_extraction(string filename)
{
  Point_range points;              //read the points from the file
  read_points( filename, points ); //and turn them into a Point_range
  return true;
}

bool test_flag_complex_construction ( string filename
                                    , double threshold 
                                    , int max_dim )
{
  Point_range points;              //read the points from the file
  read_points( filename, points ); //and turn them into a Point_range

  // Create a metric space from the points, with Euclidean metric:
  Euclidean_geometry< Point > ms;
  ms.init(points);

  // Create a NeighborGraph with the space:
  Rips_graph_naive< Euclidean_geometry< Point > > ng( ms, threshold );


  clock_t start, end;
  //Construct the Simplex Tree
  Simplex_tree< Euclidean_geometry< Point > > st ( ms );
  start = clock();
  st.insert_graph ( ng ); //insert the graph in the simplex tree as 1-skeleton
  end = clock();

  cout << "Insert the 1-skeleton in the simplex tree in "
       << (double)(end-start)/CLOCKS_PER_SEC << " s. \n";

  start = clock();
  st.expansion ( max_dim ); //expand the 1-skeleton until dimension max_dim
  end = clock();

  cout << "Expand the simplex tree in "
       << (double)(end-start)/CLOCKS_PER_SEC << " s. \n";


  cout << "Information of the Simplex Tree: " << endl;
  cout << "  Number of vertices = " << st.nb_vertices() << " ";
  cout << "    for " << points.size() << " points " << endl;
  cout << "  Number of simplices = " << st.nb_simplices() << endl;
  cout << endl << endl;
  return true;
}

bool test_simplex_tree_iterators ( string filename
                                 , double threshold
                                 , int max_dim )
{
  Point_range points;              //read the points from the file
  read_points( filename, points ); //and turn them into a Point_range

  // Create a metric space from the points, with Euclidean metric:
  Euclidean_geometry< Point > ms;
  ms.init(points);

  // Create a NeighborGraph with the space:
  Rips_graph_naive< Euclidean_geometry< Point > > ng( ms, threshold );

  //Construct the Simplex Tree
  Simplex_tree< Euclidean_geometry< Point > > st ( ms );
  st.insert_graph ( ng ); //insert the graph in the simplex tree as 1-skeleton

  st.expansion ( max_dim ); //expand the 1-skeleton until dimension max_dim

  cout << "Iterator on vertices: ";
    for( auto vertex : st.complex_vertex_range() ) { cout << vertex << " "; }
  cout << endl;

  cout << endl << endl;

  cout << "Iterator on simplices: " << endl;
  for( auto simplex : st.complex_simplex_range() ) 
  { 
    cout << "   ";
    for( auto vertex : st.simplex_vertex_range(simplex) ) { cout << vertex << " "; }
    cout << endl;
  }

  cout << endl << endl;

  cout << "Iterator on Simplices in the filtration, with [filtration value]:" << endl;
  for( auto f_simplex : st.filtration_simplex_range() )
  { cout << "   " << "[" << st.filtration(f_simplex) << "] "; 
    for( auto vertex : st.simplex_vertex_range(f_simplex) ) 
      { cout << vertex << " "; } cout << endl;
  }

  cout << endl << endl;

  cout << "Iterator on Simplices in the filtration, and their boundary simplices:" << endl;
  for( auto f_simplex : st.filtration_simplex_range() )
  {
    cout << "   " << "[" << st.filtration(f_simplex) << "] "; 
    for( auto vertex : st.simplex_vertex_range(f_simplex) ) 
      { cout << vertex << " "; } cout << endl;

    for( auto b_simplex : st.boundary_simplex_range(f_simplex) )
    {
      cout << "      " << "[" << st.filtration(b_simplex) << "] "; 
    for( auto vertex : st.simplex_vertex_range(b_simplex) ) 
      { cout << vertex << " "; } cout << endl;
    }
  }
  return true;
}

int main (int argc, char * const argv[]) 
{
  if (argc != 4) {
    std::cerr << "Usage: " << argv[0] 
    << " path_to_data_points threshold dim_expansion \n";  
    return -1;
  }
 
  test_points_extraction(argv[1]);
  test_flag_complex_construction(argv[1],atof(argv[2]),atoi(argv[3]));
  test_simplex_tree_iterators(argv[1],atof(argv[2]),atoi(argv[3]));

  return 0;
}
