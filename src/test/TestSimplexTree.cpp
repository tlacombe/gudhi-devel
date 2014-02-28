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

bool test_flag_complex_construction(string filename)
{
  Point_range points;              //read the points from the file
  read_points( filename, points ); //and turn them into a Point_range

  // Create a metric space from the points, with Euclidean metric:
  Euclidean_geometry< Point > ms;
  ms.init(points);

  // Create a NeighborGraph with the space:
  double threshold = 100; //max distance between neighbors
  Rips_graph_naive< Euclidean_geometry< Point > > ng( ms, threshold );

  //Construct the Simplex Tree
  Simplex_tree< Euclidean_geometry< Point > > st ( ms );
  st.insert_graph ( ng ); //insert the graph in the simplex tree as 1-skeleton

  int max_dim = 10;
  st.expansion ( max_dim ); //expand the 1-skeleton until dimension max_dim

  cout << "Information of the Simplex Tree: " << endl;
  cout << "  Number of vertices = " << st.nb_vertices() << " ";
  cout << "    for " << points.size() << " points " << endl;
  cout << "  Number of simplices = " << st.nb_simplices() << endl;
  cout << endl << endl;
  return true;
}

bool test_simplex_tree_iterators(string filename)
{
  Point_range points;              //read the points from the file
  read_points( filename, points ); //and turn them into a Point_range

  // Create a metric space from the points, with Euclidean metric:
  Euclidean_geometry< Point > ms;
  ms.init(points);

  // Create a NeighborGraph with the space:
  double threshold = 100; //max distance between neighbors
  Rips_graph_naive< Euclidean_geometry< Point > > ng( ms, threshold );

  //Construct the Simplex Tree
  Simplex_tree< Euclidean_geometry< Point > > st ( ms );
  st.insert_graph ( ng ); //insert the graph in the simplex tree as 1-skeleton

  int max_dim = 10;
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
  string test_points = "../test/Simplex_tree_points.dat";
  test_points_extraction(test_points);
  test_flag_complex_construction(test_points);
  test_simplex_tree_iterators(test_points);

  return 0;
}
