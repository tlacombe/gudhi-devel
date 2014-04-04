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
#include "Nearest_neighbors.h"
#include "Simplex_tree.h"

typedef std::vector< double > Point;
typedef std::vector< Point >  Point_range;
typedef int                   Vertex_handle;
typedef double                Filtration_value;

bool test_flag_complex_construction ( std::string filegraph
                                    , double      threshold 
                                    , int         max_dim )
{
  clock_t start, end;
  //Construct the Simplex Tree
  Simplex_tree<> st;
  start = clock();
  auto g = read_graph(filegraph); 
  st.insert_graph (g); //insert the graph in the simplex tree as 1-skeleton
  end = clock();

  std::cout << "Insert the 1-skeleton in the simplex tree in "
       << (double)(end-start)/CLOCKS_PER_SEC << " s. \n";

  start = clock();
  st.expansion ( max_dim ); //expand the 1-skeleton until dimension max_dim
  end = clock();

  std::cout << "Expand the simplex tree in "
       << (double)(end-start)/CLOCKS_PER_SEC << " s. \n";

  std::cout << "Information of the Simplex Tree: " << std::endl;
  std::cout << "  Number of vertices = " << st.num_vertices() << " ";
  std::cout << "  Number of simplices = " << st.num_simplices() << std::endl;
  std::cout << std::endl << std::endl;
  return true;
}

bool test_simplex_tree_iterators ( std::string filegraph
                                 , double      threshold
                                 , int         max_dim )
{
  clock_t start, end;
  //Construct the Simplex Tree
  Simplex_tree<> st;
  
  start = clock();
  auto g = read_graph(filegraph); 
  st.insert_graph (g); //insert the graph in the simplex tree as 1-skeleton
  end = clock();
  std::cout << "Insert the 1-skeleton in the simplex tree in "
       << (double)(end-start)/CLOCKS_PER_SEC << " s. \n";

  start = clock();
  st.expansion ( max_dim ); //expand the 1-skeleton until dimension max_dim
  end = clock();
  std::cout << "Expand the simplex tree in "
       << (double)(end-start)/CLOCKS_PER_SEC << " s. \n";

  std::cout << "Iterator on vertices: ";
    for( auto vertex : st.complex_vertex_range() ) { std::cout << vertex << " "; }
  std::cout << std::endl;

  std::cout << std::endl << std::endl;

  std::cout << "Iterator on simplices: " << std::endl;
  for( auto simplex : st.complex_simplex_range() ) 
  { 
    std::cout << "   ";
    for( auto vertex : st.simplex_vertex_range(simplex) ) { std::cout << vertex << " "; }
    std::cout << std::endl;
  }

  std::cout << std::endl << std::endl;

  std::cout << "Iterator on Simplices in the filtration, with [filtration value]:" << std::endl;
  for( auto f_simplex : st.filtration_simplex_range() )
  { std::cout << "   " << "[" << st.filtration(f_simplex) << "] "; 
    for( auto vertex : st.simplex_vertex_range(f_simplex) ) 
      { std::cout << vertex << " "; } std::cout << std::endl;
  }

  std::cout << std::endl << std::endl;

  std::cout << "Iterator on Simplices in the filtration, and their boundary simplices:" << std::endl;
  for( auto f_simplex : st.filtration_simplex_range() )
  {
    std::cout << "   " << "[" << st.filtration(f_simplex) << "] "; 
    for( auto vertex : st.simplex_vertex_range(f_simplex) ) 
      { std::cout << vertex << " "; } std::cout << std::endl;

    for( auto b_simplex : st.boundary_simplex_range(f_simplex) )
    {
      std::cout << "      " << "[" << st.filtration(b_simplex) << "] "; 
    for( auto vertex : st.simplex_vertex_range(b_simplex) ) 
      { std::cout << vertex << " "; } std::cout << std::endl;
    }
  }
  return true;
}

bool test_simplex_tree_to_file( std::string filegraph
                              , double threshold
                              , int max_dim
                              , std::string filecomplex)
{
  auto g = read_graph(filegraph); 
  Simplex_tree<> st;
  st.insert_graph (g); //insert the graph in the simplex tree as 1-skeleton
  st.expansion ( max_dim ); //expand the 1-skeleton until dimension max_dim
  
  std::ofstream out_ (filecomplex.c_str(), std::ios::trunc);
  for(auto sh : st.complex_simplex_range())
  {
    out_ << st.dimension(sh) << " ";
    for(auto v : st.simplex_vertex_range(sh))
      { out_ << v << " ";}
    out_ << st.filtration(sh) << " \n";
  }
  return true;
}

bool test_file_to_simplex_tree(std::string filecomplex)
{
  std::ifstream in_ (filecomplex.c_str(),std::ios::in);
  if(!in_.is_open()) { std::cerr << "Unable to open file " << filecomplex << std::endl; }

  std::vector<Vertex_handle> simplex;
  Filtration_value           fil;
  Simplex_tree<>             st;

  while(read_simplex(in_,simplex,fil)) //read all simplices in the file as a list of vertices
  {
    st.insert(simplex,fil); //insert every simplex in the simplex tree
    simplex.clear();
  }
  in_.close();  
  return true;
}

int main (int argc, char * const argv[]) 
{
  if (argc != 6) {
    std::cerr << "Usage: " << argv[0] 
      << " path_to_file_points path_to_file_graph path_to_file_complex threshold dim_expansion \n";  
    return 0;
  }
  std::string filepoints  = argv[1]; 
  std::string filegraph   = argv[2];
  std::string filecomplex = argv[3];
  double threshold        = atof(argv[4]);
  int dim_expansion       = atoi(argv[5]);

  compute_rips_graph(filepoints,filegraph,threshold);

  test_flag_complex_construction(filegraph,threshold,dim_expansion);

  test_simplex_tree_iterators(filegraph,threshold,dim_expansion);

  test_simplex_tree_to_file(filegraph,threshold,dim_expansion,filecomplex);

  test_file_to_simplex_tree(filecomplex);

  return 0;
}
