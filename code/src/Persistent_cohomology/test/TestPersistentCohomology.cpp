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
#include "Nearest_neighbors.h"
#include "Simplex_tree.h"
#include "Hasse_complex.h"
#include "Persistent_cohomology.h"
#include "Persistent_cohomology/Field_Zp.h"
#include "Persistent_cohomology/Multi_field.h"


//#include "topology/compare_diagrams.h"


typedef std::vector< double >                                        Point;
typedef std::vector< Point >                                         Point_range;
typedef Euclidean_geometry< Point >                             Metric_space;
typedef Rips_graph_naive< Metric_space >                        Neighbor_graph;
typedef Simplex_tree<>                                          Simplicial_cpx;
typedef Hasse_complex<>                                         Hasse_cpx;
typedef Persistent_cohomology< Simplicial_cpx, Field_Zp<11> >   PcoH;

typedef Persistent_cohomology< Hasse_cpx, Multi_field<2,3> >    PcoHasse;


int main (int argc, char * const argv[]) 
{
  clock_t start, end;

// Construct a filtered Rips complex:
//   string filepoints = "../data/Cy8.txt";
  std::string filepoints = "/Users/cmaria/Downloads/rings.txt";
  std::string filegraph  = "/Users/cmaria/Desktop/skel_graph.txt";
  
  double threshold = 1.5; //max distance between neighbors
  int max_dim = 3;

  compute_rips_graph( filepoints, filegraph, threshold );

  //Construct the Simplex Tree
  start = clock();
  Simplex_tree<> st;
  end = clock();
  std::cout << "Initialize the simplex tree in: "
            << (double)(end-start)/CLOCKS_PER_SEC << " s. \n";

  auto g = read_graph(filegraph);

  start = clock();
  st.insert_graph (g);   //insert the graph in the simplex tree as 1-skeleton
  end = clock();
  std::cout << "Insert the 1-skeleton in the simplex tree in: "
            << (double)(end-start)/CLOCKS_PER_SEC << " s. \n";

  start = clock();
  st.expansion ( max_dim ); //expand the 1-skeleton until dimension max_dim
  end = clock();
  std::cout << "Expand the simplex tree in: "
            << (double)(end-start)/CLOCKS_PER_SEC << " s. \n";

  std::cout << "Num simplices = " << st.num_simplices() << std::endl;

  std::cout << "                  Construction time per simplex = " << 
  ((double)(end-start)/CLOCKS_PER_SEC) / (double)st.num_simplices() << std::endl;

  start = clock();
  st.initialize_filtration();
  end = clock();
  std::cout << "Order the simplices according to the filtration in: "
            << (double)(end-start)/CLOCKS_PER_SEC << " s. \n";

  start = clock();
  Hasse_cpx hcpx(&st);
  end = clock();
  std::cout << "Convert the Simplex Tree into a Hasse Complex in: "
            << (double)(end-start)/CLOCKS_PER_SEC << " s. \n";

  start = clock();
  PcoH pcoh (st);
  end = clock();
  std::cout << "With ST: Initialize the PcoH structures in: "
            << (double)(end-start)/CLOCKS_PER_SEC << " s. \n";

  start = clock();
  pcoh.compute_persistent_cohomology (0.2);
  end = clock();

  std::cout << "With ST: Compute persistent cohomology in: " 
            << (double)(end-start)/CLOCKS_PER_SEC << " s. \n";

  std::cout << "                  PcoH time per simplex = " << 
    ((double)(end-start)/CLOCKS_PER_SEC) / (double)st.num_simplices() << std::endl;


 
  start = clock();
  PcoHasse pcoh_ (hcpx);
  end = clock();
  std::cout << "With HC: Initialize the PcoH structures in: "
            << (double)(end-start)/CLOCKS_PER_SEC << " s. \n";

  start = clock();
  pcoh_.compute_persistent_cohomology (0.2);
  end = clock();

  std::cout << "With HC: Compute persistent cohomology in: " 
            << (double)(end-start)/CLOCKS_PER_SEC << " s. \n";

  std::cout << "                  PcoH time per simplex = " << 
    ((double)(end-start)/CLOCKS_PER_SEC) / (double)st.num_simplices() << std::endl;

  std::ofstream outdiag;
  outdiag.open ("/Users/cmaria/Desktop/gudhi_diag.txt");
  pcoh.output_diagram(outdiag);
  outdiag.close();

  return 0;
}
