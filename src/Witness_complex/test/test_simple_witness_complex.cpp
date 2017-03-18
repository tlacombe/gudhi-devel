#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE "simple_witness_complex"
#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <gudhi/Simplex_tree.h>

#include <gudhi/Witness_complex.h>
#include <gudhi/Euclidean_witness_complex.h>
#include <gudhi/Witness_complex_new.h>
#include <gudhi/Witness_complex_cof.h>
#include <gudhi/Strong_witness_complex.h>
#include <gudhi/Euclidean_strong_witness_complex.h>


#include <gudhi/Kd_tree_search.h>

#include <iostream>
#include <vector>
#include <utility>

typedef Gudhi::Simplex_tree<> Simplex_tree;
// typedef typename Gudhi::Simplex_tree<>::Vertex_handle Vertex_handle;
// typedef std::vector< Vertex_handle > typeVectorVertex;
typedef CGAL::Epick_d<CGAL::Dynamic_dimension_tag> Kernel;
typedef typename Kernel::Point_d Point_d;
typedef std::vector<Point_d> Point_range;
typedef Gudhi::spatial_searching::Kd_tree_search<Kernel, Point_range> Kd_tree;
typedef Kd_tree::INS_range Nearest_landmark_range; 
typedef std::vector<Nearest_landmark_range> Nearest_landmark_table;
typedef Gudhi::witness_complex::Witness_complex_cof<Nearest_landmark_table> WitnessComplex;
typedef Gudhi::witness_complex::Strong_witness_complex<Nearest_landmark_table> StrongWitnessComplex;


/* All landmarks and witnesses are taken on the grid in the following manner.
   LWLWL  
   WW.WW  
   L...L  
   WW.WW  
   LWLWL  

   Witness complex consists of 8 vertices, 12 edges and 4 triangles
 */

BOOST_AUTO_TEST_CASE(simple_witness_complex) {
  using Nearest_landmark_range = std::vector<std::pair<std::size_t, double>>;
  using Nearest_landmark_table = std::vector<Nearest_landmark_range>;
  using Witness_complex = Gudhi::witness_complex::Witness_complex<Nearest_landmark_table>;
  using Simplex_tree = Gudhi::Simplex_tree<>;

  Simplex_tree stree;
  Nearest_landmark_table nlt;

  // Example contains 5 witnesses and 5 landmarks
  Nearest_landmark_range w0 = {std::make_pair(0, 0), std::make_pair(1, 1), std::make_pair(2, 2),
                               std::make_pair(3, 3), std::make_pair(4, 4)}; nlt.push_back(w0);
  Nearest_landmark_range w1 = {std::make_pair(1, 0), std::make_pair(2, 1), std::make_pair(3, 2),
                               std::make_pair(4, 3), std::make_pair(0, 4)}; nlt.push_back(w1);
  Nearest_landmark_range w2 = {std::make_pair(2, 0), std::make_pair(3, 1), std::make_pair(4, 2),
                               std::make_pair(0, 3), std::make_pair(1, 4)}; nlt.push_back(w2);
  Nearest_landmark_range w3 = {std::make_pair(3, 0), std::make_pair(4, 1), std::make_pair(0, 2),
                               std::make_pair(1, 3), std::make_pair(2, 4)}; nlt.push_back(w3);
  Nearest_landmark_range w4 = {std::make_pair(4, 0), std::make_pair(0, 1), std::make_pair(1, 2),
                               std::make_pair(2, 3), std::make_pair(3, 4)}; nlt.push_back(w4);

  // Weak witness complex: Euclidean version
  // EuclideanWitnessComplex eucl_witness_complex(landmarks,
  //                                              witnesses);
  // eucl_witness_complex.create_complex(complex, 0);

  // std::cout << "complex.num_simplices() = " << complex.num_simplices() << std::endl; 
  // BOOST_CHECK(complex.num_simplices() == 24);

  // eucl_witness_complex.create_complex(relaxed_complex, 8.01);

  // std::cout << "relaxed_complex.num_simplices() = " << relaxed_complex.num_simplices() << std::endl; 
  // BOOST_CHECK(relaxed_complex.num_simplices() == 239);
  // // The corner simplex {0,2,5,7} and its cofaces are missing.

  // Weak witness complex: non-Euclidean version
  WitnessComplex witness_complex(nearest_landmark_table);
  witness_complex.create_complex(complex_ne, 8.01);

  // std::cout << complex_ne << std::endl;
  
  std::cout << "complex.num_simplices() = " << complex_ne.num_simplices() << std::endl; 
  // BOOST_CHECK(complex_ne.num_simplices() == 24);

  //witness_complex.create_complex(relaxed_complex_ne, 8.01);

  // std::cout << "relaxed_complex.num_simplices() = " << relaxed_complex_ne.num_simplices() << std::endl; 
  // BOOST_CHECK(relaxed_complex_ne.num_simplices() == 239);
    
  
  // // Strong complex : Euclidean version
  // EuclideanStrongWitnessComplex eucl_strong_witness_complex(landmarks,
  //                                                           witnesses);

  // eucl_strong_witness_complex.create_complex(strong_relaxed_complex, 9.1);
  // eucl_strong_witness_complex.create_complex(strong_relaxed_complex2, 9.1, 2);
  
  // std::cout << "strong_relaxed_complex.num_simplices() = " << strong_relaxed_complex.num_simplices() << std::endl; 
  // BOOST_CHECK(strong_relaxed_complex.num_simplices() == 239);

  // std::cout << "strong_relaxed_complex2.num_simplices() = " << strong_relaxed_complex2.num_simplices() << std::endl;
  // BOOST_CHECK(strong_relaxed_complex2.num_simplices() == 92);


  // // Strong complex : non-Euclidean version
  // StrongWitnessComplex strong_witness_complex(nearest_landmark_table);

  // strong_witness_complex.create_complex(strong_relaxed_complex_ne, 9.1);
  // strong_witness_complex.create_complex(strong_relaxed_complex2_ne, 9.1, 2);
  
  // std::cout << "strong_relaxed_complex.num_simplices() = " << strong_relaxed_complex_ne.num_simplices() << std::endl; 
  // BOOST_CHECK(strong_relaxed_complex_ne.num_simplices() == 239);

  // std::cout << "strong_relaxed_complex2.num_simplices() = " << strong_relaxed_complex2_ne.num_simplices() << std::endl;
  // BOOST_CHECK(strong_relaxed_complex2_ne.num_simplices() == 92);


  // 8 vertices, 28 edges, 56 triangles
}
