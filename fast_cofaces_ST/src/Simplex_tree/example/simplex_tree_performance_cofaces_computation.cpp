/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Clément Maria
 *
 *    Copyright (C) 2014  INRIA Sophia Antipolis-Méditerranée (France)
 *
 *    This program is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *
 *    This program is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <gudhi/Rips_complex.h>
#include <gudhi/distance_functions.h>
#include <gudhi/Simplex_tree.h>
#include <gudhi/Persistent_cohomology.h>
#include <gudhi/Points_off_io.h>

#include <string>
#include <vector>
#include <limits>  // infinity

using namespace Gudhi;

typedef int Vertex_handle;
typedef double Filtration_value;
using Rips_complex = Gudhi::rips_complex::Rips_complex<Filtration_value>;
using Point = std::vector<double>;
typedef Simplex_tree<Simplex_tree_options_zigzag_persistence> opt_ST;
typedef Simplex_tree<>                                        brut_ST;

/* Compare the performance of the construction of Rips complexes and 
 * the computation of cofaces for simplex trees with Options 
 * Simplex_tree_options_zigzag_persistence and 
 * Simplex_tree_options_full_featured.
 *
 * The major difference between the two options is that in the first one, 
 * Nodes in the tree store two additionnal pointers, which allow for a fast 
 * cofaces computation.
 *
 * As a consequence, the construction of the same Rips complex with the first 
 * options is slightly slower than with the second, because of the time taken 
 * for allocating the pointers. However, the computation of cofaces becomes 
 * much faster on large simplicial complexes where simplces have small stars.
 */
int main(int argc, char * argv[]) {
  if(argc != 4) {
    std::cerr << "Usage: " << argv[0]
        << " path_to_file_points Rips_threshold max_dim \n";
    return 0;
  }
  std::string off_file_points = argv[1];
  Filtration_value threshold = atof(argv[2]);
  int dim_max = atoi(argv[3]);

  clock_t start, end;
  {
    // Construct the Simplex Tree
    std::cout << "---- Simplex tree with FAST cofaces option: \n";
    Gudhi::Points_off_reader<Point> off_reader(off_file_points);
    Rips_complex rips_complex_from_file(off_reader.get_point_cloud(), threshold, Euclidean_distance());

    start = clock();
    opt_ST  opt_st;
    rips_complex_from_file.create_complex(opt_st, dim_max);
    end = clock();
    std::cout << "     Construct Rips complex in "
        << static_cast<double>(end - start) / CLOCKS_PER_SEC << " s. \n";

    std::cout << "     Information of the Simplex Tree: " << std::endl;
    std::cout << "        Number of vertices = " << opt_st.num_vertices() << "\n";
    std::cout << "        Number of simplices = " << opt_st.num_simplices() << std::endl;

    std::cout << "     Compute the star of every simplex in the Simplex_tree:" 
              << std::endl;
    start = clock();
    unsigned long count = 0;
    for (auto f_simplex : opt_st.filtration_simplex_range()) { //for every simplex
      for(auto c_simplex : opt_st.cofaces_simplex_range(f_simplex,0)) 
      {++count;}//traverse the star
    }
    end = clock();
    std::cout << "     Information of the cofaces computation: " << std::endl;
    std::cout << "        Average star size = " << (double)count / (double)opt_st.num_simplices() << std::endl;
    std::cout << "        Total time for computation = " << static_cast<double>(end - start) / CLOCKS_PER_SEC << " s. \n";
    std::cout << "        Average time per simplex = " << (static_cast<double>(end - start) / CLOCKS_PER_SEC) / (double)opt_st.num_simplices() << " s. \n"; 
  }
    std::cout << std::endl;
  {
    // Construct the Simplex Tree
    std::cout << "---- Simplex tree with SLOW cofaces option: \n";
    Gudhi::Points_off_reader<Point> off_reader(off_file_points);
    Rips_complex rips_complex_from_file(off_reader.get_point_cloud(), threshold, Euclidean_distance());

    start = clock();
    brut_ST  brut_st;
    rips_complex_from_file.create_complex(brut_st, dim_max);
    end = clock();
    std::cout << "     Construct Rips complex in "
        << static_cast<double>(end - start) / CLOCKS_PER_SEC << " s. \n";

    std::cout << "     Information of the Simplex Tree: " << std::endl;
    std::cout << "        Number of vertices = " << brut_st.num_vertices() << "\n";
    std::cout << "        Number of simplices = " << brut_st.num_simplices() << std::endl;

    std::cout << "     Compute the star of every simplex in the Simplex_tree:" 
              << std::endl;
    start = clock();
    unsigned long count = 0;
    for (auto f_simplex : brut_st.filtration_simplex_range()) { //for every simplex
      for(auto c_simplex : brut_st.cofaces_simplex_range(f_simplex,0)) 
      {++count;}//traverse the star
    }
    end = clock();
    std::cout << "     Information of the cofaces computation: " << std::endl;
    std::cout << "        Average star size = " << (double)count / (double)brut_st.num_simplices() << std::endl;
    std::cout << "        Total time for computation = " << static_cast<double>(end - start) / CLOCKS_PER_SEC << " s. \n";
    std::cout << "        Average time per simplex = " << (static_cast<double>(end - start) / CLOCKS_PER_SEC) / (double)brut_st.num_simplices() << " s. \n"; 
  }

  return 0;
}
