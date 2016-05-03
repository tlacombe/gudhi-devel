 /*    This file is part of the Gudhi Library. The Gudhi library
  *    (Geometric Understanding in Higher Dimensions) is a generic C++
  *    library for computational topology.
  *
  *    Author(s):       Pawel Dlotko
  *
  *    Copyright (C) 2016  INRIA Saclay (France)
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


#include <gudhi/reader_utils.h>

//cubical complex include
#include <gudhi/Bitmap_cubical_complex_base.h>
#include <gudhi/Bitmap_cubical_complex_periodic_boundary_conditions_base.h>
#include <gudhi/Bitmap_cubical_complex.h>
#include <gudhi/Compute_persistence_with_phat.h>

//Rips compelx and simpelx tree:
#include <gudhi/graph_simplicial_complex.h>
#include <gudhi/distance_functions.h>
#include <gudhi/Simplex_tree.h>
#include <boost/program_options.hpp>

//standard stuff
#include <iostream>
#include <sstream>
#include <vector>
#include <string>
#include <limits>

using namespace Gudhi;
using namespace Gudhi::Cubical_complex;
using namespace std;

typedef int Vertex_handle;
typedef double Filtration_value;


BOOST_AUTO_TEST_CASE(check_dimension) 
{
  typedef std::vector<double> Point_t;
  std::vector< Point_t > points;
  read_points("plane_circle", points);

  // Compute the proximity graph of the points
  Graph_t prox_graph = compute_proximity_graph(points, threshold , euclidean_distance<Point_t>);

  // Construct the Rips complex in a Simplex Tree
  typedef Simplex_tree<Simplex_tree_options_fast_persistence> ST;
  ST st;
  // insert the proximity graph in the simplex tree
  st.insert_graph(prox_graph);
  // expand the graph until dimension dim_max
  st.expansion(dim_max);
  
  // Sort the simplices in the order of the filtration
  st.initialize_filtration();

  std::cout << "The complex contains " << st.num_simplices() << " simplices \n";
  std::cout << "   and has dimension " << st.dimension() << " \n";
  
  Compute_persistence_with_phat< ST , double > phat(&st);


  //phat::persistence_pairs pairs = phat.compute_persistence_pairs_dualized_chunk_reduction();
  phat::persistence_pairs pairs = phat.compute_persistence_pairs_twist_reduction();
  
  std::pair< std::vector< std::vector<double> > , std::vector< std::vector< std::pair<double,double> > > > persistence = phat.get_the_intervals( pairs );
  	
	
	
  BOOST_CHECK(increasing.dimension() == 2);
}


