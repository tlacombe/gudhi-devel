/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Pawel Dlotko
 *
 *    Copyright (C) 2017 Swansea University
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
#include <gudhi/Off_reader.h>

#include <gudhi/Topological_inference.h>
#include <gudhi/Bitmap_cubical_complex.h>
#include <gudhi/Persistent_cohomology.h>

// standard stuff
#include <iostream>
#include <sstream>
#include <vector>
#include <sstream>

#include <gudhi/Points_off_io.h>

using Point = std::vector<double>;
using Points_off_reader = Gudhi::Points_off_reader<Point>;


int main(int argc, char** argv) 
{	
	
  std::cout << "The parameters of the program are: \n";
  std::cout << "(1) A file in an OFF format with points coordinates, \n";
  std::cout << "(2) Dimension of a space, \n";
  std::cout << "(3) Minimum of a grid in first direction, \n";
  std::cout << "(4) Maximum of a grin in first direction, \n";
  std::cout << " ... ,\n";
  std::cout << "(i) Minimum of a grid in last direction, \n";
  std::cout << "(i+1) Maximum of a grin in last direction, \n";
  std::cout << "(i+2) resolution of a grid in the first direction,\n";
  std::cout << " ... ,\n";
  std::cout << "(2i-2) resolution of a grid in the last direction.\n";
  std::cout << "and finally the positive interger k, such that the considered\
   function is the distance to the k-th nearest neighbor.\n";

  int p = 2;
  double min_persistence = 0;

  const char* filename = argv[1];
  std::vector< std::vector<double> > point_cloud_;  
  
  Points_off_reader off_reader( argv[1] );
  
  int dimension = atoi( argv[2] );
  std::cout << "The dimension is : " << dimension << std::endl;
  
  std::vector< std::pair< double, double > > coorfinates_of_grid;
  for ( int i = 0 ; i != dimension ; ++i )
  {
	  int min_ = atof( argv[3+2*i] );
	  int max_ = atof( argv[3+2*i+1] );
	  coorfinates_of_grid.push_back( std::make_pair(min_,max_) );
	  std::cout << "Coordinates in direcion number : " << i << " are : " << min_ << " and " << max_ << std::endl;
  }	
  
  std::vector< unsigned > resolution_of_a_grid;
  for ( int i = 0 ; i != dimension ; ++i )
  {	 
	  size_t resolution_in_this_direction = (size_t)atoi( argv[3+2*dimension+i] );
	  resolution_of_a_grid.push_back( resolution_in_this_direction );
	  std::cout << "Resolution of a grid in direcion number : " << i << " is : " << resolution_in_this_direction << std::endl;
  }	
  
  
  unsigned number_of_nearest_neighbors = (unsigned)atoi( argv[3+3*dimension] );
  std::cout << "We will compute a distance to the " << number_of_nearest_neighbors << "-th nearest neighbor.\n";
  
  
  Gudhi::Topological_inference_with_cubical_complexes::Euclidan_distance_squared eu;
  Gudhi::Topological_inference_with_cubical_complexes::Distance_to_k_th_closest_point<Gudhi::Topological_inference_with_cubical_complexes::Euclidan_distance_squared> 
  f( off_reader.get_point_cloud() ,eu ,  number_of_nearest_neighbors );
  
  typedef Gudhi::Cubical_complex::Bitmap_cubical_complex_periodic_boundary_conditions_base<double> Periodic_bitmap_cubical_complex_base;
  typedef Gudhi::Cubical_complex::Bitmap_cubical_complex<Periodic_bitmap_cubical_complex_base> Periodic_bitmap_cubical_complex;
  typedef Gudhi::Topological_inference_with_cubical_complexes::Topological_inference< Periodic_bitmap_cubical_complex , double ,   
  Gudhi::Topological_inference_with_cubical_complexes::Distance_to_k_th_closest_point<Gudhi::Topological_inference_with_cubical_complexes::Euclidan_distance_squared> > topological_inference;
  
  typedef Gudhi::persistent_cohomology::Field_Zp Field_Zp;
  typedef Gudhi::persistent_cohomology::Persistent_cohomology<topological_inference, Field_Zp> Persistent_cohomology;

  topological_inference b( coorfinates_of_grid , resolution_of_a_grid , f );
  //b.write_to_file_Perseus_format("perse");

  // Compute the persistence diagram of the complex
  Persistent_cohomology pcoh(b);
  pcoh.init_coefficients(p);  // initializes the coefficient field for homology
  pcoh.compute_persistent_cohomology(min_persistence);
  
  std::stringstream ss;
  ss << filename << "_pers";

  std::ofstream out( ss.str().c_str() );
  pcoh.output_diagram(out);
  out.close();
  
  return 0;
}
