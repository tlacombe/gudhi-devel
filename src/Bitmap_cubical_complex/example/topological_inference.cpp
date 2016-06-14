/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Pawel Dlotko
 *
 *    Copyright (C) 2015  INRIA Saclay (France)
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

#include <gudhi/Topological_inference.h>
#include <gudhi/Bitmap_cubical_complex.h>
#include <gudhi/Persistent_cohomology.h>

// standard stuff
#include <iostream>
#include <sstream>
#include <vector>

int main(int argc, char** argv) {
  std::cout << "This program computes persistent homology, by using bitmap_cubical_complex class, of cubical " <<
      "complexes provided in text files in Perseus style (the only numbered in the first line is a dimension D of a" <<
      "bitmap. In the lines I between 2 and D+1 there are numbers of top dimensional cells in the direction I. Let " <<
      "N denote product of the numbers in the lines between 2 and D. In the lines D+2 to D+2+N there are " <<
      "filtrations of top dimensional cells. We assume that the cells are in the lexicographical order. See " <<
      "CubicalOneSphere.txt or CubicalTwoSphere.txt for example.\n" << std::endl;

  int p = 2;
  double min_persistence = 0;
  
  
  std::vector< std::vector<double> > point_cloud_;
  
  //point_cloud_ = read_points_from_file<double>( "circle" );
  
  point_cloud_ = read_points_from_file<double>( "2000_random_points_on_3D_Torus.csv" );
  /*
  std::vector<double> point1;
  point1.push_back( 0 );
  point1.push_back( 1);
  std::vector<double> point2;
  point2.push_back( 1 );
  point2.push_back( 2 );
  std::vector<double> point3;
  point3.push_back( 2 );
  point3.push_back( 1 );
  point_cloud_.push_back( point1 );
  point_cloud_.push_back( point2 );
  point_cloud_.push_back( point3 );
  */
  
  Distance_to_closest_point f( point_cloud_ );
  
  typedef Gudhi::Cubical_complex::Bitmap_cubical_complex_base<double> Bitmap_cubical_complex_base;
  typedef Gudhi::Cubical_complex::Bitmap_cubical_complex<Bitmap_cubical_complex_base> Bitmap_cubical_complex;
  typedef Gudhi::Cubical_complex::Topological_inference< Bitmap_cubical_complex , double , Distance_to_closest_point > topological_inference;
  
  typedef Gudhi::persistent_cohomology::Field_Zp Field_Zp;
  typedef Gudhi::persistent_cohomology::Persistent_cohomology<Bitmap_cubical_complex, Field_Zp> Persistent_cohomology;


  std::vector< std::pair< double, double > > coorfinates_of_grid;
  coorfinates_of_grid.push_back( std::make_pair(-4.0,4.0) );
  coorfinates_of_grid.push_back( std::make_pair(-4.0,4.0) );
  coorfinates_of_grid.push_back( std::make_pair(-2.0,2.0) );
  std::vector< unsigned > resolution_of_a_grid(3);
  resolution_of_a_grid[0] = 100;
  resolution_of_a_grid[1] = 100;
  resolution_of_a_grid[2] = 50;

/*
  std::vector< std::pair< double, double > > coorfinates_of_grid;
  coorfinates_of_grid.push_back( std::make_pair(0,5.0) );
  coorfinates_of_grid.push_back( std::make_pair(0.0,5.0) );  
  std::vector< unsigned > resolution_of_a_grid(2);
  resolution_of_a_grid[0] = 10;
  resolution_of_a_grid[1] = 10;
*/
  
  
  
  
  
  
   
  topological_inference b( coorfinates_of_grid , resolution_of_a_grid , f );
  
  b.write_to_file_Perseus_format("perse");

  // Compute the persistence diagram of the complex
  Persistent_cohomology pcoh(b);
  pcoh.init_coefficients(p);  // initializes the coefficient field for homology
  pcoh.compute_persistent_cohomology(min_persistence);

  std::ofstream out("top_inference");
  pcoh.output_diagram(out);
  out.close();

  return 0;
}
