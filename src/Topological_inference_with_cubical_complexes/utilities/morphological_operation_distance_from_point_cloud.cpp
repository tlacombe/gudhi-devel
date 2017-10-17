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
#include <gudhi/Morphological_operations_cubical_complex.h>

// standard stuff
#include <iostream>
#include <sstream>
#include <vector>
#include <sstream>

#include <gudhi/Points_off_io.h>
using namespace Gudhi;
using namespace Gudhi::Cubical_complex;


int main(int argc, char** argv) 
{
	std::cout << "This program takes as an input a point cloud. By using topological_inference it creates a filtered\
	cubical complex based on it. Later a tresholding of this function is performed, and a new complex with a filtraiton\
	equal to the distance from the chosen sublevel set of function is computed. At the end, the persistent homology of this\
	complex is computed." << std::endl;
	std::cout << "The program take the following input: " << std::endl;     
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
	std::cout << "(2i-1) the positive interger k, such that the considered\
	function is the distance to the k-th nearest neighbor.\n";
	std::cout << "(2i) A real number x, which will indicate the treshold level." << std::endl;
	std::cout << "(2i+1) An integer from a range {0,1,2}. Setting it to 0 will invoke dilation, to 1, erosion and to 2 both diation and erosion." << std::endl;
	std::cout << "When the operation is performed, the program will compute persistent homology of the obtained complex.\n";

    const char* filename = argv[1];
    std::vector< std::vector<double> > point_cloud_;  
  
    Gudhi::Points_off_reader< std::vector<double> > off_reader( filename );
  
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
    
    double treshold_value = atof( argv[3+3*dimension+1] );
    std::cout << "The treshold will be made at the funtion value : " << treshold_value << std::endl;
    
    int operation_type = atoi( argv[3+3*dimension+2] );
    std::cout << "The operation to be performed : ";
    if ( operation_type == 0 )std::cout << "dilation." << std::endl;
	if ( operation_type == 1 )std::cout << "erosion." << std::endl;
	if ( operation_type == 2 )std::cout << "both erosion and dilation." << std::endl;
	if ( (operation_type!=0)&&(operation_type!=1)&&(operation_type!=2) ) 
	{
		std::cout << "Wrong operation type, the program will now terminate.\n";
		return 1;
	}
  
    Gudhi::Topological_inference_with_cubical_complexes::Euclidan_distance_squared eu;
    Gudhi::Topological_inference_with_cubical_complexes::Distance_to_k_th_closest_point<Gudhi::Topological_inference_with_cubical_complexes::Euclidan_distance_squared> 
    f( off_reader.get_point_cloud() ,eu ,  number_of_nearest_neighbors );
  
    typedef Gudhi::Cubical_complex::Bitmap_cubical_complex_base<double> Bitmap_cubical_complex_base;
    typedef Gudhi::Cubical_complex::Bitmap_cubical_complex<Bitmap_cubical_complex_base> Bitmap_cubical_complex;
    typedef Gudhi::Topological_inference_with_cubical_complexes::Topological_inference< Bitmap_cubical_complex , double ,   
    Gudhi::Topological_inference_with_cubical_complexes::Distance_to_k_th_closest_point<Gudhi::Topological_inference_with_cubical_complexes::Euclidan_distance_squared> > topological_inference;
    topological_inference b( coorfinates_of_grid , resolution_of_a_grid , f );
	
    
    typedef Gudhi::Topological_inference_with_cubical_complexes::Filtration_below_certain_value<double> Predictor_type;
    typedef Gudhi::Topological_inference_with_cubical_complexes::Morphological_operations_cubical_complex<topological_inference,Predictor_type> MOCC;
    
    Predictor_type pred(treshold_value);    
    MOCC mor(&b,pred);
    if ( operation_type == 0 )mor.dilation( 1 , Gudhi::Topological_inference_with_cubical_complexes::considered_neighberhoods::all );
    if ( operation_type == 1 )mor.erosion( 1 , Gudhi::Topological_inference_with_cubical_complexes::considered_neighberhoods::all );
    if ( operation_type == 2 )mor.both_erosion_and_dilation( 1 , Gudhi::Topological_inference_with_cubical_complexes::considered_neighberhoods::all );
	
	
	typedef Gudhi::persistent_cohomology::Field_Zp Field_Zp;
	typedef Gudhi::persistent_cohomology::Persistent_cohomology<topological_inference, Field_Zp> Persistent_cohomology;
	Persistent_cohomology pcoh(b, true);
    pcoh.init_coefficients(2);  
	pcoh.compute_persistent_cohomology(0);

	std::stringstream ss;
	ss << filename << "_persistence";
	std::ofstream out(ss.str().c_str());
	pcoh.output_diagram(out);
	out.close();
		
	return 0;
}  
  
	
