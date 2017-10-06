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

class distnace_from_x_axis
{
public:
	distnace_from_x_axis(){}
	double operator()( const std::vector< double >& point )
	{
		double norm = 0;
		for ( size_t i = 1 ; i != point.size() ; ++i )
		{
			norm += point[i]*point[i];
		}
		norm = sqrt(norm);
		return norm;
	}
};


int main()
{
	//This example show how to do topological inference given a function. In this example we will be using the
	//distnace_from_x_axis function defined above. 
	
	
	//First we need to decide the rectagle that will encapsulate our cubical grid. In this case,
	//we are picking up a square [-1,1]^3
	std::vector< std::pair< double,double > > coorfinates_of_grid(3);
	coorfinates_of_grid[0] = coorfinates_of_grid[1] = coorfinates_of_grid[2] = std::pair<double,double>( -2.0 , 2.0 );
	
	//Later we need to decide on the resolution. Here our two dimensional grid will have 
	//100 maximal cubes resolution in each direction (giving 100*100=10000 grid points)
	std::vector< unsigned > resolution_of_a_grid(3);	
	resolution_of_a_grid[0] = resolution_of_a_grid[1] = resolution_of_a_grid[2] = 50;	
	
	//Here we need to decide the directions in which periodic boundary conditions are to be imposed:
	std::vector<bool> directions_in_which_periodic_b_cond_are_to_be_imposed(3);
	directions_in_which_periodic_b_cond_are_to_be_imposed[0] = true;
	directions_in_which_periodic_b_cond_are_to_be_imposed[1] = true;
	directions_in_which_periodic_b_cond_are_to_be_imposed[2] = true;
	
	distnace_from_x_axis dist;
	
	//A few typedefs to make the code easy to read:	
	//topological inference typedefs
    typedef Gudhi::Cubical_complex::Bitmap_cubical_complex_periodic_boundary_conditions_base<double> Periodic_bitmap_cubical_complex_base;
    typedef Gudhi::Cubical_complex::Bitmap_cubical_complex<Periodic_bitmap_cubical_complex_base> Periodic_bitmap_cubical_complex;
    typedef Gudhi::Topological_inference_with_cubical_complexes::Topological_inference< Periodic_bitmap_cubical_complex , double ,   
    distnace_from_x_axis > topological_inference;
  
    //computations of persistence typedefs
    typedef Gudhi::persistent_cohomology::Field_Zp Field_Zp;
    typedef Gudhi::persistent_cohomology::Persistent_cohomology<topological_inference, Field_Zp> Persistent_cohomology;

	//now we can create our topological inference object
    topological_inference b( coorfinates_of_grid , resolution_of_a_grid , dist , directions_in_which_periodic_b_cond_are_to_be_imposed );        

    //And compute its persistence diagram
    Persistent_cohomology pcoh(b,true);
    pcoh.init_coefficients(2); 
    
    //we want only the intervals of persistence >= 0.1
    pcoh.compute_persistent_cohomology(0.1);
    
    //and print it out to the screen. Note that thanks to the function which has global minimum (zero)
    //on the x axis, one of the infinite intervals in dimnesion 1 is born very close to zero, unlike
    //the other two of them.
    std::vector< std::tuple<size_t, size_t, int> > intervals = pcoh.get_persistent_pairs();   
    for ( size_t i = 0 ; i != intervals.size() ; ++i )
    {
		std::cout << "Dimension : " << b.get_dimension_of_a_cell(std::get<0>(intervals[i])) << " birth: " << b.filtration(std::get<0>(intervals[i])) 
		<< " death: " << b.filtration(std::get<1>(intervals[i])) << std::endl;
	}
	
	return 0;
}
