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
#include <gudhi/Bitmap_cubical_complex.h>
#include <gudhi/Persistent_cohomology.h>

// standard stuff
#include <iostream>
#include <string>
#include <vector>
#include <cstddef>

int main() 
{
	/**
	* This example ilustrates how to set up the values of individual cubes and 
	* later to compute persistent homology of the obtained complex.  
	**/ 
	typedef Gudhi::cubical_complex::Bitmap_cubical_complex_base<double> Bitmap_cubical_complex_base;
	typedef Gudhi::cubical_complex::Bitmap_cubical_complex<Bitmap_cubical_complex_base> Bitmap_cubical_complex;
	typedef Gudhi::persistent_cohomology::Field_Zp Field_Zp;
	typedef Gudhi::persistent_cohomology::Persistent_cohomology<Bitmap_cubical_complex, Field_Zp> Persistent_cohomology;

	// __________________
	//|0|__0__|0|__0__|0|
	//| |     | |     | | 	
	//|0|  5  |2|  6  |0|
	//| |_____| |_____| |
	//|0|__3__|1|__4__|0|
	//| |     | |     | | 	
	//|0|  7  |2|  8  |0|
	//| |_____| |_____| |
	//|0|__0__|0|__0__|0|
	//

	//For this we first create 2 dimensional 2 by 2 cubical complex set with
	//all the filtration initially set to zero. 
	std::vector<unsigned> sizes(2,2);
	std::vector<double> data(4,0);
	Bitmap_cubical_complex b(sizes,data);
	
	//now we will iterate through all the cells. We will use counters for that. 
	//Here are the counters for the cells:
	//________________________
	//|04|__14__|24|__34__|44|
	//|  |      |  |      |  | 	
	//|03|  13  |23|  33  |43|
	//|  |______|  |______|  |
	//|02|__12__|22|__32__|42|
	//|  |      |  |      |  | 	
	//|01|  11  |21|  31  |41|
	//|  |______|  |______|  |	
	//|00|__10__|20|__30__|40|
	//	
	
	std::vector<unsigned> counter(2);
	counter[0] = 0;counter[1] = 0;
	b.get_cell_data( b.compute_position_in_bitmap(counter) ) = 0;
	counter[0] = 1;counter[1] = 0;
	b.get_cell_data( b.compute_position_in_bitmap(counter) ) = 0;
	counter[0] = 2;counter[1] = 0;
	b.get_cell_data( b.compute_position_in_bitmap(counter) ) = 0;
	counter[0] = 3;counter[1] = 0;
	b.get_cell_data( b.compute_position_in_bitmap(counter) ) = 0;
	counter[0] = 4;counter[1] = 0;	
	b.get_cell_data( b.compute_position_in_bitmap(counter) ) = 0;
	counter[0] = 0;counter[1] = 1;
	b.get_cell_data( b.compute_position_in_bitmap(counter) ) = 0;
	counter[0] = 1;counter[1] = 1;
	b.get_cell_data( b.compute_position_in_bitmap(counter) ) = 7;
	counter[0] = 2;counter[1] = 1;
	b.get_cell_data( b.compute_position_in_bitmap(counter) ) = 2;
	counter[0] = 3;counter[1] = 1;
	b.get_cell_data( b.compute_position_in_bitmap(counter) ) = 8;
	counter[0] = 4;counter[1] = 1;
	b.get_cell_data( b.compute_position_in_bitmap(counter) ) = 0;
	counter[0] = 0;counter[1] = 2;
	b.get_cell_data( b.compute_position_in_bitmap(counter) ) = 0;
	counter[0] = 1;counter[1] = 2;
	b.get_cell_data( b.compute_position_in_bitmap(counter) ) = 3;
	counter[0] = 2;counter[1] = 2;
	b.get_cell_data( b.compute_position_in_bitmap(counter) ) = 1;
	counter[0] = 3;counter[1] = 2;
	b.get_cell_data( b.compute_position_in_bitmap(counter) ) = 4;
	counter[0] = 4;counter[1] = 2;
	b.get_cell_data( b.compute_position_in_bitmap(counter) ) = 0;
	counter[0] = 0;counter[1] = 3;
	b.get_cell_data( b.compute_position_in_bitmap(counter) ) = 0;
	counter[0] = 1;counter[1] = 3;
	b.get_cell_data( b.compute_position_in_bitmap(counter) ) = 5;
	counter[0] = 2;counter[1] = 3;
	b.get_cell_data( b.compute_position_in_bitmap(counter) ) = 2;
	counter[0] = 3;counter[1] = 3;
	b.get_cell_data( b.compute_position_in_bitmap(counter) ) = 6;
	counter[0] = 4;counter[1] = 3;
	b.get_cell_data( b.compute_position_in_bitmap(counter) ) = 0;
	counter[0] = 0;counter[1] = 4;
	b.get_cell_data( b.compute_position_in_bitmap(counter) ) = 0;
	counter[0] = 1;counter[1] = 4;
	b.get_cell_data( b.compute_position_in_bitmap(counter) ) = 0;
	counter[0] = 2;counter[1] = 4;
	b.get_cell_data( b.compute_position_in_bitmap(counter) ) = 0;
	counter[0] = 3;counter[1] = 4;
	b.get_cell_data( b.compute_position_in_bitmap(counter) ) = 0;
	counter[0] = 4;counter[1] = 4;
	b.get_cell_data( b.compute_position_in_bitmap(counter) ) = 0;
	
	//now in order to compute persistence we need to call:
	b.initialize_simplex_associated_to_key();
	
	//One can vizualize the complex with this code:
	std::cout << "Here is the complex (note that due to the way it is displayed, it is up side down.\n";
	for ( size_t i = 0 ; i != b.size() ; ++i )
	{
		std::cout << b.get_cell_data(i) << " ";
		if ( i%5 == 4 )std::cout << std::endl;
	}

	// Compute the persistence diagram of the complex
	std::cout << "Here is the persistence of the complex.\n";
	Persistent_cohomology pcoh(b);
	int p = 11;
	double min_persistence = 0;
	pcoh.init_coefficients(p);  // initializes the coefficient field for homology
	pcoh.compute_persistent_cohomology(min_persistence);
	pcoh.output_diagram();
  
  return 0;
}
