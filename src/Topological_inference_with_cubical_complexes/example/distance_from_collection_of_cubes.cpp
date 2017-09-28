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


// for persistence algorithm
#include <gudhi/reader_utils.h>
#include <gudhi/Bitmap_cubical_complex.h>
#include <gudhi/Bitmap_cubical_complex_base.h>
#include <gudhi/Persistent_cohomology.h>
#include <gudhi/Distance_from_collection_of_cubes.h>
#include <gudhi/reader_utils.h>
#include <gudhi/Bitmap_cubical_complex.h>
#include <gudhi/Persistent_cohomology.h>


// standard stuff
#include <iostream>
#include <sstream>
#include <vector>

using namespace Gudhi;
using namespace Gudhi::Cubical_complex;

tutaj dodac jakies sensowne wejscie. 


int main(int argc, char** argv) 
{
	typedef Gudhi::Cubical_complex::Bitmap_cubical_complex_base<double> Bitmap_cubical_complex_base;
    typedef Gudhi::Cubical_complex::Bitmap_cubical_complex<Bitmap_cubical_complex_base> Bitmap_cubical_complex;
	
	std::cerr << filenames[file_no] << std::endl;
	//Bitmap_cubical_complex* b = construct_cubical_complex_and_set_up_the_filtration_to_distance_from_selected_cubes< Bitmap_cubical_complex >( cubes_to_set , sizes );
	std::vector< std::vector<unsigned> > aaa = read_Mao_file( filenames[file_no] );
	Bitmap_cubical_complex* b = Topological_inference_with_cubical_complexes::construct_cubical_complex_and_set_up_the_filtration_to_distance_from_selected_cubes< Bitmap_cubical_complex >(aaa);//( cubes_to_set , sizes );

	typedef Gudhi::persistent_cohomology::Field_Zp Field_Zp;
	typedef Gudhi::persistent_cohomology::Persistent_cohomology<Bitmap_cubical_complex, Field_Zp> Persistent_cohomology;
	Persistent_cohomology pcoh(*b, true);
	pcoh.init_coefficients(2);  // initializes the coefficient field for homology
	pcoh.compute_persistent_cohomology(0);

	std::stringstream ss;
	ss << filenames[file_no] << "_persistence";
	std::ofstream out(ss.str().c_str());
	pcoh.output_diagram(out);
	out.close();
	delete b;
	
		
	return 0;
}  
  
