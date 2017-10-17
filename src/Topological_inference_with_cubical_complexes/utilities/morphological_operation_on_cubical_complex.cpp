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
#include <gudhi/Morphological_operations_cubical_complex.h>
#include <gudhi/reader_utils.h>
#include <gudhi/Bitmap_cubical_complex.h>
#include <gudhi/Persistent_cohomology.h>


// standard stuff
#include <iostream>
#include <sstream>
#include <vector>

using namespace Gudhi;
using namespace Gudhi::Cubical_complex;


int main(int argc, char** argv) 
{
	std::cout << "This is a morphological operation on cubical complex program. The program take the following input: " << std::endl; 
	std::cout << "(1) A file with cubical complex," << std::endl;
	std::cout << "(2) A real number x, which will indicate the treshold level." << std::endl;
	std::cout << "(3) An integer from a range {0,1,2}. Setting it to 0 will invoke dilation, to 1, erosion and to 2 both diation and erosion." << std::endl;
	std::cout << "When the operation is performed, the program will compute persistent homology of the obtained complex.\n";
	
	if ( argc != 4 )
	{
		std::cout << "Wrong number of parameters, the program will now terminate. " << std::endl;
		return 1;
	}
	
	const char* filename = argv[1];
	double treshold_value = atof( argv[2] );
	int operation_type = atoi( argv[3] );
	
	std::cout << "Here are the parameters of the program : name of the file : " 
	<< filename << ", value of the treshold : " <<  treshold_value << ", operation_type : ";
	if ( operation_type == 0 )std::cout << "dilation." << std::endl;
	if ( operation_type == 1 )std::cout << "erosion." << std::endl;
	if ( operation_type == 2 )std::cout << "both erosion and dilation." << std::endl;
	
	typedef Gudhi::Cubical_complex::Bitmap_cubical_complex_base<double> Bitmap_cubical_complex_base;
    typedef Gudhi::Cubical_complex::Bitmap_cubical_complex<Bitmap_cubical_complex_base> Bitmap_cubical_complex;
    
    Bitmap_cubical_complex b(filename);
    
    typedef Gudhi::Topological_inference_with_cubical_complexes::Filtration_below_certain_value<double> Predictor_type;
    typedef Gudhi::Topological_inference_with_cubical_complexes::Morphological_operations_cubical_complex<Bitmap_cubical_complex,Predictor_type> MOCC;
    
    Predictor_type pred(treshold_value);    
    MOCC mor(&b,pred);
    if ( operation_type == 0 )mor.dilation( 1 , Gudhi::Topological_inference_with_cubical_complexes::considered_neighberhoods::all );
    if ( operation_type == 1 )mor.erosion( 1 , Gudhi::Topological_inference_with_cubical_complexes::considered_neighberhoods::all );
    if ( operation_type == 2 )mor.both_erosion_and_dilation( 1 , Gudhi::Topological_inference_with_cubical_complexes::considered_neighberhoods::all );
	
	
	typedef Gudhi::persistent_cohomology::Field_Zp Field_Zp;
	typedef Gudhi::persistent_cohomology::Persistent_cohomology<Bitmap_cubical_complex, Field_Zp> Persistent_cohomology;
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
  
