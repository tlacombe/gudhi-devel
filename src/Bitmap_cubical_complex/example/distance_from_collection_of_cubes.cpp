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

/*
void write( std::vector<unsigned>& counter , Bitmap_cubical_complex_base<double>& b )
{
	size_t pos = b.give_position_of_top_dimensional_cell( counter );
	
	//std::cerr << "pos : " << pos << std::endl;
	
	std::vector<unsigned> couter_back = b.compute_counter_for_top_dimensional_cell( pos );
	for ( size_t i = 0 ; i != couter_back.size() ; ++i )
	{
		std::cerr << couter_back[i] << std::endl;	
	}
	
	
	std::vector<size_t> neighs = b.give_neighbouring_top_dimensional_cells( pos );
	for ( size_t i = 0 ; i != neighs.size() ; ++i )
	{
		std::cerr << neighs[i] << std::endl;	
	}
}*/

std::vector< std::vector<unsigned> > read_Mao_file( const char* filename )
{	
	bool dbg = false;
	
	std::ifstream in;
	in.open( filename );
	if ( !in.good() )
	{
		std::cerr << "The file do not exist \n";
		throw "The file do not exist \n";
	}
	
	double time_step;
	in >> time_step;
	size_t number_of_cubes;
	in >> number_of_cubes;
	
	std::vector< std::vector<unsigned> > result;
	result.reserve( number_of_cubes );
	while ( number_of_cubes )
	{
		unsigned x,y,z;
		in >> x >> y >> z;
		std::vector<unsigned> pt(3);
		pt[0] = x;
		pt[1] = y;
		pt[2] = z;
		result.push_back(pt); 
		if ( dbg )
		{
			std::cerr << x << " " << y << " " << z << std::endl;
			getchar();
		}
		--number_of_cubes;
	}
	return result;
}

int main(int argc, char** argv) 
{
	/*
	std::vector< unsigned > sizes;
	sizes.push_back(2);
	sizes.push_back(2);
	sizes.push_back(2);
	
	
	
	Bitmap_cubical_complex_base<double> b( sizes , data );
	
	std::vector<unsigned> counter(3);
	counter[0] = 1;
	counter[1] = 1;
	counter[2] = 1;
	write( counter , b );
	*/
	/*
	std::vector< double > data;
	data.push_back(1);
	data.push_back(2);
	data.push_back(3);
	data.push_back(4);
	data.push_back(5);
	data.push_back(6);
	data.push_back(7);
	data.push_back(8);
	
	std::vector< std::vector<unsigned> > cubes_to_set;
	{
		std::vector<unsigned> a;
		a.push_back(0);
		a.push_back(0);
		cubes_to_set.push_back( a );
	}
	{
		std::vector<unsigned> a;
		a.push_back(1);
		a.push_back(0);
		cubes_to_set.push_back( a );
	}
	{
		std::vector<unsigned> a;
		a.push_back(2);
		a.push_back(0);
		cubes_to_set.push_back( a );
	}
	{
		std::vector<unsigned> a;
		a.push_back(0);
		a.push_back(1);
		cubes_to_set.push_back( a );
	}
	{
		std::vector<unsigned> a;
		a.push_back(0);
		a.push_back(2);
		cubes_to_set.push_back( a );
	}
	{
		std::vector<unsigned> a;
		a.push_back(1);
		a.push_back(2);
		cubes_to_set.push_back( a );
	}
	{
		std::vector<unsigned> a;
		a.push_back(2);
		a.push_back(1);
		cubes_to_set.push_back( a );
	}
	{
		std::vector<unsigned> a;
		a.push_back(2);
		a.push_back(2);
		cubes_to_set.push_back( a );
	}

	
	std::vector< unsigned > sizes;
	sizes.push_back(3);
	sizes.push_back(3);*/
	
	
	
	std::vector< char* > filenames;
	filenames.push_back((char*)"ZmTIMp0020/ZmTIMp0020t00_njiang_2015-08-31_18-29-04_rootwork.out");
	filenames.push_back((char*)"ZmTIMp0020/ZmTIMp0020t01_njiang_2015-08-31_18-29-04_rootwork.out");
	filenames.push_back((char*)"ZmTIMp0020/ZmTIMp0020t02_njiang_2015-08-31_18-29-05_rootwork.out");
	filenames.push_back((char*)"ZmTIMp0020/ZmTIMp0020t03_njiang_2015-08-31_18-29-05_rootwork.out");
	filenames.push_back((char*)"ZmTIMp0020/ZmTIMp0020t04_njiang_2015-08-31_18-29-05_rootwork.out");
	filenames.push_back((char*)"ZmTIMp0020/ZmTIMp0020t05_njiang_2015-08-31_18-29-05_rootwork.out");
	filenames.push_back((char*)"ZmTIMp0020/ZmTIMp0020t06_njiang_2015-08-31_18-29-05_rootwork.out");
	filenames.push_back((char*)"ZmTIMp0020/ZmTIMp0020t07_njiang_2015-08-31_18-29-05_rootwork.out");
	filenames.push_back((char*)"ZmTIMp0020/ZmTIMp0020t08_njiang_2015-08-31_18-29-05_rootwork.out");
	filenames.push_back((char*)"ZmTIMp0020/ZmTIMp0020t09_njiang_2015-08-31_18-29-05_rootwork.out");
	filenames.push_back((char*)"ZmTIMp0020/ZmTIMp0020t10_njiang_2015-08-31_18-29-04_rootwork.out");
	filenames.push_back((char*)"ZmTIMp0020/ZmTIMp0020t11_njiang_2015-08-31_18-29-04_rootwork.out");
	filenames.push_back((char*)"ZmTIMp0020/ZmTIMp0020t12_njiang_2015-08-31_18-29-04_rootwork.out");
	filenames.push_back((char*)"ZmTIMp0020/ZmTIMp0020t13_njiang_2015-08-31_18-29-05_rootwork.out");
	filenames.push_back((char*)"ZmTIMp0020/ZmTIMp0020t14_njiang_2015-08-31_18-29-05_rootwork.out");
	filenames.push_back((char*)"ZmTIMp0020/ZmTIMp0020t15_njiang_2015-08-31_18-29-05_rootwork.out");
	filenames.push_back((char*)"ZmTIMp0020/ZmTIMp0020t16_njiang_2015-08-31_18-29-05_rootwork.out");
	filenames.push_back((char*)"ZmTIMp0020/ZmTIMp0020t17_njiang_2015-08-31_18-29-05_rootwork.out");
	filenames.push_back((char*)"ZmTIMp0020/ZmTIMp0020t18_njiang_2015-08-31_18-29-05_rootwork.out");
	filenames.push_back((char*)"ZmTIMp0020/ZmTIMp0020t19_njiang_2015-08-31_18-29-05_rootwork.out");
	filenames.push_back((char*)"ZmTIMp0020/ZmTIMp0020t20_njiang_2015-09-01_10-51-31_rootwork.out");
	filenames.push_back((char*)"ZmTIMp0020/ZmTIMp0020t21_njiang_2015-08-31_18-29-04_rootwork.out");
	filenames.push_back((char*)"ZmTIMp0020/ZmTIMp0020t22_njiang_2015-08-31_18-29-04_rootwork.out");
	filenames.push_back((char*)"ZmTIMp0020/ZmTIMp0020t23_njiang_2015-08-31_18-29-04_rootwork.out");
	filenames.push_back((char*)"ZmTIMp0020/ZmTIMp0020t24_njiang_2015-08-31_18-29-05_rootwork.out");
	filenames.push_back((char*)"ZmTIMp0020/ZmTIMp0020t25_njiang_2015-08-31_18-29-05_rootwork.out");
	filenames.push_back((char*)"ZmTIMp0020/ZmTIMp0020t26_njiang_2015-09-01_10-51-32_rootwork.out");
	filenames.push_back((char*)"ZmTIMp0020/ZmTIMp0020t27_njiang_2015-08-31_18-29-05_rootwork.out");
	filenames.push_back((char*)"ZmTIMp0020/ZmTIMp0020t28_njiang_2015-08-31_18-29-05_rootwork.out");
	filenames.push_back((char*)"ZmTIMp0020/ZmTIMp0020t29_njiang_2015-08-31_18-29-05_rootwork.out");
	filenames.push_back((char*)"ZmTIMp0020/ZmTIMp0020t30_njiang_2015-08-31_18-29-04_rootwork.out");
	filenames.push_back((char*)"ZmTIMp0020/ZmTIMp0020t31_njiang_2015-08-31_18-29-04_rootwork.out");
	filenames.push_back((char*)"ZmTIMp0020/ZmTIMp0020t32_njiang_2015-08-31_18-29-04_rootwork.out");
	filenames.push_back((char*)"ZmTIMp0020/ZmTIMp0020t33_njiang_2015-09-01_10-51-31_rootwork.out");
	filenames.push_back((char*)"ZmTIMp0020/ZmTIMp0020t34_njiang_2015-08-31_18-29-05_rootwork.out");
	filenames.push_back((char*)"ZmTIMp0020/ZmTIMp0020t35_njiang_2015-08-31_18-29-05_rootwork.out");
	filenames.push_back((char*)"ZmTIMp0020/ZmTIMp0020t36_njiang_2015-08-31_18-29-05_rootwork.out");
	filenames.push_back((char*)"ZmTIMp0020/ZmTIMp0020t37_njiang_2015-08-31_18-29-05_rootwork.out");
	filenames.push_back((char*)"ZmTIMp0020/ZmTIMp0020t38_njiang_2015-08-31_18-29-05_rootwork.out");
	filenames.push_back((char*)"ZmTIMp0020/ZmTIMp0020t39_njiang_2015-08-31_18-29-05_rootwork.out");
	filenames.push_back((char*)"ZmTIMp0020/ZmTIMp0020t40_njiang_2015-09-01_10-51-31_rootwork.out");
	filenames.push_back((char*)"ZmTIMp0020/ZmTIMp0020t41_njiang_2015-08-31_18-29-04_rootwork.out");

	 
	
	
	typedef Gudhi::Cubical_complex::Bitmap_cubical_complex_base<double> Bitmap_cubical_complex_base;
    typedef Gudhi::Cubical_complex::Bitmap_cubical_complex<Bitmap_cubical_complex_base> Bitmap_cubical_complex;
	
	
	for ( size_t file_no = 0 ; file_no != filenames.size() ; ++file_no )
	{
		std::cerr << filenames[file_no] << std::endl;
		//Bitmap_cubical_complex* b = construct_cubical_complex_and_set_up_the_filtration_to_distance_from_selected_cubes< Bitmap_cubical_complex >( cubes_to_set , sizes );
		std::vector< std::vector<unsigned> > aaa = read_Mao_file( filenames[file_no] );
		Bitmap_cubical_complex* b = construct_cubical_complex_and_set_up_the_filtration_to_distance_from_selected_cubes< Bitmap_cubical_complex >(aaa);//( cubes_to_set , sizes );
	
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
	}
  //std::cout << "Result in file: " << ss.str().c_str() << "\n";
	
	
	return 0;
}  
  
