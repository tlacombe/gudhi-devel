/*    This file is part of the Gudhi Library. The Gudhi library 
 *    (Geometric Understanding in Higher Dimensions) is a generic C++ 
 *    library for computational topology.
 *
 *    Author(s):       Pawel Dlotko, Clément Maria
 *
 *    Copyright (C) 2014  INRIA Saclay, Sophia Antipolis-Méditerranée (France)
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

#include <gudhi/Bitmap_cubical_complex.h>
#include <gudhi/Compute_persistence_with_phat.h>

#include <vector>

using namespace Gudhi;
using namespace Gudhi::phat_interface;
using namespace std;

typedef int Vertex_handle;
typedef double Filtration_value;

typedef Gudhi::cubical_complex::Bitmap_cubical_complex_base<double> Bitmap_cubical_complex_base;
typedef Gudhi::cubical_complex::Bitmap_cubical_complex<Bitmap_cubical_complex_base> Bitmap_cubical_complex;

int main(int argc, char * argv[]) 
{
	//In this file we present example of usage of options in PHAT. Note that some of them may be non optimal for some applications. For the best ones for the problem at hand, please consult
	//the PHAT wiki and the related papers that are mentioned in Gudhi documentation.	
	
	//here we create 3 by 3 2-dimensional cubical complex
	std::vector< double > data(9);
	data[0] = data[1] = data[2] = data[3] = data[5] = data[6] = data[7] = data[8] = 0;
	data[4] = 10;
	std::vector< unsigned > sizes(2);
	sizes[0] = sizes[1] = 3;
	Bitmap_cubical_complex b( sizes , data );
	
	//This is a default constructor, i.e. we create a boundary matrix and use twist_reduction algorithm from Phat. Note that this option may be sub--optimal. 
	//Compute_persistence_with_phat< Bitmap_cubical_complex > phat(&b);
	
	//In this constructor we are dualizing the matrix:
	Compute_persistence_with_phat< Bitmap_cubical_complex > phat(&b,true);
	
	phat.compute_persistence_pairs();
	std::pair< std::vector< std::vector<double> >, std::vector< std::vector< std::pair<double, double> > > > persistence = phat.get_the_intervals();
	
	std::cout << "Here are the Betti numbers : \n";
	for ( size_t dim = 0 ; dim != persistence.first.size() ; ++dim )
	{
		std::cout << "Betti numbers in dimension " << dim << " : ";	
		for ( size_t i = 0 ; i != persistence.first[dim].size() ; ++i )
		{
			std::cout << persistence.first[dim][i] << " ";
		}
		std::cout << std::endl;
		
	}
	
	std::cout << "Here are the perssitence intervals : \n";
	for ( size_t dim = 0 ; dim != persistence.second.size() ; ++dim )
	{
		std::cout << "Persistence diagrams in dimension " << dim << " : ";	
		for ( size_t i = 0 ; i != persistence.second[dim].size() ; ++i )
		{
			std::cout << persistence.second[dim][i].first << " " <<  persistence.second[dim][i].second << std::endl;
		}
		std::cout << std::endl;
	}
	
	return 0;
}

