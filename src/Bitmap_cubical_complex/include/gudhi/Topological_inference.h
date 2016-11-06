/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Pawel Dlotko
 *
 *    Copyright (C) 2015  INRIA Sophia-Saclay (France)
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

#ifndef TOPOLOGICAL_INFERENCE_H_
#define TOPOLOGICAL_INFERENCE_H_

#include <gudhi/Bitmap_cubical_complex_base.h>
#include <gudhi/Bitmap_cubical_complex_periodic_boundary_conditions_base.h>
#include <gudhi/Bitmap_cubical_complex.h>
#include <gudhi/functions_for_topological_inference/functions_for_topological_inference.h>

#ifdef GUDHI_USE_TBB
#include <tbb/parallel_sort.h>
#endif

#include <limits>
#include <utility>  // for pair<>
#include <algorithm>  // for sort
#include <vector>
#include <numeric>  // for iota

namespace Gudhi {

namespace Cubical_complex {

/*
 *  T is Bitmap_cubical_complex<whatever>. K is the number type of a grid (typically double). F is a function that computes value on a grid. It acts from K^n -> K
*/
template <typename T , typename K , typename F>
class Topological_inference : public T
{
public:
	/**
	 * Default constructor. 
	**/ 
	Topological_inference():T(){};
	
	/**
	 * Default constructor. 
	**/
	Topological_inference( const std::vector< std::pair<K , K> >& coorfinates_of_grid_ , const std::vector< unsigned >& resolution_of_a_grid_ , F& f  );
	
	/**
	 * Writing maximal cells to a file in Perseus format: 
	**/ 
	void write_to_file_Perseus_format( const char* filename );
	
protected:
	inline std::vector<size_t> compute_counter_for_this_maximal_cube( size_t maximal_cube_index )
	{
		bool dbg = false;
		if ( dbg )std::cout << "compute_counter_for_this_maximal_cube procedure for : " << maximal_cube_index << std::endl;
		
		std::vector<size_t> result( this->coorfinates_of_grid.size() , 0 );
		for ( size_t i = 0 ; i != this->coorfinates_of_grid.size() ; ++i )
		{
			result[i] = maximal_cube_index%this->resolution_of_a_grid[i];
			maximal_cube_index = maximal_cube_index/this->resolution_of_a_grid[i];
			if ( dbg )
			{
				std::cout << "result[" << i << "] = " << result[i] << std::endl;
			}			
		}
		return result;
	}
	
	inline std::vector< K > compute_center_of_cube_for_given_counter( const std::vector<size_t>& counter )
	{
		bool dbg = false;
		if ( dbg )
		{
			std::cerr << "Entering compute_center_of_cube_for_given_counter procedure \n";
			std::cout << "counter.size() : : " << counter.size() << std::endl;
			std::cout << "this->coorfinates_of_grid.size() : " << this->coorfinates_of_grid.size() << std::endl;
		}
		if ( counter.size() != this->coorfinates_of_grid.size() )throw "Wrong dimensionality of a counter in the procedure compute_center_of_cube_for_given_counter \n";
		std::vector< K > result( counter.size() , 0 );
		
		for ( size_t dim = 0 ; dim != counter.size() ; ++dim )
		{
			if ( counter[dim] >= this->resolution_of_a_grid[dim] )throw "The counter in some dimension extends dimensionality of a grid. The program will now terminate \n";
			
			K dx = (this->coorfinates_of_grid[dim].second - this->coorfinates_of_grid[dim].first)/this->resolution_of_a_grid[dim];
			result[dim] = dx*counter[dim] + this->coorfinates_of_grid[dim].first + dx/2;
		}	
		return result;
	}
	
	inline K compute_value_at_a_given_point( const std::vector< K >& point )
	{	
		return this->f( point );
	}
	
	//data:
	std::vector< std::pair< K , K > > coorfinates_of_grid;
	std::vector< unsigned > resolution_of_a_grid;
	F& f;
};



template <typename T , typename K , typename F>
Topological_inference<T,K,F>::Topological_inference( const std::vector< std::pair<K , K> >& coorfinates_of_grid_ , const std::vector< unsigned >& resolution_of_a_grid_ , F& f ):T(resolution_of_a_grid_), f(f)
{
	bool dbg = false;
	if ( dbg )
	{
		std::cerr << "Entering constructor of a Topological_inference object \n";
		std::cout << "coorfinates_of_grid_.size() : " << coorfinates_of_grid_.size() << std::endl;
		std::cout << "resolution_of_a_grid_.size() : " << resolution_of_a_grid_.size() << std::endl;
	}
	if ( coorfinates_of_grid_.size() != resolution_of_a_grid_.size() )throw "Incompatible sizes of coorfiantes of a grid, and the resoution of a grid in the constructore of Topological_inference. The program will now terminate \n";
	
	this->coorfinates_of_grid = coorfinates_of_grid_;
	this->resolution_of_a_grid = resolution_of_a_grid_;		
	
	size_t number_of_maximal_cubes = this->number_of_maximal_cubes();
	std::vector< K > values_on_maximal_cells( number_of_maximal_cubes , 0 );
	
	#pragma omp parallel for 
	for ( size_t i = 0 ; i < number_of_maximal_cubes ; ++i )
	{
		std::vector< size_t > counter = this->compute_counter_for_this_maximal_cube( i );
		if ( dbg )
		{
			std::cerr << "Here is the counter corresponding to the cube : " << i << std::endl;
			for ( size_t yy = 0 ; yy != counter.size() ; ++yy )
			{
				std::cout << counter[yy] << std::endl;
			}
			getchar();
		}
		std::vector< K > point = this->compute_center_of_cube_for_given_counter( counter );
		
		if ( dbg )
		{
			std::cerr << "Here is the the center of the considered cube (point) : " << i << std::endl;
			for ( size_t yy = 0 ; yy != point.size() ; ++yy )
			{
				std::cout << point[yy] << std::endl;
			}			
		}
		
		K value = this->compute_value_at_a_given_point( point );
		values_on_maximal_cells[i] = value;		
		if ( dbg )
		{
			std::cerr << "Value of a function at this point: " << value << std::endl;
			std::cin.ignore();
		}
	}
	
	std::vector<size_t> counter_v( coorfinates_of_grid_.size() , 0 );
	size_t i = 0;
		
	
	typename T::Top_dimensional_cells_iterator it = this->top_dimensional_cells_iterator_begin();
	while ( i < number_of_maximal_cubes )
	{		
		size_t index = (*it);
		this->data[index] = values_on_maximal_cells[i];					
		if ( dbg )
		{
			std::cout << "Maximal cube of an index : " << index << " gets the value : " << values_on_maximal_cells[i] << std::endl;
		}
		++i;
		++it;			
	}
	if ( dbg )std::cout << "Done with assigning values. Now will impose lower star filtration \n";
	this->impose_lower_star_filtration();
	this->initialize_simplex_associated_to_key();
		
	//std::ofstream out;
	//out.open("gnupl");
	//i = 0;
	//while ( i < number_of_maximal_cubes )
	//{		
	//	out << values_on_maximal_cells[i] << " ";					
	//	if ( i % resolution_of_a_grid_[0] == resolution_of_a_grid_[0]-1 )
	//	{
	//		out << std::endl;
	//	}
	//	++i;					
	//}
	//out.close();
}


template <typename T , typename K , typename F>
void Topological_inference<T,K,F>::write_to_file_Perseus_format( const char* filename )
{
	this->store_in_perseus_format( filename );
}

}  // namespace Cubical_complex

}  // namespace Gudhi

#endif  // BITMAP_CUBICAL_COMPLEX_H_
