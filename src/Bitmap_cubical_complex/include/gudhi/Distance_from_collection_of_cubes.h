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

#ifndef DISTANCE_FROM_COLLECTION_OF_CUBES_
#define DISTANCE_FROM_COLLECTION_OF_CUBES_

#include <gudhi/Bitmap_cubical_complex_base.h>
#include <gudhi/Bitmap_cubical_complex_periodic_boundary_conditions_base.h>
#include <gudhi/Bitmap_cubical_complex.h>
#ifdef GUDHI_USE_TBB
#include <tbb/parallel_sort.h>
#endif

#include <limits>
#include <utility>  // for pair<>
#include <algorithm>  // for sort
#include <vector>
#include <numeric>  // for iota

namespace Gudhi 
{

namespace Cubical_complex 
{
/**
 * The function presented here construct a cubical complex with a filtration. The constructed filtration is integer--valued and the value indicate the Manhattan distance from a celected collection of cubes
 **/ 

template < typename Cubical_complex >
Cubical_complex* construct_cubical_complex_and_set_up_the_filtration_to_distance_from_selected_cubes( const std::vector< std::vector< unsigned > >& top_dimensional_cubes_that_are_in_the_set , 
																									  std::vector< unsigned > sizes = std::vector< unsigned >() )
{
	bool dbg = false;
	if ( sizes.empty() )
	{
		if ( dbg )std::cout << "the ranges of the cubical complex were not set by the user. We will set them now based on top_dimensional_cubes_that_are_in_the_set vector \n";
		if ( top_dimensional_cubes_that_are_in_the_set.size() == 0 )return new Cubical_complex;//in this case, we return an empty strucutre of the desired type. 
		sizes = std::vector< unsigned >( top_dimensional_cubes_that_are_in_the_set[0].size() , 0 );		
		//in this case, we should set up the sizes by ourselves. 
		for ( size_t i = 0 ; i != top_dimensional_cubes_that_are_in_the_set.size() ; ++i )
		{
			if ( top_dimensional_cubes_that_are_in_the_set[0].size() != top_dimensional_cubes_that_are_in_the_set[i].size() )
			{
				std::cerr << "Different embedding dimensions of cubes in the procedure construct_cubical_complex_and_set_up_the_filtration_to_distance_from_selected_cubes. The program will now terminate.\n";
				throw "Different embedding dimensions of cubes in the procedure construct_cubical_complex_and_set_up_the_filtration_to_distance_from_selected_cubes. The program will now terminate.\n";
			}
			for ( size_t j = 0 ; j != top_dimensional_cubes_that_are_in_the_set[0].size() ; ++j )
			{
				if ( sizes[j] < top_dimensional_cubes_that_are_in_the_set[i][j] )sizes[j] = top_dimensional_cubes_that_are_in_the_set[i][j];
			} 
		}
		//since there is a cell that assigns maximum, we need one more top dimensional cube in each direction:
		for ( size_t i = 0 ; i != sizes.size() ; ++i )
		{
			++sizes[i];
		}
		if ( dbg )
		{
			std::cerr << "Here are the sizes of the cubical complex to create (in max cubes units) \n";
			for ( size_t i = 0 ; i != sizes.size() ; ++i )
			{
				std::cout << i << " " << sizes[i] << std::endl;
			}
		}
	}
	//compute total number of maximal cells:
	size_t number_of_maximal_cells = 1;
	for ( size_t i = 0 ; i != sizes.size() ; ++i )number_of_maximal_cells *= sizes[i];
	
	if ( dbg )std::cout << "number_of_maximal_cells : " << number_of_maximal_cells << std::endl;
	
	//now the sizes are set up, we have the data to create the complex:	
	std::vector< typename Cubical_complex::filtration_type > top_dimensional_cells( number_of_maximal_cells , std::numeric_limits< typename Cubical_complex::filtration_type >::max() );
	
	//and create the cubical complex:	
	Cubical_complex* result = new Cubical_complex( sizes , top_dimensional_cells );
	
	std::vector< size_t > next_iteration_top_dimensional_cubes;
	
	//set the value 0 for all the cubes in top_dimensional_cubes_that_are_in_the_set:
	for ( size_t cube_no = 0 ; cube_no != top_dimensional_cubes_that_are_in_the_set.size() ; ++cube_no )
	{
		if ( dbg )
		{
			std::cout << "Setting up the value of the cube : " << std::endl;
			for ( size_t i = 0 ; i != top_dimensional_cubes_that_are_in_the_set[cube_no].size() ; ++i )
			{
				std::cout << top_dimensional_cubes_that_are_in_the_set[cube_no][i] << " ";
			}
		}
		size_t position = result->give_position_of_top_dimensional_cell( top_dimensional_cubes_that_are_in_the_set[cube_no] );
		if ( dbg )
		{
			std::cout << "Its position in the bitmap : " << position << std::endl;
		}
		result->get_cell_data(position) = 0;
		std::vector< size_t > neigs = result->give_neighbouring_top_dimensional_cells( position );
		next_iteration_top_dimensional_cubes.insert( next_iteration_top_dimensional_cubes.end() , neigs.begin() , neigs.end() );
		if ( dbg )
		{
			std::cout << "And here are its top dimensional neighs : " << std::endl;
			for ( size_t i = 0 ; i != neigs.size() ; ++i )
			{
				std::cout << neigs[i] << " ";
			}
			std::cout << std::endl;
			//getchar();
		}
	}
	
	std::cout << "Initial stuff set, we can move on \n";
	
	size_t counter = 1;
	//now for every cube in the neighboorhood of the cubes that have the value already assigned:
	while ( !next_iteration_top_dimensional_cubes.empty() )
	{
		if ( dbg )std::cout << "next_iteration_top_dimensional_cubes.size() : " << next_iteration_top_dimensional_cubes.size() << std::endl;
		std::vector< size_t > next_next_iteration_top_dimensional_cubes;
		for ( size_t i = 0 ; i != next_iteration_top_dimensional_cubes.size() ; ++i )
		{
			if ( result->get_cell_data( next_iteration_top_dimensional_cubes[i] ) != std::numeric_limits<typename Cubical_complex::filtration_type>::max() )continue;//this cube was already set.
			//if we are here, that means that the value of this cube was not set yet. We will now set it to counter.
			result->get_cell_data( next_iteration_top_dimensional_cubes[i] )  = counter;
			if ( dbg )
			{
				std::cerr << "The maximal cube at the position : " << next_iteration_top_dimensional_cubes[i] << " gets the value : " << counter << std::endl;
			}
			//and add the top dimensional neighs of this cube to the next_next_iteration_top_dimensional_cubes
			std::vector< size_t > neighs = result->give_neighbouring_top_dimensional_cells( next_iteration_top_dimensional_cubes[i] );
			next_next_iteration_top_dimensional_cubes.insert( next_next_iteration_top_dimensional_cubes.end() , neighs.begin() , neighs.end() );
		}
		next_iteration_top_dimensional_cubes.swap( next_next_iteration_top_dimensional_cubes );
		++counter;
	}
	
	if ( dbg )std::cout << "We are done, now imposing the lower star filtration \n";
	
	result->impose_lower_star_filtration();
	
	return result;
}//construct_cubical_complex_and_set_up_the_filtration_to_distance_from_selected_cubes
 
}//namespace Cubical_complex 
}//namespace Gudhi 

#endif
