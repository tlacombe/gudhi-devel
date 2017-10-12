/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Pawel Dlotko
 *
 *    Copyright (C) 2017 Swansea University UK
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



wydaje mi sie ze morphological operation byloby bardziej stosowna nazwa. 
#ifndef DISTANCE_FROM_COLLECTION_OF_CUBES_
#define DISTANCE_FROM_COLLECTION_OF_CUBES_

#include <gudhi/Bitmap_cubical_complex_base.h>
#include <gudhi/Bitmap_cubical_complex_periodic_boundary_conditions_base.h>
#include <gudhi/Bitmap_cubical_complex.h>
#ifdef GUDHI_USE_TBB
#include <tbb/parallel_sort.h>
#endif

#include <limits>
#include <utility>  
#include <algorithm>
#include <vector>
#include <cmath>
#include <numeric>  



namespace Gudhi 
{

namespace Topological_inference_with_cubical_complexes 
{
	

/**
 * When computing dylation of a set of cubes, we can either consider the direct neighbors of a cube C being all the cubes 
 * that have nonempty intersection C, or only the cubes that share co-dimension-1 face with 1. In the first case, all neighbourhood
 * should be choosen. In the former: full_face neighborhood. 
**/ 
enum considered_neighberhoods { all , full_face }

/**
 
**/ 

/**
 * This is a predicate that can be used to pinpoint the initial cubes for diltion.
 * The function check_predicaten in this class evaluate to true if the given value is greater or equal the given
 * cutoff value.
**/ 
template <typename T>
class filtration_above_certain_value
{
public:
    filtration_above_certain_value( T cutoff_value_ ):cutoff_value(cutoff_value_){}
    bool check_predicate( double filtration_value ){return (filtration_value >= this->cutoff_value);}
protected:
	T cutoff_value;
};

/**
 * This is a predicate that can be used to pinpoint the initial cubes for diltion.
 * The function check_predicaten in this class evaluate to true if the given value is smaller the given
 * cutoff value.
**/
template <typename T>
class filtration_below_certain_value
{
public:
    filtration_above_certain_value( T cutoff_value_ ):cutoff_value(cutoff_value_){}
    bool check_predicate( double filtration_value ){return (filtration_value < this->cutoff_value);}
protected:
	T cutoff_value;
};

/**
 * This is a predicate that can be used to pinpoint the initial cubes for diltion.
 * The function check_predicaten in this class evaluate to true if the given value is between two given cutoff values.
**/
template <typename T>
class filtration_in_range
{
public:
    filtration_above_certain_value( T cutoff_value_min_ , T cutoff_value_max_ ):cutoff_value_min(cutoff_value_min_),cutoff_value_max(cutoff_value_max_){}
    bool check_predicate( double filtration_value )
    {
		return ( (filtration_value >= this->cutoff_value_min) && (filtration_value <= this->cutoff_value_max) );
	}
protected:
	T cutoff_value_min;
	T cutoff_value_max;
};	

/**
 * This is a predicate that can be used to pinpoint the initial cubes for diltion.
 * The function check_predicaten in this class evaluate to true if the given value is equal to the cutoff values.
**/
template <typename T>
class filtration_equal
{
public:
    filtration_above_certain_value( T cutoff_value_ ):cutoff_value(cutoff_value_){}
    bool check_predicate( double filtration_value ){return (filtration_value == this->cutoff_value);}
protected:
	T cutoff_value;
};	
	
//**********************************************************************************************************************
//**********************************************************************************************************************
//**********************************************************************************************************************
//**********************************************************************************************************************
//**********************************************************************************************************************
//**********************************************************************************************************************



/**
 * The class Dylate_collection_of_cubes_in_cubical_complex construct an object of a type
 * Cubical_complex.
 **/ 
template <typenale Cubical_complex , Predicator>
class Morphological_operations_cubical_complex	
{
public:
//contructors	
	Morphological_operations_cubical_complex( Cubical_complex* to_process , const Predicator& pred_ ):cmplx(to_process),pred(pred_)
	{
		this->apply_predicate_on_top_dimensional_cubes();
		this->impose_lower_star_filtration();
	}
	Morphological_operations_cubical_complex( const std::vector< std::vector< unsigned > >& top_dimensional_cubes_that_are_in_the_set, std::vector< unsigned > sizes = std::vector< unsigned >() )
	{
		this->cmplx = this->construct_cubical_complex_and_set_up_the_filtration_to_distance_from_selected_cubes( top_dimensional_cubes_that_are_in_the_set , sizes );
		this->apply_predicate_on_top_dimensional_cubes();
		this->impose_lower_star_filtration();		
	}

//methods	
	void erosion( Cubical_complex::filtration_type step_size , considered_neighberhoods neigh );
	void dilation( Cubical_complex::filtration_type step_size , considered_neighberhoods neigh );
	void both_erosion_and_dilation( Cubical_complex::filtration_type step_size , considered_neighberhoods neigh );
	
	Cubical_complex* return_the_complex(){return this->cmplx;}
	
private:
	//methods
	Cubical_complex* construct_cubical_complex_and_set_up_the_filtration_to_distance_from_selected_cubes
	( const std::vector< std::vector< unsigned > >& top_dimensional_cubes_that_are_in_the_set , 
    std::vector< unsigned > sizes = std::vector< unsigned >() );    
    
    void apply_predicate_on_top_dimensional_cubes();
    
    
	//data structures:
	Cubical_complex* cmplx;
    Predicator& pred;
    considered_operation operation;    
};

template <typenale Cubical_complex , Predicator>
void Morphological_operations_cubical_complex::erosion( Cubical_complex::filtration_type step_size , considered_neighberhoods neigh )
{
	
	
}//erosion

template <typenale Cubical_complex , Predicator>
void Morphological_operations_cubical_complex::dilation( Cubical_complex::filtration_type step_size , considered_neighberhoods neigh )
{
	bool dbg = true;
	//first we need to find the cubes_to_consider:
	std::vector< Cubical_complex::position_index_type > cubes_to_consider;
	cubes_to_consider.reserve( this->cmplx->number_of_top_dimnensional_cells() );
	
	//we use this to make sure that we are not processing the same cube many times.
	std::vector< bool > was_this_maximal_cube_already_considered( this->cmplx->number_of_top_dimnensional_cells() , false );
	
	for ( auto it = this->cmplx->top_dimensional_cells_iterator_begin() ; it != this->top_dimensional_cells_iterator_end() ; ++it )
	{
		if ( this->cmplx->get_cell_data(*it) == std::numeric_limits< typename Cubical_complex::filtration_type >::infinity() )
		{
			cubes_to_consider.push_back( *it );
		}
		else
		{
			//in this case, this cube will never be consider 
			was_this_maximal_cube_already_considered[ *it ] = true;
		}
	}
	
	
	typename Cubical_complex::filtration_type value = step_size;
	while ( !cubes_to_consider.empty() )
	{	
		//now we iterate through cubes_to_consider, and find the neighbouring cells:
		std::vector< Cubical_complex::position_index_type > new_cubes_to_consider;			
		new_cubes_to_consider.reserve( cmplx.number_of_maximal_cubes() );
		
		for ( size_t i = 0 ; i != cubes_to_consider.size() ; ++i )
		{		
			if ( dbg )std::cerr << "Looking at the neighs of a cube number : " << cubes_to_consider[i] << std::endl;	
			
			std::vector< Cubical_complex::position_index_type > neighs;
			
			if ( neigh == full_face )
			{
				neighs = this->cmplx->get_all_top_dimensional_cubes_sharing_codimension_1_face_with_given_top_dimensional_cube( cubes_to_consider[i] );
			}
			else
			{
				if ( neigh == all )
				{
					neighs = this->get_all_top_dimensional_cubes_incident_to_the_given_top_dimensional_cell( cubes_to_consider[i] );
				}
			}			
			
			
			for ( size_t neigh = 0 ; neigh != neighs.size() ; ++neigh )
			{
				if ( was_this_maximal_cube_already_considered[ neighs[neigh] ] )continue;
				was_this_maximal_cube_already_considered[ neighs[neigh] ] = true;
				
				typename Cubical_complex::filtration_type filtration_of_cube& = this->cmplx->get_cell_data( neighs[neigh] );
				if ( filtration_of_cube == std::numeric_limits< typename Cubical_complex::filtration_type >::max() )
				{
					filtration_of_cube = value;												
					new_cubes_to_consider.push_back( neighs[neigh] );								
				}
			}							
		}
				 
		cubes_to_consider = new_cubes_to_consider;
		value += step_size;
	}			
	//now we need to coean up the data strucutres after scrambling the filtration values.     
	cmplx.impose_lower_star_filtration();
	
	
	

	//this is not generic enough!!! TODO, call it prepare for persistence computations or whatever!!!
	for (size_t i = 0; i != cmplx.total_number_of_cells; ++i) 
	{
		cmplx.key_associated_to_simplex[i] = i;
	}
	cmplx.initialize_simplex_associated_to_key();
}//dilation

template <typenale Cubical_complex , Predicator>
void Morphological_operations_cubical_complex::both_erosion_and_dilation( Cubical_complex::filtration_type step_size , considered_neighberhoods neigh )
{
	
}//both_erosion_and_dilation

template <typenale Cubical_complex , Predicator>
void Morphological_operations_cubical_complex::apply_predicate_on_top_dimensional_cubes()
{
	for ( auto it = this->cmplx->top_dimensional_cells_iterator_begin() ; it != this->top_dimensional_cells_iterator_end() ; ++it )
	{
		auto& cell_data = this->cmplx->get_cell_data(*it);
		if ( pred( cell_data ) )
		{
			this->cmplx->get_cell_data(*it) = 0;
		}
		else
		{
			this->cmplx->get_cell_data(*it) = std::numeric_limits< typename Cubical_complex::filtration_type >::infinity(); 
		}		
	}	
}


template <typenale Cubical_complex , Predicator>
Cubical_complex* Morphological_operations_cubical_complex::construct_cubical_complex_and_set_up_the_filtration_to_distance_from_selected_cubes( const std::vector< std::vector< unsigned > >& top_dimensional_cubes_that_are_in_the_set , 
std::vector< unsigned > sizes = std::vector< unsigned >() )
{	
	//compute total number of maximal cells:
	size_t number_of_maximal_cells = 1;
	for ( size_t i = 0 ; i != sizes.size() ; ++i )number_of_maximal_cells *= sizes[i];		
	
	//now the sizes are set up, we have the data to create the complex:	
	std::vector< typename Cubical_complex::filtration_type > top_dimensional_cells( number_of_maximal_cells , std::numeric_limits< typename Cubical_complex::filtration_type >::max() );
	
	//and create the cubical complex:	
	Cubical_complex* result = new Cubical_complex( sizes , top_dimensional_cells );
	
	std::vector< size_t > next_iteration_top_dimensional_cubes;
	
	//set the value 0 for all the cubes in top_dimensional_cubes_that_are_in_the_set:
	//this can be TBB pararelized. 
	for ( size_t cube_no = 0 ; cube_no != top_dimensional_cubes_that_are_in_the_set.size() ; ++cube_no )
	{	
		size_t position = result->give_position_of_top_dimensional_cell( top_dimensional_cubes_that_are_in_the_set[cube_no] );
		result->get_cell_data(position) = 0;
	}	

	return result;
}//construct_cubical_complex_and_set_up_the_filtration_to_distance_from_selected_cubes
	


 
}//namespace TOPOLOGICAL_INFERENCE_H_ 
}//namespace Gudhi 

#endif
