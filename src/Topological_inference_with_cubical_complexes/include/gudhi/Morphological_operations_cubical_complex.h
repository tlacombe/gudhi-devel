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

#ifndef MORPHOLOGICAL_OPERATIONS_CUBICAL_COMPLEX_H_
#define MORPHOLOGICAL_OPERATIONS_CUBICAL_COMPLEX_H_

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
enum considered_neighberhoods { all , full_face };

/**
 
**/ 

/**
 * This is a predicate that can be used to pinpoint the initial cubes for diltion.
 * The function check_predicaten in this class evaluate to true if the given value is greater or equal the given
 * cutoff value.
**/ 
template <typename T>
class Filtration_above_certain_value
{
public:
    Filtration_above_certain_value( T cutoff_value_ ):cutoff_value(cutoff_value_){}
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
class Filtration_below_certain_value
{
public:
    Filtration_below_certain_value( T cutoff_value_ ):cutoff_value(cutoff_value_){}
    bool check_predicate( double filtration_value ){return (filtration_value < this->cutoff_value);}
protected:
	T cutoff_value;
};

/**
 * This is a predicate that can be used to pinpoint the initial cubes for diltion.
 * The function check_predicaten in this class evaluate to true if the given value is between two given cutoff values.
**/
template <typename T>
class Filtration_in_range
{
public:
    Filtration_in_range( T cutoff_value_min_ , T cutoff_value_max_ ):cutoff_value_min(cutoff_value_min_),cutoff_value_max(cutoff_value_max_){}
    bool check_predicate( double filtration_value )
    {
		return ( (filtration_value >= this->cutoff_value_min) && (filtration_value <= this->cutoff_value_max) );
	}
protected:
	T cutoff_value_min;
	T cutoff_value_max;
};	

/**
 * Type of morphological operation to be performed: erosion, dilation or both.
**/ 
enum Operation_type { dylation_ , erosion_ , both_ };

/**
 * This is a predicate that can be used to pinpoint the initial cubes for diltion.
 * The function check_predicaten in this class evaluate to true if the given value is equal to the cutoff values.
**/
template <typename T>
class Filtration_equal
{
public:
    Filtration_equal( T cutoff_value_ ):cutoff_value(cutoff_value_){}
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
template <typename Cubical_complex ,typename Predicator>
class Morphological_operations_cubical_complex	
{
public:
//contructors	
	Morphological_operations_cubical_complex( Cubical_complex& to_process , Predicator& pred_ ):cmplx(to_process),pred(pred_)
	{
		this->apply_predicate_on_top_dimensional_cubes();
		this->cmplx.impose_lower_star_filtration();
		//Set up stuff for perssitence computations
	}
	Morphological_operations_cubical_complex( const std::vector< std::vector< unsigned > >& top_dimensional_cubes_that_are_in_the_set, std::vector< unsigned > sizes = std::vector< unsigned >() )
	{
		this->cmplx = this->construct_cubical_complex( top_dimensional_cubes_that_are_in_the_set , sizes );
		this->apply_predicate_on_top_dimensional_cubes();
		this->cmplx.impose_lower_star_filtration();
		//Set up stuff for perssitence computations		
	}

//methods	
	void erosion( typename Cubical_complex::filtration_type step_size , considered_neighberhoods neigh );
	void dilation( typename Cubical_complex::filtration_type step_size , considered_neighberhoods neigh );
	void both_erosion_and_dilation( typename Cubical_complex::filtration_type step_size , considered_neighberhoods neigh_type );
	void morphological_operation( typename Cubical_complex::filtration_type step_size , considered_neighberhoods neigh_type , Operation_type oper_type );		
	
	 /**
     * Assuming that the only values of maximal cubes are value_of_cubes_in_set and value_of_cubes_not_in_set,
     * it exchange the values
    **/ 
    void compute_complement();
	
	
private:
	//methods
	Cubical_complex* construct_cubical_complex
	( const std::vector< std::vector< unsigned > >& top_dimensional_cubes_that_are_in_the_set , 
    std::vector< unsigned > sizes = std::vector< unsigned >() );    
    
    void apply_predicate_on_top_dimensional_cubes();      
    
	//data structures:
	Cubical_complex& cmplx;
    Predicator& pred; 
    
    static typename Cubical_complex::filtration_type value_of_cubes_in_set;
    static typename Cubical_complex::filtration_type value_of_cubes_not_in_set;
};

template <typename Cubical_complex ,typename Predicator>
typename Cubical_complex::filtration_type Morphological_operations_cubical_complex<Cubical_complex,Predicator>::value_of_cubes_in_set = 0;

template <typename Cubical_complex ,typename Predicator>
typename Cubical_complex::filtration_type Morphological_operations_cubical_complex<Cubical_complex,Predicator>::value_of_cubes_not_in_set = 
std::numeric_limits< typename Cubical_complex::filtration_type >::infinity();


template <typename Cubical_complex , typename Predicator>
void Morphological_operations_cubical_complex<Cubical_complex,Predicator>::compute_complement()
{
	for ( auto it = this->cmplx.top_dimensional_cells_iterator_begin() ; it != this->cmplx.top_dimensional_cells_iterator_end() ; ++it )
	{
		typename Cubical_complex::filtration_type& cell_data = this->cmplx.get_cell_data( *it );
		if ( cell_data == this->value_of_cubes_in_set )
		{
			cell_data = this->value_of_cubes_not_in_set;
		}
		else
		{
			if ( cell_data == this->value_of_cubes_not_in_set )
			{
				cell_data = this->value_of_cubes_in_set;
			}
		}
	}
}//compute_complement


template <typename Cubical_complex , typename Predicator>
void Morphological_operations_cubical_complex<Cubical_complex,Predicator>::morphological_operation( typename Cubical_complex::filtration_type step_size , considered_neighberhoods neigh_type , Operation_type oper_type )
{	
	bool dbg = false;
	//first we need to find the cubes_to_consider. Those are all the cubes 
	//from the set (which therefore have the value set to value_of_cubes_in_set)
	//that are neighbors of cubes that are not in the set. 
	
	//we use this to make sure that we are not processing the same cube many times.
	std::vector< bool > was_this_maximal_cube_already_considered( this->cmplx.number_cells() , false );
	
	//this is just to mark the cells in the set in the 	was_this_maximal_cube_already_considered vector.
	size_t number_of_cubes_to_consider = 0;
	for ( auto it = this->cmplx.top_dimensional_cells_iterator_begin() ; it != this->cmplx.top_dimensional_cells_iterator_end() ; ++it )
	{
		if ( this->cmplx.get_cell_data( *it ) == this->value_of_cubes_in_set )
		{			
			std::vector< typename Cubical_complex::position_index_type > neighs = (neigh_type == full_face)? 
			this->cmplx.get_all_top_dimensional_cubes_sharing_codimension_1_face_with_given_top_dimensional_cube( *it ) :
			this->cmplx.get_all_top_dimensional_cubes_incident_to_the_given_top_dimensional_cell( *it );
			
			bool does_it_have_neigs_in_the_set_complement = false;
			for ( size_t neigh = 0 ; neigh != neighs.size() ; ++neigh )
			{		
				if ( this->cmplx.get_cell_data( neighs[neigh] ) == this->value_of_cubes_not_in_set )
				{
					does_it_have_neigs_in_the_set_complement = true;
					break;
				}
			}	
			if ( does_it_have_neigs_in_the_set_complement )
			{
				was_this_maximal_cube_already_considered[ *it ] = true;
				++number_of_cubes_to_consider;
			}
		}
	}
	
	
	std::vector< typename Cubical_complex::position_index_type > cubes_to_consider;
	cubes_to_consider.reserve( number_of_cubes_to_consider );
	for ( auto it = this->cmplx.top_dimensional_cells_iterator_begin() ; it != this->cmplx.top_dimensional_cells_iterator_end() ; ++it )
	{			
		if ( was_this_maximal_cube_already_considered[ *it ] )
		{
			//this is a cell from the set, now we will search for its neighbors
			//that are not in the set.
			std::vector< typename Cubical_complex::position_index_type > neighs = (neigh_type == full_face)? 
			this->cmplx.get_all_top_dimensional_cubes_sharing_codimension_1_face_with_given_top_dimensional_cube( *it ) :
			this->cmplx.get_all_top_dimensional_cubes_incident_to_the_given_top_dimensional_cell( *it );
			
			for ( size_t neigh = 0 ; neigh != neighs.size() ; ++neigh )
			{		
				if ( !was_this_maximal_cube_already_considered[ neighs[neigh] ] )
				{
					cubes_to_consider.push_back( neighs[neigh] );
				}
			}	
		}		
	}
	
	
	typename Cubical_complex::filtration_type value = step_size;
	while ( !cubes_to_consider.empty() )
	{	
		if ( dbg )
		{
			std::cout << "Here are the cubes to consider : \n";
			for ( size_t i = 0 ; i != cubes_to_consider.size() ; ++i )
			{
				std::cout << cubes_to_consider[i] << " ";
			}
			std::cout << std::endl;
			getchar();
		}
		
		//now we iterate through cubes_to_consider, and find the neighbouring cells:
		std::vector< typename Cubical_complex::position_index_type > new_cubes_to_consider;	
		new_cubes_to_consider.reserve( this->cmplx.number_of_top_dimnensional_cells( ) );				
		
		for ( size_t i = 0 ; i != cubes_to_consider.size() ; ++i )
		{	
			//if we have already considered it, there is no need to do it again	
			if ( was_this_maximal_cube_already_considered[ cubes_to_consider[i] ] )continue;
			//if not, mark it as consider:
			was_this_maximal_cube_already_considered[ cubes_to_consider[i] ] = true;
													
			typename Cubical_complex::filtration_type& filtration_of_cube = this->cmplx.get_cell_data( cubes_to_consider[i] );
			
			if (dbg)
			{
				std::cout << "filtration_of_cube : " << filtration_of_cube <<std::endl;
				std::cout << "this->value_of_cubes_in_set : " << this->value_of_cubes_in_set << std::endl;
				std::cout << "oper_type : " << oper_type << std::endl;
				getchar();
			}
			
			if ( filtration_of_cube == this->value_of_cubes_not_in_set )
			{
				//this cube is an element of complement. If we perform dylation or both dylation and erosion, 
				//this cube will get the following value:
				if ( (oper_type==dylation_) || (oper_type==both_) )
				{
					if (dbg)std::cout << "Dilatipn or both \n";
					filtration_of_cube = value;
				}
				else
				{
					//in this case, we do not care about this cube.
					continue;
				}
			}
			else
			{
				//this cube is an element of complement. If we perform erosion or both dylation and erosion, 
				//this cube will get the following value:
				if ( (oper_type==erosion_) || (oper_type==both_) )
				{
					if (dbg)std::cout << "Erosion or both \n";
					filtration_of_cube = -value;
				}
				else
				{
					//in this case, we do not care about this cube.
					continue;
				}
			}			
		
		
		
					
			std::vector< typename Cubical_complex::position_index_type > neighs = (neigh_type == full_face)? 
			this->cmplx.get_all_top_dimensional_cubes_sharing_codimension_1_face_with_given_top_dimensional_cube( cubes_to_consider[i] ) :
			this->cmplx.get_all_top_dimensional_cubes_incident_to_the_given_top_dimensional_cell( cubes_to_consider[i] );											
			
			for ( size_t neigh = 0 ; neigh != neighs.size() ; ++neigh )
			{				
				if ( !was_this_maximal_cube_already_considered[ neighs[neigh] ] )
				{
					new_cubes_to_consider.push_back( neighs[neigh] );
				}								
			}							
		
		
		
			
		}
				 
		cubes_to_consider = new_cubes_to_consider;
		value += step_size;
	}			
	
	//now we need to coean up the data strucutres after scrambling the filtration values.     
	this->cmplx.impose_lower_star_filtration();
	this->cmplx.initialize_arrays_for_persistence_computation();
}//morphological_operation


template <typename Cubical_complex , typename Predicator>
void Morphological_operations_cubical_complex<Cubical_complex,Predicator>::dilation( typename Cubical_complex::filtration_type step_size , considered_neighberhoods neigh_type )
{
	this->morphological_operation( step_size , neigh_type , dylation_ );
}//dilation


template <typename Cubical_complex , typename Predicator>
void Morphological_operations_cubical_complex<Cubical_complex,Predicator>::erosion( typename Cubical_complex::filtration_type step_size , considered_neighberhoods neigh_type )
{
	this->morphological_operation( step_size , neigh_type , erosion_ );
}//erosion


template <typename Cubical_complex , typename Predicator>
void Morphological_operations_cubical_complex<Cubical_complex,Predicator>::both_erosion_and_dilation( typename Cubical_complex::filtration_type step_size , considered_neighberhoods neigh_type )
{
	this->morphological_operation( step_size , neigh_type , both_ );
}//both_erosion_and_dilation

/*
template <typename Cubical_complex , typename Predicator>
void Morphological_operations_cubical_complex<Cubical_complex,Predicator>::dilation( typename Cubical_complex::filtration_type step_size , considered_neighberhoods neigh_type )
{	
	bool dbg = false;
	//first we need to find the cubes_to_consider:
	std::vector< typename Cubical_complex::position_index_type > cubes_to_consider;
	cubes_to_consider.reserve( this->cmplx.number_of_top_dimnensional_cells() );
	
	//we use this to make sure that we are not processing the same cube many times.
	std::vector< bool > was_this_maximal_cube_already_considered( this->cmplx.number_cells() , false );
	
	//this is just to makr the cells in the set in the 	was_this_maximal_cube_already_considered vector.
	for ( auto it = this->cmplx.top_dimensional_cells_iterator_begin() ; it != this->cmplx.top_dimensional_cells_iterator_end() ; ++it )
	{
		if ( this->cmplx.get_cell_data( *it ) == this->value_of_cubes_in_set )
		{			
			was_this_maximal_cube_already_considered[ *it ] = true;	
		}
	}
	
	for ( auto it = this->cmplx.top_dimensional_cells_iterator_begin() ; it != this->cmplx.top_dimensional_cells_iterator_end() ; ++it )
	{			
		if ( was_this_maximal_cube_already_considered[ *it ] )
		{
			//this is a cell from the set, now we will search for its neighbors
			//that are not in the set.
			std::vector< typename Cubical_complex::position_index_type > neighs;			
			if ( neigh_type == full_face )
			{			
				neighs = this->cmplx.get_all_top_dimensional_cubes_sharing_codimension_1_face_with_given_top_dimensional_cube( *it );
			}
			else
			{
				if ( neigh_type == all )
				{
					neighs = this->cmplx.get_all_top_dimensional_cubes_incident_to_the_given_top_dimensional_cell( *it );
				}
			}
			for ( size_t neigh = 0 ; neigh != neighs.size() ; ++neigh )
			{		
				if ( !was_this_maximal_cube_already_considered[neighs[neigh]] )
				{
					cubes_to_consider.push_back( neighs[neigh] );
				}
			}	
		}		
	}
	
	
	typename Cubical_complex::filtration_type value = step_size;
	while ( !cubes_to_consider.empty() )
	{	
		if ( dbg )
		{
			std::cout << "Here are the cubes to consider : \n";
			for ( size_t i = 0 ; i != cubes_to_consider.size() ; ++i )
			{
				std::cout << cubes_to_consider[i] << " ";
			}
			std::cout << std::endl;
			getchar();
		}
		
		//now we iterate through cubes_to_consider, and find the neighbouring cells:
		std::vector< typename Cubical_complex::position_index_type > new_cubes_to_consider;	
		new_cubes_to_consider.reserve( this->cmplx.number_of_top_dimnensional_cells( ) );				
		
		for ( size_t i = 0 ; i != cubes_to_consider.size() ; ++i )
		{	
			//if we have already considered it, there is no need to do it again	
			if ( was_this_maximal_cube_already_considered[ cubes_to_consider[i] ] )continue;
			//if not, mark it as consider:
			was_this_maximal_cube_already_considered[ cubes_to_consider[i] ] = true;
													
			typename Cubical_complex::filtration_type& filtration_of_cube = this->cmplx.get_cell_data( cubes_to_consider[i] );
			filtration_of_cube = value;
			
			
			std::vector< typename Cubical_complex::position_index_type > neighs;			
			if ( neigh_type == full_face )
			{			
				neighs = this->cmplx.get_all_top_dimensional_cubes_sharing_codimension_1_face_with_given_top_dimensional_cube( cubes_to_consider[i] );
			}
			else
			{
				if ( neigh_type == all )
				{
					neighs = this->cmplx.get_all_top_dimensional_cubes_incident_to_the_given_top_dimensional_cell( cubes_to_consider[i] );
				}
			}			
	
		
			
			for ( size_t neigh = 0 ; neigh != neighs.size() ; ++neigh )
			{				
				if ( was_this_maximal_cube_already_considered[ neighs[neigh] ] )
				{
					continue;
				}				
				new_cubes_to_consider.push_back( neighs[neigh] );
			}							
		}
				 
		cubes_to_consider = new_cubes_to_consider;
		value += step_size;
	}			
	//now we need to coean up the data strucutres after scrambling the filtration values.     
	this->cmplx.impose_lower_star_filtration();
	
	this->cmplx.initialize_arrays_for_persistence_computation();
}//dilation
*/


template <typename Cubical_complex , typename Predicator>
void Morphological_operations_cubical_complex<Cubical_complex,Predicator>::apply_predicate_on_top_dimensional_cubes()
{
	for ( auto it = this->cmplx.top_dimensional_cells_iterator_begin() ; it != this->cmplx.top_dimensional_cells_iterator_end() ; ++it )
	{
		auto& cell_data = this->cmplx.get_cell_data(*it);
		if ( pred.check_predicate( cell_data ) )
		{
			//std::cout << this->value_of_cubes_in_set << " ";			
			this->cmplx.get_cell_data(*it) = this->value_of_cubes_in_set;
		}
		else
		{
			//std::cout << this->value_of_cubes_not_in_set << " ";			
			this->cmplx.get_cell_data(*it) = this->value_of_cubes_not_in_set;
		}		
	}	
}



template <typename Cubical_complex , typename Predicator>
Cubical_complex* 
Morphological_operations_cubical_complex<Cubical_complex,Predicator>::
construct_cubical_complex( const std::vector< std::vector< unsigned > >& top_dimensional_cubes_that_are_in_the_set , std::vector< unsigned > sizes )
{	
	if ( sizes.size() == 0 ) 
	{
		//in this case, we need to read off the size of the bitmap from the top dimensional cells.
		//Here by default we create mesh minimal rectangle that contains all the top dimensional cubes. 
		sizes = std::vector<unsigned>( top_dimensional_cubes_that_are_in_the_set[0].size() , 0);
		
		for ( size_t top_dim_cube = 0 ; top_dim_cube != top_dimensional_cubes_that_are_in_the_set.size() ; ++top_dim_cube )
		{
			for ( size_t pos = 0 ; pos != top_dimensional_cubes_that_are_in_the_set[top_dim_cube].size() ; ++pos )
			{
				if ( top_dimensional_cubes_that_are_in_the_set[top_dim_cube][pos] > sizes[pos] )
				{
					sizes[pos] = top_dimensional_cubes_that_are_in_the_set[top_dim_cube][pos];
				}
			}
		}
	}
	
	//compute total number of maximal cells:
	size_t number_of_maximal_cells = 1;
	for ( size_t i = 0 ; i != sizes.size() ; ++i )number_of_maximal_cells *= sizes[i];		
	
	//now the sizes are set up, we have the data to create the complex:	
	std::vector< typename Cubical_complex::filtration_type > top_dimensional_cells( number_of_maximal_cells , std::numeric_limits< typename Cubical_complex::filtration_type >::infinity() );
	
	//and create the cubical complex:	
	Cubical_complex* result = new Cubical_complex( sizes , top_dimensional_cells );
			
	//set the value 0 for all the cubes in top_dimensional_cubes_that_are_in_the_set:
	//this can be TBB pararelized. 
	for ( size_t cube_no = 0 ; cube_no != top_dimensional_cubes_that_are_in_the_set.size() ; ++cube_no )
	{	
		size_t position = result->give_position_of_top_dimensional_cell( top_dimensional_cubes_that_are_in_the_set[cube_no] );
		result->get_cell_data(position) = 0;
	}	
	//over here we do not yet impose lower star filtration, since this complex will be processed further, and this will be done later. 		

	return result;
}//construct_cubical_complex_and_set_up_the_filtration_to_distance_from_selected_cubes
	


 
}//namespace TOPOLOGICAL_INFERENCE_H_ 
}//namespace Gudhi 

#endif
