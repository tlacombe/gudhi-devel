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
    Filtration_above_certain_value():cutoff_value(0){}
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
    Filtration_below_certain_value():cutoff_value(0){}
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
    Filtration_in_range():cutoff_value_min(0),cutoff_value_max(0){}
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
 * This is a predicate that can be used to pinpoint the initial cubes for diltion.
 * The function check_predicaten in this class evaluate to true if the given value is equal to the cutoff values.
**/
template <typename T>
class Filtration_equal
{
public:
    Filtration_equal():cutoff_value(0){}
    Filtration_equal( T cutoff_value_ ):cutoff_value(cutoff_value_){}
    bool check_predicate( double filtration_value ){return (filtration_value == this->cutoff_value);}
protected:
	T cutoff_value;
};	


/**
 * This is a preicate that always return true. It will be used aw a default one 
 * for the class Morphological_operations_cubical_complex.
**/
class Always_true
{
public:
    Always_true(){}    
    bool check_predicate( double filtration_value ){return true;}
};
	
//**********************************************************************************************************************
//**********************************************************************************************************************
//**********************************************************************************************************************
//**********************************************************************************************************************
//**********************************************************************************************************************
//**********************************************************************************************************************

/**
 * Type of morphological operation to be performed: erosion, dilation or both.
**/ 
enum Operation_type { dylation_ , erosion_ , both_ };

/**
 * The class Dylate_collection_of_cubes_in_cubical_complex. It either construct or proces 
 * an object of a type Cubical_complex by using Predicator function. By doing so, the whole
 * complex is treshold into the cubes that belongs to a 'set', and those that do now (in which
 * case we will say that they are in the 'complement'). The tresholding process happens when calling 
 * the constructor of Morphological_operations_cubical_complex. In case the cubical complex is created outside,
 * it will be modified by that class Morphological_operations_cubical_complex.
 * Omce an object of Morphological_operations_cubical_complex class is constructd, one can call three possible
 * operations on it: erosion, dilation or both of them. 
 * \ingroup Topological_inference_with_cubical_complexes
 **/ 
template <typename Cubical_complex ,typename Predicator = Always_true>
class Morphological_operations_cubical_complex	
{
public:
	/**
	 * This constructor take as a parameter the cubical complex to_process and a predicator that
	 * will be used to modify the cubical complex. All the maximal cubes the filtration of which
	 * satisfy the predicate will be assigned to the set. All the cubes that do not satisfy the predicate
	 * will be set to the set's complement.
	**/ 
	Morphological_operations_cubical_complex( Cubical_complex* to_process , Predicator& pred_ ):cmplx(to_process),pred(pred_)
	{
		this->apply_predicate_on_top_dimensional_cubes();
		this->cmplx->impose_lower_star_filtration();
		this->was_cubical_complex_created_by_this_class = false;
	}
	
	/**
	 * This constructor take as a parameter a collection of indices of top dimensonal cubes that are in the set, 
	 * and the vector<unsigned> that determine the size of the complex. It creates the cubica complex based on that. 
	 * The collection of top dimensional cubes is given by std::vector< std::vector< unsigned > > where the internal
	 * std::vector< unsigned > is a counter that deterine the position of each top dimensional cube. TODO
	**/ 
	Morphological_operations_cubical_complex( const std::vector< std::vector< unsigned > >& top_dimensional_cubes_that_are_in_the_set, 
	std::vector< unsigned > sizes = std::vector< unsigned >(),
	std::vector< bool > directions_in_which_to_impose_periodic_b_cond = std::vector< bool >()
	)
	{	this->cmplx = this->construct_cubical_complex( top_dimensional_cubes_that_are_in_the_set , sizes , directions_in_which_to_impose_periodic_b_cond );		
		this->cmplx->impose_lower_star_filtration();
		this->was_cubical_complex_created_by_this_class = true;
		//Set up stuff for perssitence computations		
	}
	
	/**
	 * Destructor.
	**/ 
	~Morphological_operations_cubical_complex()
	{
		if ( this->was_cubical_complex_created_by_this_class )
		{
			delete this->cmplx;
			//to make sure that this object is not disposed twice. 
			this->was_cubical_complex_created_by_this_class = false;
		}
	}
	
	/**
	 * Return this complex. Can be used for persistent homology computations.
	**/ 
	Cubical_complex* give_me_the_complex()
	{
		return this->cmplx;
	}
	

//methods	
    /**
     * This is a implementation of erosion operation. When it is performed, the maximal cubes that belong
     * to a set, will obtain new filtration values. Value 1 will be given to all the cubes which,
     * in the considered_neighberhoods, have a neighbor in the set's complement. Value 2 will be given to 
     * the cubes which are neighbors of cubes with filtration value 1 (assigned in the previous step). And
     * so on. 
    **/ 
	void erosion( typename Cubical_complex::filtration_type step_size , considered_neighberhoods neigh );

    /**
     * This is a implementation of dilation operation. When it is performed, the maximal cubes that belong
     * to a set's complement, will obtain new filtration values. Value 1 will be given to all the cubes 
     * from the set's complement which, in the considered_neighberhoods, have a neighbor in the set. Value 2 
     * will be given to the cubes from the set's complement which are neighbors of cubes with filtration value 1 
     * (assigned in the previous step). And so on.
    **/ 
	void dilation( typename Cubical_complex::filtration_type step_size , considered_neighberhoods neigh );
	
	/**
     * This procedure perform both erosion and dilation. For a desciption of erosion and dilation, please refere to the 
     * suitable methods. 
     * The procedure start from the top dimensional cubes in the set, that are neighbors of cubes in the set's complement.
     * Startign from them, bothe erosion and dilation is performed. 
     * If you would rather start from the sollection of cubes in the set complement that have neighbors in the set, please
     * use compute_complement() function to swap set and set complement. 
    **/ 
	void both_erosion_and_dilation( typename Cubical_complex::filtration_type step_size , considered_neighberhoods neigh_type );	
	
	 /**
     * Assuming that the only values of maximal cubes are value_of_cubes_in_set and value_of_cubes_not_in_set,
     * it exchange the values
    **/ 
    void compute_complement();
	
	
protected:
	//methods
	/**
	 * A general method that implements erosion, dilation, and both of them. 
	**/ 
	void morphological_operation( typename Cubical_complex::filtration_type step_size , considered_neighberhoods neigh_type , Operation_type oper_type );		
	
	/**
	 * A proedure that create a cubical complex. It is called by one of the consructors of the 
	 * Morphological_operations_cubical_complex class.
	**/ 
	Cubical_complex* construct_cubical_complex
	( const std::vector< std::vector< unsigned > >& top_dimensional_cubes_that_are_in_the_set , 
    std::vector< unsigned > sizes = std::vector< unsigned >() , std::vector< bool > directions_in_which_to_impose_periodic_b_cond = std::vector< bool >() );    
    
    /**
     * A function that iterate through all top dimensional cells and apply predicate on them. By doing so, it divide
     * all top dimensional cell to those which belong to the set (in which case they will get the new filtration value
     * value_of_cubes_in_set and of the set complement) and to the set complement (in which case they will get the 
     * filtration value value_of_cubes_not_in_set).
    **/ 
    void apply_predicate_on_top_dimensional_cubes();      
    
	//data structures:
	Cubical_complex* cmplx;
    Predicator pred;     
    bool was_cubical_complex_created_by_this_class;
    
    //Default values of elements in the set, and in its coplement. 
    static typename Cubical_complex::filtration_type value_of_cubes_in_set;
    static typename Cubical_complex::filtration_type value_of_cubes_not_in_set;
};

//Setting up the default values of cells in the set.
template <typename Cubical_complex ,typename Predicator>
typename Cubical_complex::filtration_type Morphological_operations_cubical_complex<Cubical_complex,Predicator>::value_of_cubes_in_set = 0;

//Setting up the default values of cells in the set's complement.
template <typename Cubical_complex ,typename Predicator>
typename Cubical_complex::filtration_type Morphological_operations_cubical_complex<Cubical_complex,Predicator>::value_of_cubes_not_in_set = 
std::numeric_limits< typename Cubical_complex::filtration_type >::infinity();


template <typename Cubical_complex , typename Predicator>
void Morphological_operations_cubical_complex<Cubical_complex,Predicator>::compute_complement()
{
	for ( auto it = this->cmplx->top_dimensional_cells_iterator_begin() ; it != this->cmplx->top_dimensional_cells_iterator_end() ; ++it )
	{
		typename Cubical_complex::filtration_type& cell_data = this->cmplx->get_cell_data( *it );
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
	std::vector< bool > was_this_maximal_cube_already_considered( this->cmplx->number_cells() , false );
	
	//this is just to mark the cells in the set in the 	was_this_maximal_cube_already_considered vector.
	size_t number_of_cubes_to_consider = 0;
	for ( auto it = this->cmplx->top_dimensional_cells_iterator_begin() ; it != this->cmplx->top_dimensional_cells_iterator_end() ; ++it )
	{
		if ( this->cmplx->get_cell_data( *it ) == this->value_of_cubes_in_set )
		{			
			std::vector< typename Cubical_complex::position_index_type > neighs = (neigh_type == full_face)? 
			this->cmplx->get_all_top_dimensional_cubes_sharing_codimension_1_face_with_given_top_dimensional_cube( *it ) :
			this->cmplx->get_all_top_dimensional_cubes_incident_to_the_given_top_dimensional_cell( *it );
			
			bool does_it_have_neigs_in_the_set_complement = false;
			for ( size_t neigh = 0 ; neigh != neighs.size() ; ++neigh )
			{		
				if ( this->cmplx->get_cell_data( neighs[neigh] ) == this->value_of_cubes_not_in_set )
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
	for ( auto it = this->cmplx->top_dimensional_cells_iterator_begin() ; it != this->cmplx->top_dimensional_cells_iterator_end() ; ++it )
	{			
		if ( was_this_maximal_cube_already_considered[ *it ] )
		{
			//this is a cell from the set, now we will search for its neighbors
			//that are not in the set.
			std::vector< typename Cubical_complex::position_index_type > neighs = (neigh_type == full_face)? 
			this->cmplx->get_all_top_dimensional_cubes_sharing_codimension_1_face_with_given_top_dimensional_cube( *it ) :
			this->cmplx->get_all_top_dimensional_cubes_incident_to_the_given_top_dimensional_cell( *it );
			
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
		new_cubes_to_consider.reserve( this->cmplx->number_of_top_dimnensional_cells( ) );				
		
		for ( size_t i = 0 ; i != cubes_to_consider.size() ; ++i )
		{	
			//if we have already considered it, there is no need to do it again	
			if ( was_this_maximal_cube_already_considered[ cubes_to_consider[i] ] )continue;
			//if not, mark it as consider:
			was_this_maximal_cube_already_considered[ cubes_to_consider[i] ] = true;
													
			typename Cubical_complex::filtration_type& filtration_of_cube = this->cmplx->get_cell_data( cubes_to_consider[i] );
			
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
			this->cmplx->get_all_top_dimensional_cubes_sharing_codimension_1_face_with_given_top_dimensional_cube( cubes_to_consider[i] ) :
			this->cmplx->get_all_top_dimensional_cubes_incident_to_the_given_top_dimensional_cell( cubes_to_consider[i] );											
			
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
	this->cmplx->impose_lower_star_filtration();
	this->cmplx->initialize_arrays_for_persistence_computation();
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


template <typename Cubical_complex , typename Predicator>
void Morphological_operations_cubical_complex<Cubical_complex,Predicator>::apply_predicate_on_top_dimensional_cubes()
{
	for ( auto it = this->cmplx->top_dimensional_cells_iterator_begin() ; it != this->cmplx->top_dimensional_cells_iterator_end() ; ++it )
	{
		auto& cell_data = this->cmplx->get_cell_data(*it);
		if ( pred.check_predicate( cell_data ) )
		{
			//std::cout << this->value_of_cubes_in_set << " ";			
			this->cmplx->get_cell_data(*it) = this->value_of_cubes_in_set;
		}
		else
		{
			//std::cout << this->value_of_cubes_not_in_set << " ";			
			this->cmplx->get_cell_data(*it) = this->value_of_cubes_not_in_set;
		}		
	}	
}



template <typename Cubical_complex , typename Predicator>
Cubical_complex*
Morphological_operations_cubical_complex<Cubical_complex,Predicator>::
construct_cubical_complex( const std::vector< std::vector< unsigned > >& top_dimensional_cubes_that_are_in_the_set , std::vector< unsigned > sizes , 
					       std::vector< bool > directions_in_which_to_impose_periodic_b_cond )
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
	if ( directions_in_which_to_impose_periodic_b_cond.size() == 0 )
	{
		directions_in_which_to_impose_periodic_b_cond = std::vector< bool >( sizes.size() , false );
	}
	
	//compute total number of maximal cells:
	size_t number_of_maximal_cells = 1;
	for ( size_t i = 0 ; i != sizes.size() ; ++i )
	{
		number_of_maximal_cells *= sizes[i];		
	}
	
	//now the sizes are set up, we have the data to create the complex:	
	std::vector< typename Cubical_complex::filtration_type > top_dimensional_cells( number_of_maximal_cells , this->value_of_cubes_not_in_set );
	//we create the complex:
	Cubical_complex* result = new Cubical_complex( sizes , top_dimensional_cells , directions_in_which_to_impose_periodic_b_cond );		
	
	//now for every top dimensional cube, given as a coutner, compute its index in the 
	//complex result:
	for ( size_t i = 0 ; i != top_dimensional_cubes_that_are_in_the_set.size() ; ++i )
	{
		std::vector<unsigned> bitmap_counter = result->convert_top_dimensional_cube_counter_to_bitmap_counter( top_dimensional_cubes_that_are_in_the_set[i] );
		typename Cubical_complex::position_index_type pos = result->compute_position_in_bitmap(bitmap_counter);
		result->get_cell_data(pos) = this->value_of_cubes_in_set;				
	}	
		
	//Over here we do not call the methods 	cmplx->impose_lower_star_filtration() or
	//cmplx->initialize_arrays_for_persistence_computation(), since this will be done in the
	//constructor. Since this method is private, the user will not be able to invoke it
	//and forget to call those two methods. If you are, in any case, changing this procedure,
	//keep in mind that the two methods lited above have to be invoked in order to have a 
	//consistent data structure. 
		
	return result;
}//construct_cubical_complex_and_set_up_the_filtration_to_distance_from_selected_cubes
	


 
}//namespace TOPOLOGICAL_INFERENCE_H_ 
}//namespace Gudhi 

#endif
