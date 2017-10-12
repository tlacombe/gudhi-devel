/*    This file is part of the Gudhi Library. The Gudhi lhibrary
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Pawel Dlotko
 *
 *    Copyright (C) 2017, Swansea University UK
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
#include <tbb/parallel_for.h>
#endif

#include <limits>
#include <utility>
#include <algorithm>
#include <vector>
#include <numeric>

namespace Gudhi {

namespace Topological_inference_with_cubical_complexes {

/**
 * \brief Topological inference class.
 *
 * \ingroup Topological_inference_with_cubical_complexes
 *
 * \details
 * This is a topological inference class. Given a function, and a rectangular domain in R^n
 * it computes the values of the function on the top dimensional cubes of the domain. Later
 * it compute persistence of a function obtained in this way.
 * The file gudhi/functions_for_topological_inference/functions_for_topological_inference.h contains
 * a wide range of functions that can be used here. Those include functions define on point clouds.
 * As for the template parameters, T is a class of cubical complex that implements the following operations:
 * impose_lower_star_filtration();
 * initialize_simplex_associated_to_key();
 * set_the_value_of_top_dimensional_cell( std::vector<unsigned> , value );		
 * store_in_perseus_format( const char* ); 
 * For instance, the Bitmap_cubical_complex<whatever> satisy all the requirements. 
 * K is valye stored in the grid (typically double). 
 * F is a function that computes value on a grid points. It acts from K^n -> K. Please consult 
 * gudhi/functions_for_topological_inference/functions_for_topological_inference.h
 * for examples. 
**/
template <typename T , typename K , typename F>
class Topological_inference : public T
{
public:
	/**
	 * Default constructor. 
	**/ 
	Topological_inference():T(){};
	
	/**
	 * Constructor taking as parameters:
	 * std::vector< std::pair<K , K> >& coordinates_of_grid_ - a vector of pairs of points
	 * such that the pair of points in the position i of the vector is a projection
	 * of the cubical complex on the i-th direction. For example, for a cubical complex
	 * meshing the box [-1,1]x[-2,2] the vector contain two pairs: (-1,1) and (-2,2).
	 * The second parameter of the constructor is std::vector< unsigned >& resolution_of_a_grid_ 
	 * The vector need to have the same size as the first parameter of the function. It determines
	 * the number of maximal cubes in each direction. The last parmeteris a function that is to be computed
	 * on top dimensional cubes. 
	 * Please use this consructor together with T being non-periodic cubical complexes.  
	**/
	Topological_inference( const std::vector< std::pair<K , K> >& coordinates_of_grid_ , const std::vector< unsigned >& resolution_of_a_grid_ , F& f );
	
	
	/**
	 * Constructor taking as parameters:
	 * std::vector< std::pair<K , K> >& coordinates_of_grid_ - a vector of pairs of points
	 * such that the pair of points in the position i of the vector is a projection
	 * of the cubical complex on the i-th direction. For example, for a cubical complex
	 * meshing the box [-1,1]x[-2,2] the vector contain two pairs: (-1,1) and (-2,2).
	 * The second parameter of the constructor is std::vector< unsigned >& resolution_of_a_grid_ 
	 * The vector need to have the same size as the first parameter of the function. It determines
	 * the number of maximal cubes in each direction. The last but one parmeteris a function that 
	 * is to be computed on top dimensional cubes. The last parameter is determines the vector of
	 * directions in which periodic boundary conditions are to be imposed. 
	 * Please use this consructor together with T being a periodic cubical complexes.  
	**/ 
	Topological_inference( const std::vector< std::pair<K , K> >& coordinates_of_grid_ , const std::vector< unsigned >& resolution_of_a_grid_ , F& f , 
						   const std::vector<bool>& directions_in_which_periodic_b_cond_are_to_be_imposed );
	
	
	/**
	 * A procedure that store maximal cells to a file in Perseus format: 
	**/ 
	void write_to_file_Perseus_format( const char* filename );
	
	/**
	 * For debugging purposes, it write the filtration of top dimensional cell to a file in the way that
	 * the newlines are in the same places, where the cibical complex terminates in x direction. 
	**/ 	
	void write_to_file_with_newlines_at_the_ends_of_structure( const char* filename );
	
protected:	
    /**
     * Let us assume that a cubical complex store only top dimensional cubes. 
     * Index gives a lexicographical position of a cube in the complex. This procedure
     * compute the counter for that index, i.e. the coordinates of the cube.
    **/ 
    inline std::vector< unsigned > compute_counter_for_maximal_cube( size_t maximal_cube_index )
	{
		std::vector< unsigned > result( this->coordinates_of_grid.size() , 0 );
		for ( size_t i = 0 ; i != this->coordinates_of_grid.size() ; ++i )
		{
			result[i] = maximal_cube_index%this->resolution_of_a_grid[i];
			maximal_cube_index = maximal_cube_index/this->resolution_of_a_grid[i];				
		}
		return result;
	}

	/**
	 * Given a counter produced by compute_counter_for_maximal_cube, the following
	 * procedure find a center of this top dimensional cube. This information is later
	 * used to compute the valu of the function on that top dimensional cube.
	**/ 
	inline std::vector< K > compute_center_of_cube_for_given_counter( const std::vector< unsigned >& counter )
	{
		bool dbg = false;
		if ( dbg )
		{
			std::cerr << "Entering compute_center_of_cube_for_given_counter procedure \n";
			std::cout << "counter.size() : : " << counter.size() << std::endl;
			std::cout << "this->coordinates_of_grid.size() : " << this->coordinates_of_grid.size() << std::endl;
		}
		if ( counter.size() != this->coordinates_of_grid.size() )throw "Wrong dimensionality of a counter in the procedure compute_center_of_cube_for_given_counter \n";
		std::vector< K > result( counter.size() , 0 );
			
		for ( size_t dim = 0 ; dim != counter.size() ; ++dim )
		{
			if ( counter[dim] >= this->resolution_of_a_grid[dim] )throw "The counter in some dimension extends dimensionality of a grid. The program will now terminate \n";					
			result[dim] = dx_vector[dim]*counter[dim] + this->coordinates_of_grid[dim].first + dx_vector[dim]/2.0;
		}	
				
		return result;
	}
	
	/**
	 * Function that construct topological_inference object, used in all constructors. 
	**/ 
	void construct_topological_inference_object();
	
	//data:
	std::vector< std::pair< K , K > > coordinates_of_grid;
	std::vector< unsigned > resolution_of_a_grid;
	std::vector< K > dx_vector;
	F& f;
};



template <typename T , typename K , typename F>
Topological_inference<T,K,F>::Topological_inference( const std::vector< std::pair<K , K> >& coordinates_of_grid_ , const std::vector< unsigned >& resolution_of_a_grid_ , F& f )
												 :T(resolution_of_a_grid_), f(f)
{	
	this->coordinates_of_grid = coordinates_of_grid_;
	this->resolution_of_a_grid = resolution_of_a_grid_;	
	this->construct_topological_inference_object();
}//Topological_inference


template <typename T , typename K , typename F>
Topological_inference<T,K,F>::Topological_inference( const std::vector< std::pair<K , K> >& coordinates_of_grid_ , 
												     const std::vector< unsigned >& resolution_of_a_grid_ , F& f , 
													 const std::vector<bool>& directions_in_which_periodic_b_cond_are_to_be_imposed ):
													  T(resolution_of_a_grid_,directions_in_which_periodic_b_cond_are_to_be_imposed), 
													  f(f)
{
	this->coordinates_of_grid = coordinates_of_grid_;
	this->resolution_of_a_grid = resolution_of_a_grid_;	
	this->construct_topological_inference_object();
}//Topological_inference


template <typename T , typename K , typename F>
void Topological_inference<T,K,F>::construct_topological_inference_object( )
{
	bool dbg = false;
	if ( dbg )
	{
		std::cerr << "Entering constructor of a Topological_inference object \n";
		std::cout << "coordinates_of_grid_.size() : " << this->coordinates_of_grid.size() << std::endl;
		std::cout << "resolution_of_a_grid_.size() : " << this->resolution_of_a_grid.size() << std::endl;
	}
	if ( this->coordinates_of_grid.size() != this->resolution_of_a_grid.size() )throw "Incompatible sizes of coorfiantes of a grid, and the resoution of a grid in the constructore of Topological_inference. The program will now terminate \n";
		
	size_t number_of_maximal_cubes = 1;
	for ( size_t i = 0 ; i != this->resolution_of_a_grid.size() ; ++i ) number_of_maximal_cubes *= this->resolution_of_a_grid[i];
	
	this->dx_vector = std::vector< K >(this->coordinates_of_grid.size());
	for ( size_t dim = 0 ; dim != this->coordinates_of_grid.size() ; ++dim )
	{				
		this->dx_vector[dim] = (this->coordinates_of_grid[dim].second - this->coordinates_of_grid[dim].first)/this->resolution_of_a_grid[dim];
	}
	
	
#ifdef GUDHI_USE_TBB    
	tbb::parallel_for(size_t(0), number_of_maximal_cubes, [=](size_t i) 
#else 
     for ( size_t i = 0 ; i < number_of_maximal_cubes ; ++i )
#endif		
	{
		std::vector< unsigned > counter = this->compute_counter_for_maximal_cube( i );		
		std::vector< K > point = this->compute_center_of_cube_for_given_counter( counter );				
		K value = this->f( point );
		//values_on_maximal_cells[i] = value;		
		this->set_the_value_of_top_dimensional_cell( counter , value );			
	}	
#ifdef GUDHI_USE_TBB    
	);
#endif
	
	std::vector<size_t> counter_v( this->coordinates_of_grid.size() , 0 );	
	this->impose_lower_star_filtration();
	this->initialize_simplex_associated_to_key();	
}



template <typename T , typename K , typename F>
void Topological_inference<T,K,F>::write_to_file_Perseus_format( const char* filename )
{
	this->store_in_perseus_format( filename );
}

template <typename T , typename K , typename F>
void Topological_inference<T,K,F>::write_to_file_with_newlines_at_the_ends_of_structure( const char* filename )
{
	std::ofstream out(filename);
	unsigned counter = 0;
	for ( auto it = this->top_dimensional_cells_iterator_begin() ; it != this->top_dimensional_cells_iterator_end() ; ++it )
	{
		out << this->get_cell_data(*it) << " ";				
		if ( counter%this->resolution_of_a_grid[0] == this->resolution_of_a_grid[0]-1 )out << std::endl;
		++counter;	
	}
	out.close();
}//write_to_file_with_newlines_at_the_ends_of_structure

}  // namespace Topological_inference_with_cubical_complexes

}  // namespace Gudhi

#endif  //TOPOLOGICAL_INFERENCE_H_
	
