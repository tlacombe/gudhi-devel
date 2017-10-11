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


TODO:
3) Rozwazyc zarowno erozje jak dylacje -- przemyslec jak to zrobic algoritymicznie.


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
	






Jezeli rozwazymy zarowno erozje jak i dylacje, to trzeba zmienic nazwe klasy.
Btw, czy mozna te dwie rzeczy robic jednoczenie??
Jezeli tak, to przedzialy ktore maja dodatnie konce byyby z dylacji, a te ujemne z erozji. 
Byc moze mozna znowu emu azdefiniwac:
enum considered_operation { erosion, dilation, both };???
Trzeba wtedy w dokumentacji zaznaczyc, ze erosion dilation nie zawsze to to samo co dilation erosion.
/**
 * The class Dylate_collection_of_cubes_in_cubical_complex construct an object of a type
 * Cubical_complex.
 **/ 
template <typenale Cubical_complex , Predicator>
class Dylate_collection_of_cubes_in_cubical_complex	
{
public:
//how to relate constructors of this class to methods we have at the moment??
	Dylate_collection_of_cubes_in_cubical_complex(){}
	Dylate_collection_of_cubes_in_cubical_complex( Cubical_complex* to_process , const Predicator& pred );
	
private:
	void dilation( Cubical_complex* result , std::vector< size_t >& next_iteration_top_dimensional_cubes , Cubical_complex::filtration_type diameter_of_cube = 1)
	{
		
		typename Cubical_complex::filtration_type value = diameter_of_cube;
		while ( !cubes_to_consider.empty() )
		{	
			//now we iterate through cubes_to_consider, and find the neighbouring cells:
			std::vector< size_t > new_cubes_to_consider;			
			new_cubes_to_consider.reserve( cmplx.number_of_maximal_cubes() );
			
			for ( size_t i = 0 ; i != cubes_to_consider.size() ; ++i )
			{		
				if ( dbg )std::cerr << "Looking at the neighs of a cube number : " << cubes_to_consider[i] << std::endl;	
				
				tutaj w zaleznosci od sasiedztwa ktore wybieramy trzeba wolac rozne funkcje!!!
				std::vector< size_t > neighs = cmplx.give_neighbouring_top_dimensional_cells( cubes_to_consider[i] );
				for ( size_t neigh = 0 ; neigh != neighs.size() ; ++neigh )
				{
					if ( cmplx.data[ neighs[neigh] ] == std::numeric_limits< typename Cubical_complex::filtration_type >::max() )
					{
						cmplx.data[ neighs[neigh] ] = value;												
						new_cubes_to_consider.push_back( neighs[neigh] );								
					}
				}			
			}
			
			 
			cubes_to_consider = new_cubes_to_consider;
			value += diameter_of_cube;
		}	
			
		//now we need to coean up the data strucutres after scrambling the filtration values.     
		cmplx.impose_lower_star_filtration();
		for (size_t i = 0; i != cmplx.total_number_of_cells; ++i) 
		{
		    cmplx.key_associated_to_simplex[i] = i;
		}
		cmplx.initialize_simplex_associated_to_key();
		
		if ( dbg )std::cout << "We are done, now imposing the lower star filtration \n";
		
		
		upewnic sie, ze ta procedura dziala niezaleznie od tego jakie sa poczatkowe wartosci
		Tj. dodac odpowiednia flage
		result->impose_lower_star_filtration();
	}
	
	
	
	Widac, ze tutaj mamy dwie mozliwosci dzialania
	1) Mamy juz kompleks + predykat + sasiedztwo i chemy zrobic erozje / dylacje. 
	2) Nie mamy kompleksu, ale mamy liste kostek maksymalnych. Moze w tej sytuacji po prostu zarzadac od uzytkownika, by stworzyl nam taki kompleks
	i my zaplimy tutaj kostki maksymalne i zrobmy dylacje? A moze po prostu lepiej jest tak jak teraz, ze przyjmujemy parametry,a kompls stworzymy my tutaj
	lokalnie. Teraz jestem za 2g aopcja, tak jak jest. 
	Cubical_complex* construct_cubical_complex_and_set_up_the_filtration_to_distance_from_selected_cubes( const std::vector< std::vector< unsigned > >& top_dimensional_cubes_that_are_in_the_set , 
																										  std::vector< unsigned > sizes = std::vector< unsigned >() )
	{
		bool dbg = false;

		//compute total number of maximal cells:
		size_t number_of_maximal_cells = 1;
		for ( size_t i = 0 ; i != sizes.size() ; ++i )number_of_maximal_cells *= sizes[i];		
		
		//now the sizes are set up, we have the data to create the complex:	
		std::vector< typename Cubical_complex::filtration_type > top_dimensional_cells( number_of_maximal_cells , std::numeric_limits< typename Cubical_complex::filtration_type >::max() );
		
		//and create the cubical complex:	
		Cubical_complex* result = new Cubical_complex( sizes , top_dimensional_cells );
		
		std::vector< size_t > next_iteration_top_dimensional_cubes;
		
		//set the value 0 for all the cubes in top_dimensional_cubes_that_are_in_the_set:
		for ( size_t cube_no = 0 ; cube_no != top_dimensional_cubes_that_are_in_the_set.size() ; ++cube_no )
		{	
			size_t position = result->give_position_of_top_dimensional_cell( top_dimensional_cubes_that_are_in_the_set[cube_no] );
			result->get_cell_data(position) = 0;
			std::vector< size_t > neigs = result->give_neighbouring_top_dimensional_cells( position );
			next_iteration_top_dimensional_cubes.insert( next_iteration_top_dimensional_cubes.end() , neigs.begin() , neigs.end() );			
		}
		
		this->dylate_cubes_in_the_set( result , next_iteration_top_dimensional_cubes );
	
		return result;
	}//construct_cubical_complex_and_set_up_the_filtration_to_distance_from_selected_cubes

};

template <typenale Cubical_complex , Predicator>
Dylate_collection_of_cubes_in_cubical_complex::Dylate_collection_of_cubes_in_cubical_complex( Cubical_complex* to_process , const Predicator& pred )
{
	for ( auto it = cmplx->top_dimensional_cells_iterator_begin() ; it != top_dimensional_cells_iterator_end() ; ++it )
	{
		auto& cell_data = cmplx.get_cell_data(*it);
		if ( pred( cell_data ) )
		{
			cell_data = 0;			
		}
		else
		{
			cell_data = std::numeric_limits< typename Cubical_complex::filtration_type >::infinity();
		}		
	}	
}

 
}//namespace TOPOLOGICAL_INFERENCE_H_ 
}//namespace Gudhi 

#endif
