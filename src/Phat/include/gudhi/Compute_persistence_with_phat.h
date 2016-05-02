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

#pragma once


#include "phat/compute_persistence_pairs.h"
#include "phat/representations/vector_vector.h"
#include "phat/algorithms/standard_reduction.h"
#include "phat/algorithms/chunk_reduction.h"
#include "phat/algorithms/row_reduction.h"
#include "phat/algorithms/twist_reduction.h"
#include "phat/algorithms/spectral_sequence_reduction.h"


namespace Gudhi
{


/*
 * A procedure that outputs a persistence to a files in a format alternative to a standard Gudhi format. It create separated files for Betti numbers, and separated files for persistence in each dimension. 
 * The persistence intervals in each dimension are not sorted.  
 * 
*/
template <typename K>
void writeBettiNumbersAndPersistenceIntervalsToFile( const char* prefix , std::pair< std::vector<std::vector< K > > , std::vector< std::vector< std::pair<K,K> > > > resutsFromPhat )
{
    std::ostringstream filenameStr;
    filenameStr << prefix << "_bettiNumbers";
    std::string str = filenameStr.str();
    const char* filename = str.c_str();
    std::ofstream out;
    out.open( filename );
    for ( size_t dim = 0 ; dim != resutsFromPhat.first.size() ; ++dim )
    {
        out << "Dimension : " << dim << std::endl;
        for ( size_t i = 0 ; i != resutsFromPhat.first[dim].size() ; ++i )
        {
            out << resutsFromPhat.first[dim][i] << std::endl;
        }
        out << std::endl;
    }
    out.close();


    for ( size_t dim = 0 ; dim != resutsFromPhat.second.size() ; ++dim )
    {        
        if ( resutsFromPhat.second[dim].size() == 0 )continue;
        std::ostringstream filenameStr;
        filenameStr << prefix << "_persistence_" << dim;
        std::string str = filenameStr.str();
        const char* filename = str.c_str();
        std::ofstream out1;
        out1.open( filename );
        for ( size_t i = 0 ; i != resutsFromPhat.second[dim].size() ; ++i )
        {
            out1 << resutsFromPhat.second[dim][i].first << " " << resutsFromPhat.second[dim][i].second << std::endl;
        }
        out1.close();
    }
}//writeBettiNumbersAndPersistenceIntervalsToFile


/*
 * A class that create Phat boundary matrix based on Gudhi data structure and call one of four possible methods from Phat to compute ZZ/2 persisteince. 
 * The template parameters are: 
 * T, a Gudhi data strucutre with the filtered complex (it can be for instance a simplex tree Simplex_tree or Bitmap_cubical_complex
 * The second parameter is a type of filtration (most typically double). 
*/
template <typename T , typename K>
class Compute_persistence_with_phat
{
public:
    /*
     * The only constructor of a Compute_persistence_with_phat class. It takes as an input a pointer to a Gudhi data structure, and based on it, create a phat boundary matrix. 
     * No computations of persistence are run by a constructor. To run persistence computations you need to choose one of four possibel methods to compute persistence in phat:
     * (1) compute_persistence_pairs_dualized_chunk_reduction();
     * (2) compute_persistence_pairs_twist_reduction();
     * (3) compute_persistence_pairs_standard_reduction();
     * (4) compute_persistence_pairs_spectral_sequence_reduction();
    */
    Compute_persistence_with_phat( T* data_structure_ );
    
    
    /*
     * A function that call Phat function compute_persistence_pairs_dualized_chunk_reduction.
    */
    phat::persistence_pairs compute_persistence_pairs_dualized_chunk_reduction();
    
    /*
     * A function that call Phat function compute_persistence_pairs_twist_reduction.
    */
    phat::persistence_pairs compute_persistence_pairs_twist_reduction();
    
    /*
     * A function that call Phat function compute_persistence_pairs_standard_reduction.
    */
    phat::persistence_pairs compute_persistence_pairs_standard_reduction();
    
    /*
     * A function that call Phat function compute_persistence_pairs_spectral_sequence_reduction.
    */
    phat::persistence_pairs compute_persistence_pairs_spectral_sequence_reduction();


	/*
	 * This function can be called only when the Phat boundary matrix is reduced by one of four possible Phat functions to compute persistence. Once the matrix is reduced, this function return 
	 * Betti numbers (as a vector of birth times of infinite generators) and finite persistence intervls.
	 * The data returned from this function is a pair. 
	 * The first element of the pair is: std::vector< std::vector<K> >, this is a graded (by a dimension) vector of beginning of infinite persistence intervals.
	 * The second element of the pair is: std::vector< std::vector< std::pair<K,K> > > >, this is a graded (by a dimension) vector of pairs<K,K>. Each such a pair is a beginning and end of a persistence interval.
	*/
    std::pair< std::vector< std::vector<K> > , std::vector< std::vector< std::pair<K,K> > > >  get_the_intervals( phat::persistence_pairs pairs );
    
	 /*
	  * This function store a boundary matrix in a phat format to a given file. For some reason this format is not recognized by load_ascii procedure in Phat. Therefore one can store the matrix with this procedure
	  * but will not be later able to read it with load_ascii. If you want to read it later, please use save_for_reading. 
	 */
     void save_ascii( std::string filename )
     {
         this->boundary_matrix.save_ascii( filename );
     }
     
     /*
	  * This function store a boundary matrix in a phat format to a given file. The matrix can be later read by load_ascii finction.
	 */
     void save_for_reading( std::string filename )
     {
         std::ofstream out;
         out.open( filename.c_str() );
         out << (int)this->boundary_matrix.get_num_cols() << std::endl;
        
         for( phat::index col_idx = 0; col_idx < this->boundary_matrix.get_num_cols(); col_idx++ ) 
         {
             out  << (int)this->boundary_matrix.get_dim( col_idx ) << " ";
			 
            std::vector< phat::index > temp_col;
            this->boundary_matrix.get_col( col_idx, temp_col );            
            for( phat::index idx = 0; idx < (phat::index)temp_col.size(); idx++ )
            {
                out << temp_col[ idx ] << " ";
			}
			out << std::endl;
		 }
    
         
         out.close();
     }
     
     /*
	  * This function load a boundary matrix from a file and stores it as a boundary matrix in this data structure. 
	 */
     void load_ascii( std::string filename )
     {
            this->boundary_matrix.load_ascii( filename );
     }
     
     /*
	  * This function is used for a debugging purposes only. It write the boundary matrix on a screen.
	 */
     void print_bd_matrix();

     
private:
    phat::boundary_matrix< phat::vector_vector > boundary_matrix;
    T* data_structure;
};

template <typename T , typename K>
void Compute_persistence_with_phat<T,K>::print_bd_matrix()
{
    std::cout << "The boundary matrix has " << this->boundary_matrix.get_num_cols() << " columns: " << std::endl;
    for( phat::index col_idx = 0; col_idx < this->boundary_matrix.get_num_cols(); col_idx++ ) {
        std::cout << "Colum " << col_idx << " represents a cell of dimension " << (int)this->boundary_matrix.get_dim( col_idx ) << ". ";
        if( !this->boundary_matrix.is_empty( col_idx ) ) {
            std::vector< phat::index > temp_col;
            this->boundary_matrix.get_col( col_idx, temp_col );
            std::cout << "Its boundary consists of the cells";
            for( phat::index idx = 0; idx < (phat::index)temp_col.size(); idx++ )
                std::cout << " " << temp_col[ idx ];
        }
        std::cout << std::endl;
    }
}

template <typename T , typename K>
phat::persistence_pairs Compute_persistence_with_phat<T,K>::compute_persistence_pairs_dualized_chunk_reduction()
{
    phat::persistence_pairs pairs;
    phat::compute_persistence_pairs_dualized< phat::chunk_reduction >( pairs, this->boundary_matrix );
    return pairs;
}

template <typename T , typename K>
phat::persistence_pairs Compute_persistence_with_phat<T,K>::compute_persistence_pairs_twist_reduction()
{
    phat::persistence_pairs pairs;
    phat::compute_persistence_pairs< phat::twist_reduction >( pairs, this->boundary_matrix );
    return pairs;
}

template <typename T , typename K>
phat::persistence_pairs Compute_persistence_with_phat<T,K>::compute_persistence_pairs_standard_reduction()
{	
    phat::persistence_pairs pairs;
    phat::compute_persistence_pairs< phat::standard_reduction >( pairs, this->boundary_matrix );
    return pairs;
}

template <typename T , typename K>
phat::persistence_pairs Compute_persistence_with_phat<T,K>::compute_persistence_pairs_spectral_sequence_reduction()
{
    phat::persistence_pairs pairs;
    phat::compute_persistence_pairs< phat::spectral_sequence_reduction  >( pairs, this->boundary_matrix );
    return pairs;
}

template <typename T , typename K>
Compute_persistence_with_phat<T,K>::Compute_persistence_with_phat( T* data_structure_ ):data_structure( data_structure_ )
{
    bool dbg = false;
   
    this->boundary_matrix.set_num_cols( this->data_structure->num_simplices() );
    
    
    //due to interface of phom in Gudhi it seems that I have to fill-in the keys for all simplices here. They are sorting according to filtration, but keys are not filled in (for some reason).
    //in the same time when doing the enumeration, we wills set up the dimensions and setting up the boundary of cells.
    typename T::Filtration_simplex_range range = this->data_structure->filtration_simplex_range();
    size_t position = 0;
    std::vector< phat::index > temp_col;
    for ( typename T::Filtration_simplex_iterator it = range.begin() ; it != range.end() ; ++it )
    {
		//enumeration
		this->data_structure->assign_key( *it , position );
	
		
		//setting up the dimension:
		this->boundary_matrix.set_dim( position, this->data_structure->dimension( *it ) );
		
		//setting up the boundary of this cell:
		typename T::Boundary_simplex_range boundary_range = this->data_structure->boundary_simplex_range( *it );
        for ( typename T::Boundary_simplex_iterator bd = boundary_range.begin() ; bd != boundary_range.end() ; ++bd )
        {
            temp_col.push_back( this->data_structure->key( *bd ) );
        }     
        //we do not know if the boundary elements are sorted according to filtration, that is why I am enforcing it here:
        std::sort( temp_col.begin() , temp_col.end() );
        this->boundary_matrix.set_col( this->data_structure->key( *it ) , temp_col );
        temp_col.clear();         
		
		++position;
	}
	std::cout << "Done enumerating the simplices \n";
}

template <typename T , typename K>
std::pair< std::vector< std::vector<K> > , std::vector< std::vector< std::pair<K,K> > > >  Compute_persistence_with_phat<T,K>::get_the_intervals( phat::persistence_pairs pairs )
{
    bool dbg = false;
    //in order to find the birth times of the infinite homology classes, we need to know which elements are not paired. To search for them, we will use this vector:
    std::vector<bool> isTheElementPaired( this->data_structure->num_simplices() , false );

    //now it is time to recover the finite persistence pairs and the Betti numbers:
    std::vector< std::vector< std::pair<K,K> > > finitePersistencePairs( this->data_structure->dimension() );
    for( phat::index idx = 0; idx < pairs.get_num_pairs(); idx++ )
    {
        typename T::Simplex_key positionOfBeginOfInterval = pairs.get_pair( idx ).first;
        typename T::Simplex_key positionOfEndOfInterval = pairs.get_pair( idx ).second;

        typename T::Simplex_handle first_simplex = this->data_structure->simplex(positionOfBeginOfInterval);
        typename T::Simplex_handle second_simplex = this->data_structure->simplex(positionOfEndOfInterval);

        typename T::Filtration_value valueFirst = this->data_structure->filtration( first_simplex );
        typename T::Filtration_value valueSecond = this->data_structure->filtration( second_simplex );

        if ( valueFirst > valueSecond ){std::swap( valueFirst , valueSecond );}

        unsigned dimFirst = this->data_structure->dimension(first_simplex);
        unsigned dimSecond = this->data_structure->dimension(second_simplex);
        unsigned dim = std::min( dimFirst , dimSecond );


        //we are ignoring trivial barcodes
        if ( valueFirst != valueSecond )
        {
            finitePersistencePairs[ dim ].push_back( std::make_pair(valueFirst , valueSecond) );
            if ( dbg ){std::cerr << "Adding barcode : " << valueFirst << "," << valueSecond << std::endl;}
        }

        isTheElementPaired[ pairs.get_pair( idx ).first ] = true;
        isTheElementPaired[ pairs.get_pair( idx ).second ] = true;
    }


    std::vector< std::vector<K> > birthTimesOfInfinitePersistnceClasses(this->data_structure->dimension()+1 );
    for ( size_t i = 0 ; i != this->data_structure->dimension()+1 ; ++i )
    {
        std::vector<K> v;
        birthTimesOfInfinitePersistnceClasses[i] = v;
    }
    for ( size_t i = 0 ; i != isTheElementPaired.size() ; ++i )
    {
        if ( isTheElementPaired[i] == false )
        {
            //i-th element is not paired, therefore it gives an infinite class
            typename T::Simplex_handle simplex = this->data_structure->simplex(i);
            birthTimesOfInfinitePersistnceClasses[ this->data_structure->dimension( simplex ) ].push_back( this->data_structure->filtration(simplex) );
        }
    }

    //sorting finite persistence pairs: 
    for ( size_t dim = 0 ; dim != finitePersistencePairs.size() ; ++dim )
    {
        std::sort( finitePersistencePairs[dim].begin() , finitePersistencePairs[dim].end() );
    }
    return std::make_pair( birthTimesOfInfinitePersistnceClasses , finitePersistencePairs );
}//Compute_persistence_with_phat


/*
 * This function takes as an input the output of a Compute_persistence_with_phat procedure and stores it to a file in a strandard Gudhi format.
*/ 
template <typename K>
void write_intervas_to_file_Gudhi_format( std::pair< std::vector< std::vector<K> > , std::vector< std::vector< std::pair<K,K> > > > intervals , const char* filename , size_t dimension_cup = -1 )
{
	std::	ofstream out;
	out.open( filename );
	//first goes the infinite intervals:
	int number_of_infinite_intervas = 0;
	size_t last_dimension_to_consider = intervals.first.size();
	//we need this line, since if we stop creating compelx at some dimension, we will have a lot of holo simplices over there which will show up as infinite generators here. We do not want to display them.
	if ( dimension_cup != -1 )last_dimension_to_consider = dimension_cup;
	for ( size_t dim = 0  ; dim != last_dimension_to_consider ; ++dim )
	{
		number_of_infinite_intervas += intervals.first[dim].size();
	}
	std::vector< std::pair< K , size_t > > beginnings_of_infinite_intervals;
	beginnings_of_infinite_intervals.reserve( number_of_infinite_intervas );
	
	for ( size_t dim = 0  ; dim != last_dimension_to_consider ; ++dim )
	{
		for ( size_t i = 0 ; i != intervals.first[dim].size() ; ++i )
		{
			beginnings_of_infinite_intervals.push_back( std::make_pair( intervals.first[dim][i] , dim ) );
		}
	}
	
	//now we need to sort beginnings_of_infinite_intervals according to the first coordinate:
	std::sort(beginnings_of_infinite_intervals.begin(), beginnings_of_infinite_intervals.end(), [](const std::pair<K , size_t>& lhs, const std::pair<K , size_t>& rhs) {return lhs.second > rhs.second; } );
	
	//and now we output the sorded pairs to a file:
	for ( size_t i = 0 ; i != beginnings_of_infinite_intervals.size() ; ++i )
	{
		out << "2  " << beginnings_of_infinite_intervals[i].second << " " <<  beginnings_of_infinite_intervals[i].first << " inf" << std::endl; 
	} 	
	
	
	//now it is time to deal with the finite intervals. We do very similar trick as above for them.
	int number_of_finite_intervas = 0;
	for ( size_t dim = 0  ; dim != intervals.second.size() ; ++dim )
	{
		number_of_finite_intervas += intervals.second[dim].size();
	}
	std::vector< std::pair< std::pair<K,K> , size_t > > finite_intervals;
	finite_intervals.reserve( number_of_finite_intervas );
	
	for ( size_t dim = 0  ; dim != intervals.second.size() ; ++dim )
	{
		for ( size_t i = 0 ; i != intervals.second[dim].size() ; ++i )
		{
			finite_intervals.push_back( std::make_pair( intervals.second[dim][i] , dim ) );
		}
	}
	
	//and now we need to sort the finite_intervals vector according to the length of intervals, i.e. according to finite_intervals[i].first.second - finite_intervals[i].first.first.
	std::sort(finite_intervals.begin(), finite_intervals.end(), [](const std::pair<std::pair<K,K> , size_t>& lhs, const std::pair<std::pair<K,K> , size_t>& rhs) {return lhs.first.second - lhs.first.first > rhs.first.second - rhs.first.first; } );

	//and now we should output them to a file:
	for ( size_t i = 0 ; i != finite_intervals.size() ; ++i )
	{
		out << "2  " << finite_intervals[i].second << " " <<  finite_intervals[i].first.first << " " << finite_intervals[i].first.second << std::endl; 
	}
	
	out.close();
}

}//namespace Gudhi
