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


namespace Gudhi
{


//the only aim of this class is to have a ability to compute persistence with phat.
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


template <typename T , typename K>
class Compute_persistence_with_phat
{
public:
    Compute_persistence_with_phat( T* data_structure_ );
    std::pair< std::vector< std::vector<K> > , std::vector< std::vector< std::pair<K,K> > > >  get_the_intervals( phat::persistence_pairs pairs );
    
    
    //this constructor shoud be removed, since when for a class constructed with it we will call get_the_intervals, it will all crash. 
    Compute_persistence_with_phat()
    {
		this->data_structure = 0;
	}

    phat::persistence_pairs compute_persistence_pairs_dualized_chunk_reduction();
    phat::persistence_pairs compute_persistence_pairs_twist_reduction();
    phat::persistence_pairs compute_persistence_pairs_standard_reduction();//for some reason, this do not always work.
    //phat::persistence_pairs compute_persistence_pairs_spectral_sequence_reduction();

     void save_ascii( std::string filename )
     {
         this->boundary_matrix.save_ascii( filename );
     }
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
     void load_ascii( std::string filename )
     {
            this->boundary_matrix.load_ascii( filename );
     }
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

//template <typename T , typename K>
//phat::persistence_pairs Compute_persistence_with_phat<T,K>::compute_persistence_pairs_spectral_sequence_reduction()
//{
//    phat::persistence_pairs pairs;
//    phat::compute_persistence_pairs< phat::spectral_sequence_reduction >( pairs, this->boundary_matrix );
//    return pairs;
//}

template <typename T , typename K>
Compute_persistence_with_phat<T,K>::Compute_persistence_with_phat( T* data_structure_ ):data_structure( data_structure_ )
{
    bool dbg = true;
    
    if ( dbg )
    {
		std::cerr << "Number of simplices : " << this->data_structure->num_simplices() << std::endl;
		getchar();
	}
    this->boundary_matrix.set_num_cols( this->data_structure->num_simplices() );

    //setting up the dimensions of cells:
    for ( size_t i = 0 ; i != this->data_structure->num_simplices() ; ++i )
    {
		if ( dbg )std::cerr << "dimension of a simplex number : " << i << " " << this->data_structure->dimension( this->data_structure->simplex(i) ) << std::endl;
        this->boundary_matrix.set_dim( i, this->data_structure->dimension( this->data_structure->simplex(i) ) );
    }
    
    std::cerr << "OUT \n:";


    //now it is time to set up the boundary matrix:
    int aa = 0;
    
    typename T::Filtration_simplex_range range = this->data_structure->filtration_simplex_range();
    std::vector< phat::index > temp_col;
    for ( typename T::Filtration_simplex_iterator it = range.begin() ; it != range.end() ; ++it )
    {
		if ( dbg )
		{
			std::cerr << " this->data_structure->key( *it )  : " <<  this->data_structure->key( *it )  << " and its boundary : " << std::endl;
		}
		
        typename T::Boundary_simplex_range boundary_range = this->data_structure->boundary_simplex_range( *it );
        for ( typename T::Boundary_simplex_iterator bd = boundary_range.begin() ; bd != boundary_range.end() ; ++bd )
        {
			if ( dbg )
			{
				std::cerr << this->data_structure->key( *bd ) << std::endl;
			}
            temp_col.push_back( this->data_structure->key( *bd ) );
        }
        getchar();
        //we do not know if the boundary elements are sorted according to filtration, that is why I am enforcing it here:
        std::sort( temp_col.begin() , temp_col.end() );
        this->boundary_matrix.set_col( this->data_structure->key( *it ) , temp_col );
        temp_col.clear();
        
        
        
        std::cerr << "aa : " << aa << std::endl;
        ++aa;
    }
    
    /*
    std::vector< phat::index > temp_col;
    for ( size_t i = 0 ; i != this->data_structure->num_simplices() ; ++i )
    {
        typename T::Boundary_simplex_range boundary_range = this->data_structure->boundary_simplex_range( i );
        for ( typename T::Boundary_simplex_iterator bd = boundary_range.begin() ; bd != boundary_range.end() ; ++bd )
        {
            temp_col.push_back( this->data_structure->key( *bd ) );
        }
        //we do not know if the boundary elements are sorted according to filtration, that is why I am enforcing it here:
        std::sort( temp_col.begin() , temp_col.end() );
        this->boundary_matrix.set_col( i , temp_col );
        temp_col.clear();
    }
	*/
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
std::vector< std::pair< std::vector< std::pair<size_t , size_t> > , std::vector< std::pair<size_t , size_t> > > > return_cubes_paired_in_peristence( phat::persistence_pairs pairs )
{
	std::vector< std::pair< std::vector< std::pair<size_t , size_t> > , std::vector< std::pair<size_t , size_t> > > > result;
	result.reserve( pairs.get_num_pairs() );
	
	for( phat::index idx = 0; idx < pairs.get_num_pairs(); idx++ )
	{
		size_t positionOfBeginOfInterval = pairs.get_pair( idx ).first;
        size_t positionOfEndOfInterval = pairs.get_pair( idx ).second;
        
        std::vector<unsigned> counter1 = b->compute_counter_for_given_cell( positionOfBeginOfInterval );
        std::vector<unsigned> counter2 = b->compute_counter_for_given_cell( positionOfEndOfInterval );
        
        std::vector< std::pair<size_t , size_t> > frist_cube;
        frist_cube.reserve( counter1.size() );
        for ( size_t i = 0 ; i != counter1.size() ; ++i )
        {
			if ( counter1[i]%2 == 1 )
			{
				first_cube.push_back( std::make_pair(counter1[i]%2 , counter1[i]%2+1) );
			}
			else
			{
				first_cube.push_back( std::make_pair(counter1[i]%2 , counter1[i]%2) );
			}
		}
		
		std::vector< std::pair<size_t , size_t> > second_cube;
        second_cube.reserve( counter2.size() );
        for ( size_t i = 0 ; i != counter2.size() ; ++i )
        {
			if ( counter2[i]%2 == 1 )
			{
				second_cube.push_back( std::make_pair(counter2[i]%2 , counter2[i]%2+1) );
			}
			else
			{
				second_cube.push_back( std::make_pair(counter2[i]%2 , counter2[i]%2) );
			}
		}
		
		result.push_back( std::make_pair( frist_cube , second_cube ) );
	}
	
	
	return result;
}//return_cubes_paired_in_peristence
*/


}//namespace Gudhi
