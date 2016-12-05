/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Pawel Dlotko
 *
 *    Copyright (C) 2016  INRIA Sophia-Saclay (France)
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

#ifndef COMPUTE_PERSISTENCE_WITH_PHAT_H_
#define COMPUTE_PERSISTENCE_WITH_PHAT_H_

//collumn representations in phat
#include <phat/representations/vector_vector.h>
#include <phat/representations/vector_heap.h>
#include <phat/representations/vector_set.h>
#include <phat/representations/vector_list.h>
#include <phat/representations/sparse_pivot_column.h>
#include <phat/representations/heap_pivot_column.h>
#include <phat/representations/full_pivot_column.h>
#include <phat/representations/bit_tree_pivot_column.h>

//algorithms
#include <phat/compute_persistence_pairs.h>
#include <phat/algorithms/standard_reduction.h>
#include <phat/algorithms/chunk_reduction.h>
#include <phat/algorithms/row_reduction.h>
#include <phat/algorithms/twist_reduction.h>
#include <phat/algorithms/spectral_sequence_reduction.h>


namespace Gudhi {

namespace phat_interface {
	
	using namespace std;

/**
 * \class Compute_persistence_with_phat
 * \brief Interface class to compute PHAT persistence from GUDHI data structures.
 * 
 * \ingroup phat_interface
 * 
 * \details
 * A class that creates PHAT boundary matrix based on GUDHI data structure and call one of four possible methods from
 * PHAT to compute \f$Z_{2}\f$ persistence. 
 * The template parameters are: 
 * GudhiFilteredComplex, a GUDHI data structure with the filtered complex (it can be for instance a simplex tree Simplex_tree, a
 * Bitmap_cubical_complex, and so on). A type of filtration can be recovered by taking GudhiFilteredComplex::Filtration_value. 
 * 
 * For more information on PHAT, please refer to \cite phat2017jsc.
 */
template <typename GudhiFilteredComplex , typename PHAT_reduction_algorithm = phat::twist_reduction , typename PHAT_collumn_representation = phat::bit_tree_pivot_column >
class Compute_persistence_with_phat {
 public:

  /** \brief Compute_persistence_with_phat constructor from a GUDHI data structure.
   *
   * Creates a PHAT boundary matrix based on the GUDHI data structure.
   * Notice that no computations of persistence are run by the constructor. To run persistence computations you need to
   * call a method
   * The second parameter of a constructor is a boolean variable, dualized which default value is set to false. In the dafault case the boundary matrix of a complex is created.
   * If the variable dualized is set to true, a coboundary matrix will be created.
   * Note also that the descructor of this class will not delete an object GudhiFilteredComplex* data_structure. The user is responsible for doing this.
   * 
   * @param[in] data_structure GUDHI data structure
   */
  Compute_persistence_with_phat(GudhiFilteredComplex* data_structure, bool dualized = false)
  : data_structure_(data_structure) {
	try
	{		
		this->dualized = dualized;
		size_t number_of_cells = this->data_structure_->num_simplices();		
		this->boundary_matrix_.set_num_cols( number_of_cells );		
						
		//We are creating a boundary matrix
		typename GudhiFilteredComplex::Filtration_simplex_range range = this->data_structure_->filtration_simplex_range();
		size_t position = 0;
		std::vector< phat::index > temp_col;
		for (typename GudhiFilteredComplex::Filtration_simplex_iterator it = range.begin(); it != range.end(); ++it) 
		{
		  //enumeration
		  this->data_structure_->assign_key(*it, position);
	  
	  
		  // setting up the dimension:
		  this->boundary_matrix_.set_dim(position, this->data_structure_->dimension(*it));			  			
	  
		  // setting up the boundary of this cell:
		  typename GudhiFilteredComplex::Boundary_simplex_range boundary_range = this->data_structure_->boundary_simplex_range(*it);
		  for (typename GudhiFilteredComplex::Boundary_simplex_iterator bd = boundary_range.begin(); bd != boundary_range.end(); ++bd) {
			temp_col.push_back(this->data_structure_->key(*bd));
		  }
		  // we do not know if the boundary elements are sorted according to filtration, that is why I am enforcing it here:
		  
		  
		  //Question: temp_col is a short vector. Does it make sense to use tbb::pararel_sort over here?
		  #ifdef GUDHI_USE_TBB
		  tbb::parallel_sort(temp_col.begin(), temp_col.end());
		  #else
		  std::stable_sort(temp_col.begin(), temp_col.end());
		  #endif		 				
		  
		  this->boundary_matrix_.set_col(position, temp_col);
		  temp_col.clear();
	  
		  ++position;
		}		
		if ( dualized )
		{
			cerr << "We are dualizing the matrix \n";
			phat::dualize ( this->boundary_matrix_ );
		}
	
    }
	catch (...)
	{
		  std::cout << "The constructor of a class Compute_persistence_with_phat thrown an exception/ Most probably the structure you are trying to use do not satisfy all the requirements of Gudhi FilteredComplex concept. \n";
		  throw "The constructor of a class Compute_persistence_with_phat thrown an exception/ Most probably the structure you are trying to use do not satisfy all the requirements of Gudhi FilteredComplex concept. \n";
	}
  }  

  /** \brief A function that call PHAT function compute_persistence_pairs. The PHAT algorithm used for the reduction of the matrix is specified as a second template of the Compute_persistence_with_phat 
   * class. The options avaialble at the moment are:
   * The "standard" algorithm (see [1], p.153)
   * The "row" algorithm from [2] (called pHrow in that paper)
   * The "twist" algorithm, as described in [3] (default algorithm)
   * The "chunk" algorithm presented in [4]
   * The "spectral sequence" algorithm (see [1], p.166) //TODO -- FIX THE REFEENCES!!!
   * Note that the function returns a structure od phat::persistence_pairs. That allows later to get vector of birth--death pairs (by using get_the_intervals method) as well as
   * getting the information about the cells that are paired. 
   *
   * @return The persistence pairs.
   */
  void compute_persistence_pairs() {
    phat::persistence_pairs pairs;
    phat::compute_persistence_pairs< PHAT_reduction_algorithm >(pairs, this->boundary_matrix_);
    if( this->dualized ) dualize_persistence_pairs( pairs, this->boundary_matrix_.get_num_cols() );
    this->pairs_ = pairs;
  }
 
  /** \brief Once the matrix is reduced, this function returns Betti numbers (as a vector of birth times of infinite
   * generators) and finite persistence intervals.
   *
   * \details
   * This function can be called only when the PHAT boundary matrix is reduced by one of four possible PHAT functions
   * to compute persistence.
   * 
   * @return The data returned from this function is a pair. 
   * The first element of the pair is: std::vector< std::vector<GudhiFilteredComplex::Filtration_value> >, this is a graded (by a dimension) vector of
   * beginning of infinite persistence intervals.
   * The second element of the pair is: std::vector< std::vector< std::pair<GudhiFilteredComplex::Filtration_value,GudhiFilteredComplex::Filtration_value> > > >, this is a graded (by a
   * dimension) vector of pairs<GudhiFilteredComplex::Filtration_value,GudhiFilteredComplex::Filtration_value>. Each such a pair is a beginning and end of finite persistence interval.
   */
  std::pair< std::vector< std::vector<typename GudhiFilteredComplex::Filtration_value> >, std::vector< std::vector< std::pair<typename GudhiFilteredComplex::Filtration_value , typename GudhiFilteredComplex::Filtration_value> > > >
  get_the_intervals() 
  {
	 try
	 {  
		bool dbg = false;				
		//first we need to check if the persistence pairs have been computed:
		if ( this->pairs_.get_num_pairs() == 0 )
		{
			std::cerr << "Persistence pairs have not been computed. Please call a function compute_persistence_pairs() before calling a function get_the_intervals(). The program will now termiante. \n";			
			throw "Persistence pairs have not been computed. Please call a function compute_persistence_pairs() before calling a function get_the_intervals(). The program will now termiante. \n";			
		}
		
		
		// in order to find the birth times of the infinite homology classes, we need to know which elements are not paired.
		// To search for them, we will use this vector:
		std::vector<bool> isTheElementPaired(this->data_structure_->num_simplices(), false);
	  
		// now it is time to recover the finite persistence pairs and the Betti numbers:
		std::vector< std::vector< std::pair<typename GudhiFilteredComplex::Filtration_value,typename GudhiFilteredComplex::Filtration_value> > > finitePersistencePairs(this->data_structure_->dimension());
		for (phat::index idx = 0; idx < this->pairs_.get_num_pairs(); idx++) {
		  typename GudhiFilteredComplex::Simplex_key positionOfBeginOfInterval = this->pairs_.get_pair(idx).first;
		  typename GudhiFilteredComplex::Simplex_key positionOfEndOfInterval = this->pairs_.get_pair(idx).second;
	  
		  typename GudhiFilteredComplex::Simplex_handle first_simplex = this->data_structure_->simplex(positionOfBeginOfInterval);
		  typename GudhiFilteredComplex::Simplex_handle second_simplex = this->data_structure_->simplex(positionOfEndOfInterval);
	  
		  typename GudhiFilteredComplex::Filtration_value valueFirst = this->data_structure_->filtration(first_simplex);
		  typename GudhiFilteredComplex::Filtration_value valueSecond = this->data_structure_->filtration(second_simplex);
	  
		  if (valueFirst > valueSecond) {
			std::swap(valueFirst, valueSecond);
		  }
	  
		  unsigned dimFirst = this->data_structure_->dimension(first_simplex);
		  unsigned dimSecond = this->data_structure_->dimension(second_simplex);
		  unsigned dim = std::min(dimFirst, dimSecond);
	  
	  
		  // we are ignoring trivial barcodes
		  if (valueFirst != valueSecond) {
			finitePersistencePairs[ dim ].push_back(std::make_pair(valueFirst, valueSecond));
			if (dbg) {
			  std::cerr << "Adding barcode : " << valueFirst << "," << valueSecond << std::endl;
			}
		  }
	  
		  isTheElementPaired[ this->pairs_.get_pair(idx).first ] = true;
		  isTheElementPaired[ this->pairs_.get_pair(idx).second ] = true;
		}
	  
	  
		std::vector< std::vector<typename GudhiFilteredComplex::Filtration_value> > birthTimesOfInfinitePersistnceClasses(this->data_structure_->dimension() + 1);
		for (size_t i = 0; i != size_t(this->data_structure_->dimension() + 1); ++i) {
		  std::vector<typename GudhiFilteredComplex::Filtration_value> v;
		  birthTimesOfInfinitePersistnceClasses[i] = v;
		}
		for (size_t i = 0; i != isTheElementPaired.size(); ++i) {
		  if (isTheElementPaired[i] == false) {
			// i-th element is not paired, therefore it gives an infinite class
			typename GudhiFilteredComplex::Simplex_handle simplex = this->data_structure_->simplex(i);
			birthTimesOfInfinitePersistnceClasses[this->data_structure_->dimension(simplex)].push_back(
				this->data_structure_->filtration(simplex));
		  }
		}
	  
		// sorting finite persistence pairs: 
		for (size_t dim = 0; dim != finitePersistencePairs.size(); ++dim) {
		  #ifdef GUDHI_USE_TBB
		  tbb::parallel_sort( finitePersistencePairs[dim].begin(), finitePersistencePairs[dim].end() );
		  #else
		  std::stable_sort( finitePersistencePairs[dim].begin(), finitePersistencePairs[dim].end() );
		  #endif      
		}
		return std::make_pair(birthTimesOfInfinitePersistnceClasses, finitePersistencePairs);
	}
	catch (...)
	{
		std::cout << "The method get_the_intervals() thrown an exception. The most probable reason for that is because the data structure you are trying to use do not satisfy all teh requirements of GudhiFilteredComplex concept.  \n";
		throw "The method get_the_intervals() thrown an exception. The most probable reason for that is because the data structure you are trying to use do not satisfy all teh requirements of GudhiFilteredComplex concept.  \n";		
	}	
  }


  /** \brief This function load a boundary matrix from a file and stores it as a boundary matrix in this data
   * structure.
   */
  void load_ascii(std::string filename) {
    this->boundary_matrix_.load_ascii(filename);
  }

  /** \brief This function store a boundary matrix in a PHAT format to a given file. The matrix can be later read by
   * load_ascii function.
   */
  void save_ascii(std::string filename) {
    std::ofstream out;
    out.open(filename.c_str());
    out << (long int) this->boundary_matrix_.get_num_cols() << std::endl;

    for (phat::index col_idx = 0; col_idx < this->boundary_matrix_.get_num_cols(); col_idx++) {
      out << (long int) this->boundary_matrix_.get_dim(col_idx) << " ";

      std::vector< phat::index > temp_col;
      this->boundary_matrix_.get_col(col_idx, temp_col);
      for (phat::index idx = 0; idx < (phat::index)temp_col.size(); idx++) {
        out << temp_col[ idx ] << " ";
      }
      out << std::endl;
    }


    out.close();
  }

  /** \brief This function is used for a debugging purposes only. It write the boundary matrix on standard output.
   */
  void print_boundary_matrix() {
    std::cout << "The boundary matrix has " << this->boundary_matrix_.get_num_cols() << " columns: " << std::endl;
    for (phat::index col_idx = 0; col_idx < this->boundary_matrix_.get_num_cols(); col_idx++) {
      std::cout << "Colum " << col_idx << " represents a cell of dimension "
          << (long int) this->boundary_matrix_.get_dim(col_idx) << ". ";
      if (!this->boundary_matrix_.is_empty(col_idx)) {
        std::vector< phat::index > temp_col;
        this->boundary_matrix_.get_col(col_idx, temp_col);
        std::cout << "Its boundary consists of the cells";
        for (phat::index idx = 0; idx < (phat::index)temp_col.size(); idx++)
          std::cout << " " << temp_col[ idx ];
      }
      std::cout << std::endl;
    }
  }



 private:
  phat::boundary_matrix< PHAT_collumn_representation > boundary_matrix_;
  GudhiFilteredComplex* data_structure_;
  phat::persistence_pairs pairs_;
  bool dualized;
};


/** \brief This function takes as an input the output of a Compute_persistence_with_phat procedure and stores it to a
 * file in a standard GUDHI format.
 */
template <typename K>
void write_intervals_to_file_Gudhi_format(std::pair< std::vector< std::vector<K> >,
                                         std::vector< std::vector< std::pair<K, K> > > > intervals,
                                         const char* filename, size_t dimension_cup = std::numeric_limits< size_t >::max()) {
  std::ofstream out;
  out.open(filename);
  //first goes the infinite intervals:
  int number_of_infinite_intervas = 0;
  size_t last_dimension_to_consider = intervals.first.size();
  // we need this line, since if we stop creating complex at some dimension, we will have a lot of holo simplices over
  // there which will show up as infinite generators here. We do not want to display them.
  if (dimension_cup != std::numeric_limits< size_t >::max())last_dimension_to_consider = dimension_cup;
  for (size_t dim = 0; dim != last_dimension_to_consider; ++dim) {
    number_of_infinite_intervas += intervals.first[dim].size();
  }
  std::vector< std::pair< K, size_t > > beginnings_of_infinite_intervals;
  beginnings_of_infinite_intervals.reserve(number_of_infinite_intervas);

  for (size_t dim = 0; dim != last_dimension_to_consider; ++dim) {
    for (size_t i = 0; i != intervals.first[dim].size(); ++i) {
      beginnings_of_infinite_intervals.push_back(std::make_pair(intervals.first[dim][i], dim));
    }
  }

  // now we need to sort beginnings_of_infinite_intervals according to the first coordinate:
  #ifdef GUDHI_USE_TBB
  tbb::parallel_sort(beginnings_of_infinite_intervals.begin(), beginnings_of_infinite_intervals.end()
  [](const std::pair<K, size_t>& lhs, const std::pair<K, size_t>& rhs) {
    return lhs.second > rhs.second; } );
  #else
  std::stable_sort(beginnings_of_infinite_intervals.begin(), beginnings_of_infinite_intervals.end(),
  [](const std::pair<K, size_t>& lhs, const std::pair<K, size_t>& rhs) {
    return lhs.second > rhs.second; });
  #endif
  //std::sort(beginnings_of_infinite_intervals.begin(), beginnings_of_infinite_intervals.end(),
  //          [](const std::pair<K, size_t>& lhs, const std::pair<K, size_t>& rhs) {
  //  return lhs.second > rhs.second; });

  // and now we output the sorted pairs to a file:
  for (size_t i = 0; i != beginnings_of_infinite_intervals.size(); ++i) {
    out << "2  " << beginnings_of_infinite_intervals[i].second << " " << beginnings_of_infinite_intervals[i].first
        << " inf" << std::endl;
  }


  // now it is time to deal with the finite intervals. We do very similar trick as above for them.
  int number_of_finite_intervas = 0;
  for (size_t dim = 0; dim != intervals.second.size(); ++dim) {
    number_of_finite_intervas += intervals.second[dim].size();
  }
  std::vector< std::pair< std::pair<K, K>, size_t > > finite_intervals;
  finite_intervals.reserve(number_of_finite_intervas);

  for (size_t dim = 0; dim != intervals.second.size(); ++dim) {
    for (size_t i = 0; i != intervals.second[dim].size(); ++i) {
      finite_intervals.push_back(std::make_pair(intervals.second[dim][i], dim));
    }
  }

  // and now we need to sort the finite_intervals vector according to the length of intervals, i.e. according to
  // finite_intervals[i].first.second - finite_intervals[i].first.first.
  #ifdef GUDHI_USE_TBB
  tbb::parallel_sort(finite_intervals.begin(), finite_intervals.end(), [](const std::pair<std::pair<K, K>, size_t>& lhs,
      const std::pair<std::pair<K, K>, size_t>& rhs) {
    return lhs.first.second - lhs.first.first > rhs.first.second - rhs.first.first; });
  #else
  std::stable_sort(finite_intervals.begin(), finite_intervals.end(), [](const std::pair<std::pair<K, K>, size_t>& lhs,
      const std::pair<std::pair<K, K>, size_t>& rhs) {
    return lhs.first.second - lhs.first.first > rhs.first.second - rhs.first.first; });
  #endif  
  //std::sort(finite_intervals.begin(), finite_intervals.end(), [](const std::pair<std::pair<K, K>, size_t>& lhs,
  //    const std::pair<std::pair<K, K>, size_t>& rhs) {
  //  return lhs.first.second - lhs.first.first > rhs.first.second - rhs.first.first; });



  // and now we should output them to a file:
  for (size_t i = 0; i != finite_intervals.size(); ++i) {
    out << "2  " << finite_intervals[i].second << " " << finite_intervals[i].first.first << " "
        << finite_intervals[i].first.second << std::endl;
  }

  out.close();
}

/** \brief A procedure that outputs a persistence to a files in a format alternative to a standard Gudhi format. It
 * creates separated files for Betti numbers, and separated files for persistence in each dimension. 
 * The persistence intervals in each dimension are not sorted.
 */
template <typename K>
void write_betti_numbers_and_persistence_intervals_to_file(const char* prefix,
                                                           std::pair< std::vector<std::vector< K > >,
                                                           std::vector< std::vector< std::pair<K, K> > > > resutsFromPhat) {
  std::ostringstream filenameStr;
  filenameStr << prefix << "_bettiNumbers";
  std::string str = filenameStr.str();
  const char* filename = str.c_str();
  std::ofstream out;
  out.open(filename);
  for (size_t dim = 0; dim != resutsFromPhat.first.size(); ++dim) {
    out << "Dimension : " << dim << std::endl;
    for (size_t i = 0; i != resutsFromPhat.first[dim].size(); ++i) {
      out << resutsFromPhat.first[dim][i] << std::endl;
    }
    out << std::endl;
  }
  out.close();


  for (size_t dim = 0; dim != resutsFromPhat.second.size(); ++dim) {
    if (resutsFromPhat.second[dim].size() == 0)continue;
    std::ostringstream filenameStr;
    filenameStr << prefix << "_persistence_" << dim;
    std::string str = filenameStr.str();
    const char* filename = str.c_str();
    std::ofstream out1;
    out1.open(filename);
    for (size_t i = 0; i != resutsFromPhat.second[dim].size(); ++i) {
      out1 << resutsFromPhat.second[dim][i].first << " " << resutsFromPhat.second[dim][i].second << std::endl;
    }
    out1.close();
  }
}

}  // phat_interface

}  // namespace Gudhi

#endif  // COMPUTE_PERSISTENCE_WITH_PHAT_H_
