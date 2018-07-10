/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Clément Maria
 *
 *    Copyright (C) 2018 Inria
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



// cmd f "todo" for transverse row system


#ifndef _ZIGZAG_PERSISTENT_HOMOLOGY_
#define _ZIGZAG_PERSISTENT_HOMOLOGY_

#include <cmath>
#include <iomanip>

#include <boost/tuple/tuple.hpp>
#include <boost/intrusive/set.hpp>
#include <boost/pending/disjoint_sets.hpp>
#include <boost/intrusive/list.hpp>
#include <boost/pool/object_pool.hpp>
#include <boost/progress.hpp>

namespace Gudhi {

namespace zigzag_persistence {

template <typename FilteredComplex>
class matrix_chain;

struct zzmat_h_tag; // for horizontal traversal in the persistence matrix
struct zzmat_v_tag; // for vertical traversal in the persistence matrix
typedef boost::intrusive::list_base_hook 
< boost::intrusive::tag < zzmat_h_tag >
, boost::intrusive::link_mode < boost::intrusive::auto_unlink >
>      base_hook_zzmat_h; //allows .unlink()
typedef boost::intrusive::list_base_hook 
< boost::intrusive::tag < zzmat_v_tag >
, boost::intrusive::link_mode < boost::intrusive::auto_unlink >
>      base_hook_zzmat_v; //faster hook, less safe todo really?

//----------------------------------------------------------------------------------
/* 
 * Cell for the persistence matrix. Contains a key for the simplex index, and
 * horizontal and vertical hooks for connections within sparse rows and columns.
 */
template < typename FilteredComplex >
class Zigzag_persistence_cell 
	: public base_hook_zzmat_h, public base_hook_zzmat_v {
public:
    typedef typename FilteredComplex::Simplex_key              Simplex_key;
    typedef matrix_chain< FilteredComplex >                    Matrix_chain;

    Zigzag_persistence_cell( Simplex_key     key
			     , Matrix_chain * self_chain )
	: key_(key)
	, self_chain_(self_chain) {}

    Simplex_key     key_;
    Matrix_chain * self_chain_;
};

//----------------------------------------------------------------------------------
// /*
//  * Sparse column in the persistence matrix. Essentially a wrapper for a 
//  * boost::intrusive::list.
//  */
// template < typename FilteredComplex >
// class Zigzag_persistence_column
// : public boost::intrusive::set_base_hook 
//              < boost::intrusive::link_mode< boost::intrusive::normal_link > > 
// {
// public:
//   typedef typename FilteredComplex::Simplex_key              Simplex_key;   
//   typedef Zigzag_persistence_cell < FilteredComplex >        Cell;

//   typedef boost::intrusive::list < Cell 
//                                  , boost::intrusive::constant_time_size<false> 
//                                  , boost::intrusive::base_hook< base_hook_zzmat_v > 
//                                  >               Col_type;
// /** \brief Creates an empty column.*/
//   Zigzag_persistence_column () { col_ = Col_type(); }
//   ~Zigzag_persistence_column () {};

//   Col_type           col_;
// };

//----------------------------------------------------------------------------------
/*
 * Chain in the persistence matrix. It has access to a
 * - column that represent the chain as a sum of simplices,
 * - a row
 * - is paired with another chain, indicating its type F G H
 * - has a direct access to its lowest index.
 */
template < typename FilteredComplex >
class matrix_chain {
public:
    typedef typename FilteredComplex::Simplex_key              Simplex_key;
    typedef typename FilteredComplex::Simplex_handle           Simplex_handle;
    typedef typename FilteredComplex::Filtration_value         Filtration_value;
    // Encoding Matrix types:
    // // Column type
    // typedef Zigzag_persistence_column < FilteredComplex >  Column; // contains 1 set_hook
    // // Cell type
    // typedef typename Column::Cell                     Cell;   // contains 2 list_hooks
    typedef Zigzag_persistence_cell < FilteredComplex > Cell; // contains 2 list_hooks
    typedef boost::intrusive::list < Cell
    , boost::intrusive::constant_time_size<false>
    , boost::intrusive::base_hook< base_hook_zzmat_v >
    >                    Column;//vertical list of cells
    // Remark: constant_time_size must be false because base_hook_cam_h has auto_unlink link_mode
    typedef boost::intrusive::list < Cell
    , boost::intrusive::constant_time_size<false>
    , boost::intrusive::base_hook< base_hook_zzmat_h >
    >                  Row;  //horizontal list of cells
    matrix_chain()
	: column_(NULL), row_(NULL), paired_col_(NULL)
	, birth_(-3),
	  // birth_fil_(0),
	  lowest_idx_(-1) {}

    matrix_chain( Column         * c
		  , Row            * r
		  , matrix_chain   * p_col
		  , Simplex_key      b
		  // , Filtration_value fil
		  , Simplex_key      l)
	: column_(c), row_(r), paired_col_(p_col)
	, birth_(b),
	  // birth_fil_(fil),
	  lowest_idx_(l) {}

    matrix_chain * paired_col() { return paired_col_; }
    void assign_paired_col( matrix_chain * other_col ) { paired_col_ = other_col; }
    Simplex_key birth() { return birth_; }
    void assign_birth( Simplex_key b ) { birth_ = b; }
    void assign_birth(matrix_chain *other)
    { birth_ = other->birth_;
	// birth_fil_ = other->birth_fil_;
    }
    // void assign_fil(Filtration_value fil) { birth_fil_ = fil; }
    // private:
    Column             * column_      ; //col at index i
    Row                * row_         ; //row at index i
    matrix_chain       * paired_col_  ; //\in F -> NULL, \in H -> g, \in G -> h
    Simplex_key          birth_       ; //\in F -> b, \in H -> -2 \in G -> -1
    // Filtration_value     birth_fil_   ; //filtration value of simplex with key birth_
    Simplex_key          lowest_idx_  ; //lowest_idx_ = i

};

// Print a Simplex_tree in os.
template <typename... T>
std::ostream& operator<<(std::ostream& os, matrix_chain<T...>& chain) {
    os << " : [k=" << chain.lowest_idx_ <<"] ";
    os << "[b=" << chain.birth_ <<"] ";
    if(chain.paired_col_ != NULL) { os << "[pc=" << chain.paired_col_->lowest_idx_ <<"]     ";}
    else { os << "[pc=" << "F" <<"]     ";}
    os <<"| ";
    for(auto &cell : *(chain.column_)) { os << cell.key_ << " "; }
    os << " | ";
    return os;
}

//----------------------------------------------------------------------------------
/**
 * \class Zigzag_persistence Zigzag_persistence gudhi/Zigzag_persistence.h
 * \brief Methods to compute the zigzag persistent homology of a zigzag
 * filtered complex.
 *
 * The type ZigzagComplex::Simplex_key counts the number of insertions and
 * deletions of simplices, which may be large in zigzag persistence and require
 * more than 32 bits of storage. The type used (int, long, etc) should be chosen in
 * consequence. Simplex_key must be signed.
 */
template < typename ZigzagComplex >
class Zigzag_persistence {
public:
    typedef ZigzagComplex                                 Complex_ds;
    // Data attached to each simplex to interface with a Property Map.
    typedef typename Complex_ds::Simplex_key              Simplex_key;//must be signed
    typedef typename Complex_ds::Simplex_handle           Simplex_handle;
    typedef typename Complex_ds::Filtration_value         Filtration_value;
    typedef matrix_chain< Complex_ds >                    Matrix_chain;
    typedef typename Matrix_chain::Cell                   Cell;
    typedef typename Matrix_chain::Column                 Column;
    typedef typename Matrix_chain::Row                    Row;

private:
    /*Structure to store persistence intervals. By convention, interval [b;d] are
     *closed for finite indices b and d, and open for left-infinite and/or
     *right-infinite endpoints.*/
    struct interval_t {
	interval_t() {}
	interval_t(int dim, double b, double d) : dim_(dim), b_(b), d_(d) {}
	double length() {
	    if (b_ == d_) return 0; //o.w. inf - inf would return nan.
	    if (b_ > d_) return b_ - d_;
	    else return d_ - b_;
	    //return abs(d_ - b_);
	}
	int    dim_;	//dimension
	double b_;	//birth index
	double d_;	//death index
    };
    //Comparison function to sort intervals by length in the output persistence diagram
    struct cmp_intervals_by_log_length {
	cmp_intervals_by_log_length(){}
	bool operator()( interval_t p, interval_t q)
	{
	    if ((double)p.d_ > (double)p.b_){
		if ((double)q.d_ > (double)q.b_) {
		    return log2((double)p.d_) - log2((double)p.b_) > log2((double)q.d_) - log2((double)q.b_);
		} else {
		    return log2((double)p.d_) - log2((double)p.b_) > log2((double)q.b_) - log2((double)q.d_);
		}
	    } else {
		if ((double)q.d_ > (double)q.b_) {
		    return log2((double)p.b_) - log2((double)p.d_) > log2((double)q.d_) - log2((double)q.b_);
		} else {
		    return log2((double)p.b_) - log2((double)p.d_) > log2((double)q.b_) - log2((double)q.d_);
		}
	    }
	    /*return abs(log2((double)p.d_) - log2((double)p.b_)) >
		    abs(log2((double)q.d_) - log2((double)q.b_));*/
	}
    };
    struct cmp_intervals_by_length {
	cmp_intervals_by_length(){}
	bool operator()( interval_t p, interval_t q)
	{
	    if(p.dim_ != q.dim_) {return p.dim_ < q.dim_;}
	    if(p.length() != q.length()) { return p.length() > q.length(); } //longest first
	    if(p.b_ != q.b_) {return p.b_ < q.b_;}
	    // if(p.d_ != q.d_) return p.d_ < q.d_;
	    return p.d_ < q.d_;
	}
    };

public:


    Zigzag_persistence( Complex_ds &cpx )
	: cpx_(&cpx)
	, lowidx_to_matidx_()
	, matrix_()
	, birth_ordering_()
	, persistence_diagram_()
	, num_arrow_(0)
	// , prev_fil_(0)
	// , curr_fil_(0)
	, filtration_values_()
    {}

private:
    // /* Swap the columns .col_ stored in two matrix_chains. */
    // void col_swap(matrix_chain * curr_col, matrix_chain * other_col)
    // {
    //   std::swap(curr_col->column_ , other_col->column_ );
    //   for(auto &cell : curr_col->column_->col_)  { cell.self_chain_ = curr_col; }
    //   for(auto &cell : other_col->column_->col_) { cell.self_chain_ = other_col;}
    // }
    /*
 * Set c1 <- c1 + c2, assuming canonical order of indices induced by the order in
 * the vertical lists. self1 is the matrix_chain whose column is c1, for self
 * reference of the new cells.
 */
    void plus_equal_column(Matrix_chain * self1, Column * c1, Column * c2)
    {
	auto it1 = c1->begin();   auto it2 = c2->begin();
	while(it1 != c1->end() && it2 != c2->end())
	{
	    if(it1->key_ < it2->key_) { ++it1; }
	    else {
		if(it1->key_ > it2->key_) {
		    Cell * new_cell = new Cell(it2->key_, self1);
		    c1->insert(it1, *new_cell); //col link, in order
		    lowidx_to_matidx_[it2->key_]->row_->push_back(*new_cell);//row link,no order
		    ++it2;
		}
		else { //it1->key_ == it2->key_
		    auto tmp_it = it1;    ++it1; ++it2;
		    Cell * tmp_ptr = &(*tmp_it);
		    tmp_it->base_hook_zzmat_h::unlink(); //unlink from row
		    c1->erase(tmp_it); //remove from col
		    delete tmp_ptr;
		}
	    }
	}
	while(it2 != c2->end()) { //if it1 reached the end of its column, but not it2
	    Cell * new_cell = new Cell(it2->key_,self1);
	    lowidx_to_matidx_[it2->key_]->row_->push_back(*new_cell); //row links
	    c1->push_back(*new_cell);
	    ++it2;
	}
    }


    /*
 * Set c1 <- c1 + c2, assuming canonical order of indices induced by the order in
 * the vertical lists. self1 is the matrix_chain whose column is c1, for self
 * reference of the new cells.
 */
    // void plus_equal_column(matrix_chain * self1, Column * c1, std::set<Simplex_key> &c2)
    // {
    // auto it1 = c1->begin();   auto it2 = c2->begin();
    // while(it1 != c1->end() && it2 != c2->end())
    // {
    //   if(it1->key_ < it2->key_) { ++it1; }
    //   else {
    //     if(it1->key_ > it2->key_) {
    //       Cell * new_cell = new Cell(it2->key_, self1);
    //       c1->insert(it1, *new_cell); //col link, in order
    //       lowidx_to_matidx_[it2->key_]->row_->push_back(*new_cell);//row link,no order
    //       ++it2;
    //     }
    //     else { //it1->key_ == it2->key_
    //       auto tmp_it = it1;    ++it1; ++it2;
    //       Cell * tmp_ptr = &(*tmp_it);
    //       tmp_it->base_hook_zzmat_h::unlink(); //unlink from row
    //       c1->erase(tmp_it); //remove from col
    //       delete tmp_ptr;
    //     }
    //   }
    // }
    // while(it2 != c2->end()) { //if it1 reached the end of its column, but not it2
    //   Cell * new_cell = new Cell(it2->key_,self1);
    //   lowidx_to_matidx_[it2->key_]->row_->push_back(*new_cell); //row links
    //   c1->push_back(*new_cell);
    //   ++it2;
    // }
    // }

    /*
 * Set c1 <- c1 + c2, assuming canonical order of indices induced by the order in
 * the vertical lists. self1 is the matrix_chain whose column is c1, for self
 * reference of the new cells.
 */
    // void plus_equal_setcol(matrix_chain * self1, Column * c1, std::set<Simplex_key> & c2)
    // {
    // auto it1 = c1->begin();   auto it2 = c2->begin();
    // while(it1 != c1->end() && it2 != c2->end())
    // {
    //   if(it1->key_ < it2->key_) { ++it1; }
    //   else {
    //     if(it1->key_ > it2->key_) {
    //       Cell * new_cell = new Cell(it2->key_, self1);
    //       c1->insert(it1, *new_cell); //col link, in order
    //       lowidx_to_matidx_[it2->key_]->row_->push_back(*new_cell);//row link,no order
    //       ++it2;
    //     }
    //     else { //it1->key_ == it2->key_
    //       auto tmp_it = it1;    ++it1; ++it2;
    //       Cell * tmp_ptr = &(*tmp_it);
    //       tmp_it->base_hook_zzmat_h::unlink(); //unlink from row
    //       c1->erase(tmp_it); //remove from col
    //       delete tmp_ptr;
    //     }
    //   }
    // }
    // while(it2 != c2->end()) { //if it1 reached the end of its column, but not it2
    //   Cell * new_cell = new Cell(it2->key_,self1);
    //   lowidx_to_matidx_[it2->key_]->row_->push_back(*new_cell); //row links
    //   c1->push_back(*new_cell);
    //   ++it2;
    // }
    // }




    /* Maintains the <=b order of indices.*/
    // struct birth_vector {
    //   birth_vector() : inv_b_vec(), maxb(0), minb(0) {}

    //   void add_birth_forward() { //forward arrow
    //     inv_b_vec.push_back(maxb+1); ++maxb;
    //   }
    //   void add_birth_backward() { //backward arrow
    //     inv_b_vec.push_back(minb-1); --minb;
    //   }
    //   int operator[](Simplex_key k) { return inv_b_vec[k]; }
    //   //give priority to a Simplex_key inv_b_vec[key] == the number of smaller birth
    //   std::vector<int> inv_b_vec;
    //   int              maxb;
    //   int              minb;
    // };

    /* Maintains the birth ordering <=b. Contains an std::map of size the number of
 * non-zero rows of the homology matrix.
 */
    struct birth_ordering {
	//example quiver indices    empty_cpx -> 0 -> 1 -> 2 <- 3 <- 4 -> 5 <- 6 etc
	birth_ordering() : birth_to_pos_(), max_birth_pos_(0), min_birth_pos_(-1) {}

	void add_birth_forward(Simplex_key key) { //amortized constant time
	    birth_to_pos_.emplace_hint(birth_to_pos_.end(), key, max_birth_pos_);
	    ++max_birth_pos_;
	}
	void add_birth_backward(Simplex_key key) { //amortized constant time
	    birth_to_pos_.emplace_hint(birth_to_pos_.end(), key, min_birth_pos_);
	    --min_birth_pos_;
	}
	//when the row at index key is removed from the homology matrix
	void remove_birth(Simplex_key key) { birth_to_pos_.erase(key); }
	//increasing birth order <=b, true iff b1 <b b2
	bool birth_order(Simplex_key k1, Simplex_key k2) {
	    return birth_to_pos_[k1] < birth_to_pos_[k2];
	}
	//decreasing birth order <=b, true iff b1 >b b2
	bool reverse_birth_order(Simplex_key k1, Simplex_key k2) {
	    // std::cout << "X\n";
	    return birth_to_pos_[k1] > birth_to_pos_[k2];
	}

	void display() {
	    for(auto pp : birth_to_pos_) {
		std::cout << "(" << pp.first << ";" << pp.second << ") ";
	    }
	    std::cout << std::endl;
	}

	// private:
	//given a birth index, starting at 0, returns the position in the birth ordering <=b
	//over the entire quiver.
	std::map< Simplex_key, Simplex_key > birth_to_pos_;
	Simplex_key                          max_birth_pos_;
	Simplex_key                          min_birth_pos_;
    };

public:
    /**
  * \brief Compute the zigzag persistent homology of a zigzag filtered complex,
  * using the reflection and transposition algorithm of [Maria, Oudot '14].
  *
  * matrix_, originally empty, maintains the set of chains, with a partition F.G.H,
  * representing a compatible homology basis as in [Maria, Oudot '14].
  *
  * Each simplex type in the complex stores a .key_ field that stores the index of
  * its insertion in the zigzag filtration.
  *
  * The algorithm maintains a compatible homology basis for the zigzag filtration
  * \emptyset = K_0 <-> (...) <-> K_i <- . <- . <- ... <- \emptyset
  * where the prefix from K_0 to K_i is equal to the i-prefix of the input zigzag
  * filtration given by cpx_->filtration_simplex_range(), and the suffix (from K_i
  * on to the right end) is a sequence of simplex removals. Due to the structure of * reflection diamonds, the removal are in reverse order of the insertions.
  *
  * Consequently, using cpx_->key(zzsh) as indexing for the matrix rows/cells,
  * with the natural order on integers, makes our homology matrix matrix_ upper
  * triangular for the suffix K_i <- ... <- 0, seen as a standard persistence
  * filtration. At K_i, the natural order on integers is also equivalent to the
  * death-order <=d (because all arrows in the suffix are backward).
  *
  *
  */
    void compute_zigzag_persistence()
    { //compute index persistence, interval are closed, i.e., [b,d) is stored as [b,d-1]
	// std::vector< std::pair< Simplex_key, Filtration_value > > filtration_values;
	Filtration_value                                          prev_fil_, curr_fil_;

	assert(num_arrow_ == 0); //start with an 'empty' complex.
	auto zzrg = cpx_->filtration_simplex_range();
	auto zzit = zzrg.begin();

	double nber_forward = 0;
	double nber_backward = 0;
	double nber_halfpair = 0;
	std::chrono::time_point<std::chrono::system_clock> start, end;
	std::chrono::time_point<std::chrono::system_clock> start2, end2;
	std::chrono::duration<double> enlapsed_sec;
	double enlapsed_total = 0.0;

	/* TO DO num_arrows must be == cpx_->key(zzsh) for non-contiguous keys! */

	prev_fil_ = cpx_->filtration(*zzit);
	filtration_values_.emplace_back(cpx_->key(*zzit), prev_fil_);

	while( zzit != zzrg.end() )
	{
	    start2 = std::chrono::system_clock::now();
	    //std::cout << cpx_->key(*zzit) << " key\n";

	    /*for (Matrix_chain mc : matrix_){
		std::cout << mc.column_->empty() << " " << mc.row_->empty() << " ";
		if (mc.paired_col_ == NULL) std::cout << "NULL ";
		else std::cout << "paired ";
		std::cout << mc.birth_ << " " << mc.lowest_idx_ << "\n";
	    }
	    std::cout << "\n";*/

	    //num_arrow_ = cpx_->key(*zzit);

	    //if(num_arrow_ % 100000 == 0) std::cout << num_arrow_ << "\n";
	    // display_mat();
	    // std::cout << std::endl;
	    // if(zzit.arrow_direction()) std::cout << "+ ";
	    // else std::cout << "- ";
	    // for(auto v : cpx_->simplex_vertex_range(*zzit)) {
	    //   std::cout << v << " ";
	    // }
	    // std::cout << "      k" << cpx_->key(*zzit)  << "  f" << cpx_->filtration(*zzit) <<  "\n";
	    // std::cout << std::endl;



	    //gives only critical simplices for insertion, and potentially maximal non-critical (i.e., paired) simplices at deletion.
	    //keys must be already assigned by the filtration_simplex_iterator
	    if(zzit.arrow_direction()) //forward arrow
	    {
		//std::cout << "forward arrow\n";
		start = std::chrono::system_clock::now();
		++nber_forward;
		curr_fil_ = cpx_->filtration(*zzit);//check whether the filt val has changed
		if(curr_fil_ != prev_fil_)
		{ prev_fil_ = curr_fil_;
		    filtration_values_.emplace_back(cpx_->key(*zzit), prev_fil_);}

		forward_arrow(*zzit);
		end = std::chrono::system_clock::now();
		enlapsed_sec = end-start;
		//std::cout << "forward arrow: " << enlapsed_sec.count() << " sec.\n";
	    }
	    else //backward arrow
	    {
		++nber_backward;
		//std::cout << "backward arrow\n";
		if(!cpx_->is_critical(*zzit)) { //if the simplex is critical
		    //matrix A becomes matrix A U \{\tau,sigma\}
		    ++nber_halfpair;
		    start = std::chrono::system_clock::now();
		    make_pair_critical(*zzit);
		    end = std::chrono::system_clock::now();
		    enlapsed_sec = end-start;
		    //std::cout << "halfpair arrow: " << enlapsed_sec.count() << " sec.\n";
		}
		start = std::chrono::system_clock::now();
		backward_arrow(*zzit);
		end = std::chrono::system_clock::now();
		enlapsed_sec = end-start;
		//std::cout << "backward arrow: " << enlapsed_sec.count() << " sec.\n";
	    }
	    // }
	    // else {
	    //   auto tmpsh = *zzit; ++zzit;
	    //   if(!cpx_->is_pair(tmpsh,*zzit)) {std::cout << "Error \n"; return; }
	    // }
	    end2 = std::chrono::system_clock::now();
	    enlapsed_sec = end2-start2;
	    enlapsed_total += enlapsed_sec.count();

	    ++zzit;
	    ++num_arrow_; //same as simplex index, count all simplices, even non criticals
	}
	if(!matrix_.empty()) {std::cout << "Remains " << matrix_.size() << " columns.\n";}
	std::cout << "total: " << enlapsed_total << " sec.\n";

	//Compute the right-open intervals, for the remaining columns in the matrix.
	// for(auto & col : matrix_)
	// {
	//   switch(col.birth()) {
	//     case -2:{ break; } //in H
	//     case -1:{ //in G       //todo reverse birth and death?
	//       persistence_diagram_.emplace_back( cpx_->dimension(col.lowest_idx_)
	//                                        , cpx_->filtration(col.paired_col()->lowest_idx_)
	//                                        , cpx_->filtration(col.lowest_idx_));
	//       break;
	//     }
	//     default:{
	//       persistence_diagram_.emplace_back( cpx_->dimension(col.lowest_idx_)
	//                                        , cpx_->filtration(col.birth_)
	//                                        , log2(0));//todo gestion infinity in GUDHI
	//       break;
	//     }
	//   }
	// }

	//std::cout << "Total number of arrows: " << num_arrow_+1 << std::endl;
	std::cout << std::setprecision(std::numeric_limits<double>::digits10 + 1) << "Number of forward arrows: " << nber_forward << std::endl;
	std::cout << "Number of backward arrows: " << nber_backward << std::endl;
	std::cout << "Number of halfpair arrows: " << nber_halfpair << std::setprecision(5) << std::endl;

    }

    /* sh is a maximal simplex paired with a simplex tsh
 * Morse pair (tau,sigma)
 *
 * sh must be paired with tsh, i.e., tsh = cpx_->morse_pair(zzsh); return the
 * handle for tau, and cpx_->is_critical(zzsh) is false. The Morse iterator is in
 * charge of modifying this (make the simplices critical) after the
 * Zigzag_persistence update has been done.
 *
 *
 * key(tsh) < key(sh) must be true.
 *
 * cpx_->boundary_simplex_range(zzsh) must iterate along the boundary of sigma in
 * the Morse complex A' where sigma and tau have become critical
 *
 * cpx_->coboundary_simplex_range(tsh) must iterate along the cofaces of tsh in the
 * Morse complex A' where sigma and tau have become critical, i.e., those critical
 * faces \nu such that [\nu : \tau]^{A'} != 0.
 *
 */
    void make_pair_critical(Simplex_handle zzsh)
    {
	//std::cout << cpx_->key(zzsh) << " is half pair\n";

	/*for (Matrix_chain mc : matrix_){
	    std::cout << mc.column_->empty() << " " << mc.row_->empty() << " ";
	    if (mc.paired_col_ == NULL) std::cout << "NULL ";
	    else std::cout << "paired ";
	    std::cout << mc.birth_ << " " << mc.lowest_idx_ << "\n";
	}
	std::cout << "\n";*/

	auto tsh = cpx_->morse_pair(zzsh);//Morse pair (*tsh, *sh)
	//new column and row for sigma
	Column * new_col_s = new Column();
	Row  * new_row_s   = new Row();
	matrix_.emplace_front( new_col_s
			       , new_row_s
			       , (Matrix_chain *)0
			       , -2 //in H, paired w/ the new column for tau
			       // , curr_fil_
			       , cpx_->key(zzsh) );
	auto chain_it_s = matrix_.begin();
	//Add the bottom coefficient for zzsh
	Cell * new_cell_s = new Cell(cpx_->key(zzsh), &(*chain_it_s));
	new_col_s->push_back( *new_cell_s ); //zzsh has largest idx of all
	new_row_s->push_back( *new_cell_s );
	//Update the map 'index idx -> chain with lowest index idx' in matrix_
	lowidx_to_matidx_[cpx_->key(zzsh)] = chain_it_s;

	//new column and row for tau
	Column * new_col_t = new Column();
	Row *  new_row_t   = new Row();
	matrix_.emplace_front( new_col_t
			       , new_row_t
			       , (Matrix_chain *)0
			       , -1 //in G, paired w/ the new column for sigma
			       // , curr_fil_
			       , cpx_->key(tsh) );
	auto chain_it_t = matrix_.begin();
	//Update the map 'index idx -> chain with lowest index idx' in matrix_
	lowidx_to_matidx_[cpx_->key(tsh)] = chain_it_t;

	//pair the two new columns
	chain_it_t->assign_paired_col(&(*chain_it_s));
	chain_it_s->assign_paired_col(&(*chain_it_t));



	if(chain_it_s->lowest_idx_ != cpx_->key(zzsh)) {std::cout << "Error lowest key\n";}


	//fill up col_tau with \partial \sigma in new Morse complex
	std::set< Simplex_key > col_bsh; //set maintains the order on indices
	for( auto b_sh : cpx_->paired_simplex_boundary_simplex_range(zzsh) )//<-\partial in Morse complex
	{ col_bsh.insert(cpx_->key(b_sh)); }
	//copy \partial sigma in the new row&column
	for( auto idx : col_bsh ) //in increasing idx order
	{ //add all indices in col_tau with canonical order enforced by set col_bsh
	    Cell * new_cell = new Cell(idx, &(*chain_it_t));
	    new_col_t->push_back( *new_cell ); //insertion in column
	    lowidx_to_matidx_[idx]->row_->push_back( *new_cell ); //insertion in row
	}

	//update the row for sigma. First record all possible modified columns
	std::map<Matrix_chain *, int> modif_chain;
	for(auto c_sh : cpx_->paired_simplex_coboundary_simplex_range(tsh)) {//[*c_sh:*t_sh]^{A'} != 0
	    //all chains with != 0 index at c_sh
	    /*if(lowidx_to_matidx_.find(cpx_->key(c_sh)) == lowidx_to_matidx_.end()) {
		std::cout << "AAA\n";
		std::cout << (cpx_->key(c_sh)) << "  " << cpx_->key(zzsh) << "  " << cpx_->key(tsh) << " <--- \n";
	    }*/
	    if (cpx_->key(c_sh) != cpx_->key(zzsh)){
		for(auto cell : *(lowidx_to_matidx_[cpx_->key(c_sh)]->row_)) {
		    auto res_insert = modif_chain.emplace(cell.self_chain_,1);
		    if(!res_insert.second) {//already there
			++(res_insert.first->second); //one more occurrence of the chain
		    }
		}
	    }
	}
	//all chains appearing an odd number of times
	for(auto modif : modif_chain) {
	    if((modif.second % 2) == 1) { //sum_{nu \in chain} [nu:tau]^{A'} = 1
		//Add the bottom coefficient for zzsh to rectify new boundary
		Cell * new_cell = new Cell(cpx_->key(zzsh), modif.first);
		// modif.first->column_->push_back( *new_cell ); //zzsh doesn't have largest idx of all
		insert_cell(modif.first->column_, new_cell);
		new_row_s->push_back( *new_cell );//row for sigma
		//if chain in H, modify the paired c_g
		/*if(modif.first->birth() == -2) {
		    auto chain_g = modif.first->paired_col(); //c_g to modify <- add c_tau
		    plus_equal_column(chain_g, chain_g->column_, new_col_t);
		}*/
	    }//else sum == 0
	}

	/*std::cout << "End make_pair\n";
	for (Matrix_chain mc : matrix_){
	    std::cout << mc.column_->empty() << " " << mc.row_->empty() << " ";
	    if (mc.paired_col_ == NULL) std::cout << "NULL ";
	    else std::cout << "paired ";
	    std::cout << mc.birth_ << " " << mc.lowest_idx_ << "\n";
	    std::cout << "column: ";
	    for (Cell cell : *(mc.column_)) std::cout << cell.key_ << " ";
	    std::cout << "\n";
	}
	std::cout << "\n";*/
    }

    //insert a cell in the middle of a column
    void insert_cell(Column *c, Cell *cell) {
	auto it = c->begin();
	while(it != c->end() && it->key_ < cell->key_) {++it;}
	c->insert(it,*cell);
    }


    Filtration_value index_to_filtration(Simplex_key k) {
	auto it =
		std::upper_bound( filtration_values_.begin(), filtration_values_.end()
				  , std::pair<Simplex_key, Filtration_value>(k, std::numeric_limits<double>::infinity() )
				  , []( std::pair<Simplex_key, Filtration_value> p1
				  , std::pair<Simplex_key, Filtration_value> p2) {
		return p1.first < p2.first; }
		);
	if(it->first != k) {return (--it)->second;}
	else {return it->second;}
    }

    /** \brief Output the persistence diagram in ostream.
   *
   * The file format is the following:
   *    p1*...*pr   dim b d
   *
   * where "dim" is the dimension of the homological feature,
   * b and d are respectively the birth and death of the feature and
   * p1*...*pr is the product of prime numbers pi such that the homology
   * feature exists in homology with Z/piZ coefficients.
   */
    void output_diagram(std::ostream& ostream = std::cout) {

	int num_intervals = 5;

	std::cout << "Filtration values: ";
	for(auto pp : filtration_values_) {
	    std::cout << "[ " << pp.first << " ; " << pp.second << " ]  ";
	} std::cout << std::endl;

	std::vector< interval_t > tmp_diag;
	tmp_diag.reserve(persistence_diagram_.size());
	for(auto bar : persistence_diagram_) {
	    tmp_diag.emplace_back(bar.dim_,index_to_filtration(bar.b_),index_to_filtration(bar.d_));
	}
	cmp_intervals_by_length cmp;
	std::stable_sort(tmp_diag.begin(), tmp_diag.end(), cmp);

	if(tmp_diag.empty()) {return;}

	int curr_dim = tmp_diag.begin()->dim_;
	int curr_num_intervals = num_intervals;

	for(auto bar : tmp_diag) {
	    if(curr_dim != bar.dim_) {
		std::cout << "----------------------------------------- dim " << bar.dim_
			  << " \n";
		curr_num_intervals = num_intervals; curr_dim = bar.dim_;
	    }
	    if(curr_num_intervals > 0) {
		--curr_num_intervals;
		std::cout << bar.dim_ << "   " << bar.b_ << " " << bar.d_ <<
			     "      " << bar.length() << " \n";
	    }
	}
    }

    void output_log2_diagram(std::ostream& ostream = std::cout) {

	int num_intervals = 5;

	// std::cout << "Filtration values: ";
	// for(auto pp : filtration_values_) {
	//   std::cout << "[ " << pp.first << " ; " << pp.second << " ]  ";
	// } std::cout << std::endl;

	std::vector< interval_t > tmp_diag;
	tmp_diag.reserve(persistence_diagram_.size());
	for(auto bar : persistence_diagram_) {
	    tmp_diag.emplace_back(bar.dim_,log2(index_to_filtration(bar.b_)),log2(index_to_filtration(bar.d_)));
	}
	cmp_intervals_by_length cmp;
	std::stable_sort(tmp_diag.begin(), tmp_diag.end(), cmp);

	if(tmp_diag.empty()) {return;}

	int curr_dim = tmp_diag.begin()->dim_;
	int curr_num_intervals = num_intervals;

	for(auto bar : tmp_diag) {
	    if(curr_dim != bar.dim_) {
		std::cout << "----------------------------------------- dim " << bar.dim_
			  << " \n";
		curr_num_intervals = num_intervals; curr_dim = bar.dim_;
	    }
	    if(curr_num_intervals > 0) {
		--curr_num_intervals;
		std::cout << bar.dim_ << "   " << bar.b_ << " " << bar.d_ <<
			     "      " << bar.length() << " \n";
	    }
	}
    }


    /** \brief Output the persistence diagram in ostream.
   *
   * The file format is the following:
   *    p1*...*pr   dim b d
   *
   * where "dim" is the dimension of the homological feature,
   * b and d are respectively the birth and death of the feature and
   * p1*...*pr is the product of prime numbers pi such that the homology
   * feature exists in homology with Z/piZ coefficients.
   */
    void output_index_diagram(std::ostream& ostream = std::cout) {
	std::vector< interval_t > tmp_diag;
	tmp_diag.reserve(persistence_diagram_.size());
	for(auto bar : persistence_diagram_) {
	    tmp_diag.emplace_back(bar.dim_,bar.b_,bar.d_);
	}
	// cmp_intervals_by_length cmp;
	// std::stable_sort(tmp_diag.begin(), tmp_diag.end(), cmp);

	for(auto bar : tmp_diag) {
	    std::cout << bar.dim_ << "   " << bar.b_ << " " << bar.d_ << "      " << bar.length() << " \n";
	}
    }


    /** \brief Output the persistence diagram in ostream.
  *
  * The file format is the following:
  *    p1*...*pr   dim b d
  *
  * where "dim" is the dimension of the homological feature,
  * b and d are respectively the birth and death of the feature and
  * p1*...*pr is the product of prime numbers pi such that the homology
  * feature exists in homology with Z/piZ coefficients.
  */
    // void output_diagram(std::ostream& ostream = std::cout) {
    //   persistence_diagram_.sort(cmp_intervals_by_log_length());
    //   for(auto inter_ : persistence_diagram_) {
    //       ostream << inter_.dim_ << " "
    //               << log2((double)inter_.d_) << " " << log2((double)inter_.b_)
    //               << std::endl;
    //   }
    // }



    void display_mat() {
	std::cout << "---------------------------beg \n";
	for(auto & col : matrix_)
	{
	    // if(col.paired_col_ == NULL && col.birth() != col.lowest_idx_) {
	    //   std::cout << "birth != lowest_idx !!\n";
	    // } <---- birth != lowest_idx
	    std::cout << " : [k=" << col.lowest_idx_ <<"] ";
	    std::cout << "[b=" << col.birth_ <<"] ";
	    if(col.paired_col_ != NULL) { std::cout << "[pc=" << col.paired_col_->lowest_idx_ <<"]     ";}
	    else { std::cout << "[pc=" << "F" <<"]     ";}
	    std::cout <<"| ";
	    for(auto &cell : *(col.column_)) { std::cout << cell.key_ << " "; }
	    std::cout << " |" << std::endl;
	}
	std::cout << "----------- \n";
	// std::cout << "Bv: "; for(auto b : birth_vector_.inv_b_vec) { std::cout << b << " "; } std::cout << "\n";
	// std::cout << "Bv: "; for(auto pp : birth_ordering_.birth_to_pos_) { std::cout << "("<<pp.first<<";"<<pp.second<<") ";}
	// std::cout << "    maxb = " << birth_ordering_.max_birth_pos_ << " , minb = " << birth_ordering_.min_birth_pos_ << "\n";
	std::cout << "Pd: "; for(auto bd: persistence_diagram_)
	{ std::cout << "["<<bd.b_<<";"<<bd.d_<< "  d" << bd.dim_ << "] "; } std::cout << "\n";
	// std::cout << "l_to_m: "; for(auto ltm : lowidx_to_matidx_)
	//                       { std::cout << ltm.first << "->" << ltm.second->lowest_idx_ << " ";} std::cout << "\n";
	for(auto ltm : lowidx_to_matidx_) {
	    if(ltm.first != ltm.second->lowest_idx_) {std::cout << "ltm != lowest_idx ??\n";}
	}
	std::cout << "---------------------------end \n";
    }




private: 
    /**
  * \brief Computes the boundary cycle of the new simplex zzsh, and express it as a
  * sum of cycles. If all cycles are boundary cycles, i.e., columns with G-index
  * in the matrix, then [\partial zzsh] = 0 and we apply an injective diamond to
  * the zigzag module. Otherwise, we keep reducing with boundary and live cycles,
  * i.e., columns with (F \cup G)-indices, and then apply a surjective diamond to
  * the zigzag module.
  */
    void forward_arrow( Simplex_handle zzsh )
    { //maintain the <=b order
	birth_ordering_.add_birth_forward(cpx_->key(zzsh));//num_arrow_);
	//Reduce the boundary of zzsh in the basis of cycles.
	//Compute the simplex keys of the simplices of the boundary of zzsh.
	std::set< Simplex_key > col_bsh; //set maintains the order on indices
	for( auto b_sh : cpx_->boundary_simplex_range(zzsh) )
	{ col_bsh.insert(cpx_->key(b_sh)); }
	// {  col_bsh.insert( b_sh ); } //todo prob here!

	// std::cout << "col_bsh : ";
	// for(auto k : col_bsh) {std::cout << k << " ";}
	// std::cout << std::endl;


	//If empty boundary (e.g., if zzsh is a vertex in a simplicial complex)
	if(col_bsh.empty())
	{
	    // std::cout << "     --inj shortcut\n";
	    //New row and column with a bottom-right non-zero element, at index key(zzsh)
	    Column * new_col  = new Column();
	    Row    * new_row  = new Row();
	    matrix_.emplace_front( new_col
				   , new_row
				   , (Matrix_chain *)0 //in F, paired with NULL
				   , cpx_->key(zzsh)
				   // , curr_fil_
				   , cpx_->key(zzsh) );
	    auto chain_it = matrix_.begin();
	    Cell   * new_cell = new Cell(cpx_->key(zzsh), &(*chain_it) );
	    new_col->push_back( *new_cell ); //zzsh must have largest idx of the column
	    new_row->push_back( *new_cell );
	    //Update the map 'index idx -> chain with lowest index idx' in matrix_
	    lowidx_to_matidx_[cpx_->key(zzsh)] = chain_it;
	    return;
	}

	// col_bsh.rbegin()) is idx of lowest element in col_bsh, because it is a set.
	Matrix_chain * col_low = &(*lowidx_to_matidx_[*(col_bsh.rbegin())]);

	// std::cout << "7=" << *(col_bsh.rbegin()) << " " << "7=" << col_low->lowest_idx_ <<"\n";

	// auto col_low   =  &(*lowidx_to_matidx_[*(col_bsh.rbegin())]); //<-- prob here
	auto paired_idx = col_low->paired_col(); //col with which col_low is paired
	std::vector< Matrix_chain * > chains_in_H; //for corresponding indices in H
	std::vector< Matrix_chain * > chains_in_G;

	//Reduce col_bsh with boundary cycles, i.e., indices in G.
	std::pair< typename std::set< Simplex_key >::iterator, bool > res_insert;
	while( paired_idx != NULL ) {
	    chains_in_H.push_back(paired_idx);//keep the col_h with which col_g is paired
	    chains_in_G.push_back(col_low);   //keep the col_g
	    for(auto &cell : *(col_low->column_)) { //Recuce with the column col_g
		res_insert = col_bsh.insert(cell.key_);
		if( !res_insert.second ) { col_bsh.erase(res_insert.first); } //1+1 = 0
		//o.w. insertion has succeeded.
	    }
	    //If col_bsh is entirely reduced, \partial zzsh is a boundary cycle.
	    if(col_bsh.empty()) {
		// if(cpx_->dimension(zzsh) >= max_dim_) {return;} we need max_dim creators
		injective_reflection_diamond(zzsh, chains_in_H);
		return;
	    }
	    //Continue the reduction
	    col_low     =  &(*lowidx_to_matidx_[*(col_bsh.rbegin())]);//curr low index col
	    paired_idx  =  col_low->paired_col();//col with which col_low is paired

	}


	// std::cout << "Columns in G: \n";
	// for(auto col : chains_in_G)
	// {
	//   std::cout << " : [k=" << col->lowest_idx_ <<"] ";
	//   std::cout << "[b=" << col->birth_ <<"] ";
	//   if(col->paired_col_ != NULL) { std::cout << "[pc=" << col->paired_col_->lowest_idx_ <<"]     ";}
	//   else { std::cout << "[pc=" << "F" <<"]     ";}
	//   std::cout <<"| ";
	//   for(auto &cell : *(col->column_)) { std::cout << cell.key_ << " "; }
	//   std::cout << " |" << std::endl;
	// }


	//Continue reducing with boundary and 'live' cycles, i.e., indices in G U F.
	std::vector< Matrix_chain * > chains_in_F;
	while(true)
	{
	    if (col_low->column_->empty()) {
		throw std::out_of_range("col_low->column_->empty() " + std::to_string(col_low->birth_) + " "
					+ std::to_string(col_low->lowest_idx_) + " " + std::to_string(*(col_bsh.rbegin())));
	    }

	    if(paired_idx == NULL) { chains_in_F.push_back(col_low); } //col_low \in F
	    else                   { chains_in_H.push_back(paired_idx); } //paired_idx \in H
	    //Reduce with the column col_g or col_f
	    for(auto &cell : *(col_low->column_)) {
		//std::cout << "col_low->column non empty\n";
		res_insert = col_bsh.insert(cell.key_);
		if( !res_insert.second ) { col_bsh.erase(res_insert.first); } //1+1 = 0
		//o.w. insertion has succeeded.
	    }
	    //If col_bsh is entirely reduced, i.e. col_bsh == \emptyset.

	    /*for (Simplex_key k : col_bsh) std::cout << k << " ";
	    std::cout << "\n";*/

	    if(col_bsh.empty()){
		//std::cout << "col_low->column empty\n";
		surjective_reflection_diamond(zzsh, chains_in_F, chains_in_H);
		return;
	    }
	    //Else, keep reducing.
	    col_low     =  &(*lowidx_to_matidx_[*(col_bsh.rbegin())]); //curr low index col
	    paired_idx  =  col_low->paired_col_;//col with which col_low is paired
	}
    }


    /**
 * \brief Compute an injective diamond in the zigzag module, by inserting a new
 * column for the chain   zzsh + \sum col_h, for all col_h in chains_in_H, and a
 * new row for the simplex zzsh.
 */
    inline
    void injective_reflection_diamond ( Simplex_handle                  zzsh
					, std::vector< Matrix_chain * > & chains_in_H
					)//, int dim_sh)
    {
	// std::cout << "     --inj\n";
	//Compute the chain   zzsh + \sum col_h, for col_h \in chains_in_H
	std::set< Simplex_key > col_bsh;
	std::pair< typename std::set< Simplex_key >::iterator, bool > res_insert;
	//produce the sum of all col_h in chains_in_H
	for( Matrix_chain * idx_h : chains_in_H ) {
	    for(auto &cell : *(idx_h->column_) ) {
		res_insert = col_bsh.insert(cell.key_);
		if( !res_insert.second ) { col_bsh.erase(res_insert.first); }
	    }
	}
	//Create a new row&column in the matrix
	Column * new_col = new Column();
	Row  * new_row   = new Row();
	matrix_.emplace_front( new_col
			       , new_row
			       , (Matrix_chain *)0 //in F
			       , cpx_->key(zzsh)
			       // , curr_fil_
			       , cpx_->key(zzsh) );
	auto chain_it = matrix_.begin();
	//copy the sum of col_h in chains_in_H in the new row&column
	for( auto idx : col_bsh ) //in increasing idx order
	{ //add all indices in col_new with canonical order enforced by set col_bsh
	    Cell * new_cell = new Cell(idx, &(*chain_it));
	    new_col->push_back( *new_cell ); //insertion in column
	    lowidx_to_matidx_[idx]->row_->push_back( *new_cell ); //insertion in row
	}
	//Add the bottom coefficient for zzsh
	Cell * new_cell = new Cell(cpx_->key(zzsh), &(*chain_it));
	new_col->push_back( *new_cell ); //zzsh has largest idx of all
	new_row->push_back( *new_cell );
	//Update the map 'index idx -> chain with lowest index idx' in matrix_
	lowidx_to_matidx_[cpx_->key(zzsh)] = chain_it;
    }

    /**
 * The vector chains_in_F is sorted by decreasing lowest index values in the
 * columns corresponding to the chains, due to its computation in the reduction of
 * \partial zzsh in forward_arrow(...). It is equivalent to decreasing death index
 * order w.r.t. the <d ordering.
 */
    inline
    void surjective_reflection_diamond( Simplex_handle zzsh
					, std::vector< Matrix_chain * > & chains_in_F
					, std::vector< Matrix_chain * > & chains_in_H )
    {
	// std::cout << "     --surj\n";

	// std::cout << "Columns in H: \n";
	// for(auto col : chains_in_H)
	// {
	//   std::cout << " : [k=" << col->lowest_idx_ <<"] ";
	//   std::cout << "[b=" << col->birth_ <<"] ";
	//   if(col->paired_col_ != NULL) { std::cout << "[pc=" << col->paired_col_->lowest_idx_ <<"]     ";}
	//   else { std::cout << "[pc=" << "F" <<"]     ";}
	//   std::cout <<"| ";
	//   for(auto &cell : *(col->column_)) { std::cout << cell.key_ << " "; }
	//   std::cout << " |" << std::endl;
	// }

	// std::cout << "Columns in F: \n";
	// for(auto col : chains_in_F)
	// {
	//   std::cout << " : [k=" << col->lowest_idx_ <<"] ";
	//   std::cout << "[b=" << col->birth_ <<"] ";
	//   if(col->paired_col_ != NULL) { std::cout << "[pc=" << col->paired_col_->lowest_idx_ <<"]     ";}
	//   else { std::cout << "[pc=" << "F" <<"]     ";}
	//   std::cout <<"| ";
	//   for(auto &cell : *(col->column_)) { std::cout << cell.key_ << " "; }
	//   std::cout << " |" << std::endl;
	// }


	// if(chains_in_F.size() > 3) {
	//   std::cout << "chains_in_F\n";
	//   for(auto chain : chains_in_F) {
	//     std::cout << *chain << std::endl;
	//   }
	//   std::cout << "-----end chains_in_F \n";
	// }

	//fp is the largest death index for <=d
	//Set col_fp: col_fp <- col_f1+...+col_fp (now in G); preserves lowest idx
	auto chain_fp = *(chains_in_F.begin()); //col_fp, with largest death <d index.

	for(auto other_col_it = chains_in_F.begin()+1;
	    other_col_it != chains_in_F.end(); ++other_col_it)
	{ plus_equal_column(chain_fp, chain_fp->column_, (*other_col_it)->column_); }

	// if(chains_in_F.size() > 3) {
	//   std::cout << "kernel = sum of all chains_in_F\n";
	//   std::cout << *chain_fp << std::endl;
	// }

	//Pair the col_fi, i = 1 ... p-1, according to the reflection diamond principle
	//Order the fi by reverse birth ordering <=_b           //true iff b(k1) > b(k2)
	auto cmp_birth = [this](Simplex_key k1, Simplex_key k2)->bool
	{ return birth_ordering_.reverse_birth_order(k1,k2); };

	// //i by >d value, contains at step i all b_j, j > i, and maybe b_i if not stolen
	// std::map< Simplex_key, matrix_chain * //available birth
	//         , decltype(cmp_birth) > available_birth_to_fidx(cmp_birth);
	// //for f1 to f_{p} (i by <=d), insertion in available_birth_to_fidx sorts by >=b
	// for(auto chain_f : chains_in_F)
	// { available_birth_to_fidx[ chain_f->birth() ] = chain_f; }

	//i by >d value, contains at step i all b_j, j > i, and maybe b_i if not stolen
	std::set< Simplex_key //available birth
		, decltype(cmp_birth) > available_birth(cmp_birth);
	//for f1 to f_{p} (i by <=d), insertion in available_birth_to_fidx sorts by >=b
	for(auto chain_f : chains_in_F) { available_birth.insert(chain_f->birth()); }

	// if(chains_in_F.size() > 1) {
	//   std::cout << "all birth by decreasing <b order\n";
	//   for(auto bb : available_birth) { std:: cout << bb << " ";}
	//   std::cout << std::endl;
	// }

	// available_birth_to_fidx[ chain_fp->birth() ] = chain_fp; //contains p elements, the birth of fp must be available to others
	//test if line above is necessary
	if(available_birth.find(chain_fp->birth()) == available_birth.end()) {
	    std::cout << "Miss chain_fp in available_birth when performing surjective_diamond \n"; }
	// if(available_birth_to_fidx.empty()) {std::cout << "available_birth_to_fidx empty.\n";}

	auto maxb_it = available_birth.begin();//max birth cycle
	auto maxb = *maxb_it; //max birth value, for persistence diagram
	available_birth.erase(maxb_it); //remove max birth cycle (stolen)


	// auto it_stop = chains_in_F.rend(); --it_stop;
	// for(auto chain_f_it  = chains_in_F.rbegin(); //by increasing death
	//         chain_f_it != it_stop; ++chain_f_it ) //chain_fp = begin() has max death

	for(auto chain_f_it  = chains_in_F.rbegin(); //by increasing death
	    *chain_f_it != chain_fp; ++chain_f_it ) //chain_fp = begin() has max death
	{ //find which reduced col has this birth

	    // std::cout << "   current before: " << *(*chain_f_it) << std::endl;

	    auto birth_it = available_birth.find((*chain_f_it)->birth());
	    if(birth_it == available_birth.end()) { //birth not available anymore

		// std::cout << "birth stolen\n";

		if (available_birth.empty()) {std::cout << "Should not be empty\n";}

		auto max_avail_b_it = available_birth.begin();
		Simplex_key max_avail_b = *max_avail_b_it;//max available birth
		//add all chains with smaller <d death and larger <b birth than max_avail_b
		for(auto chain_passed_it =  chains_in_F.rbegin();//all with smaller <d death
		    chain_passed_it != chain_f_it; ++chain_passed_it) {//but
		    if(cmp_birth((*chain_passed_it)->birth(), max_avail_b)) {//larger <b birth
			plus_equal_column( (*chain_f_it), (*chain_f_it)->column_,
					   (*chain_passed_it)->column_ );
			// std::cout << "X\n";
		    }
		}
		(*chain_f_it)->assign_birth(max_avail_b); //give new birth
		available_birth.erase(max_avail_b_it); //remove birth from availability
	    }
	    else { available_birth.erase(birth_it); } //birth not available anymore

	    // std::cout << "   current after:  " << *(*chain_f_it) << std::endl;
	}

	//Compute the new column zzsh + \sum col_h, for col_h in chains_in_H
	std::set< Simplex_key > col_bsh;
	std::pair< typename std::set< Simplex_key >::iterator, bool > res_insert;
	for(auto other_col : chains_in_H) { //Compute (\sum col_h) in a set
	    for(auto &cell : *(other_col->column_)) {
		res_insert = col_bsh.insert(cell.key_);
		if( !res_insert.second ) { col_bsh.erase(res_insert.first); } //1+1=0
	    }
	}
	//Create a new column with the new cycle value
	Column * new_col = new Column();
	Row    * new_row = new Row();
	matrix_.emplace_front( new_col
			       , new_row
			       , chain_fp //kills chain_fp == col_f1+...+col_fp from above
			       , -2       //belongs to H
			       // , 0
			       , cpx_->key(zzsh));
	auto chain_it = matrix_.begin();
	//Insert (\sum col_h) in matrix_
	for( Simplex_key idx : col_bsh ) //add all indices in col_new with canonical order
	{
	    Cell * new_cell = new Cell(idx, &(*chain_it));
	    new_col->push_back( *new_cell );  //add in column
	    lowidx_to_matidx_[idx]->row_->push_back( *new_cell ); //add in row
	}
	//New cell for zzsh
	Cell * new_cell = new Cell(cpx_->key(zzsh), &(*chain_it));
	new_col->push_back( *new_cell ); // zzsh has largest idx of all
	new_row->push_back( *new_cell );
	lowidx_to_matidx_[cpx_->key(zzsh)] = chain_it;
	chain_fp->assign_paired_col( &(*chain_it) );//pair chain_fp with the new chain
	chain_fp->assign_birth(-1); //now belongs to G now -> right interval [m-1,g]

	//Update persistence diagram with left interval [fil(b_max) ; fil(m))
	// // if(cpx_->filtration(max_birth) != cpx_->filtration(zzsh)) {
	// // if(max_birth_fil_ != prev_fil_) {
	//   persistence_diagram_.emplace_back( cpx_->dimension(zzsh)-1
	//                                    , max_birth_fil_//cpx_->filtration(max_birth)
	//                                    , curr_fil_);//cpx_->filtration(zzsh));
	// // }
	persistence_diagram_.emplace_back( cpx_->dimension(zzsh)-1
					   , maxb//cpx_->filtration(max_birth)
					   , cpx_->key(zzsh));//-1);//cpx_->filtration(zzsh));

    }


    //cpx_->key(zzsh) is the key of the simplex we remove, not a new one
    void backward_arrow( Simplex_handle zzsh )
    { //maintain the <=b order
	birth_ordering_.add_birth_backward(cpx_->key(zzsh));//num_arrow_);

	// std::cout << "backward_arrow \n";

	auto curr_col_it           = lowidx_to_matidx_.find(cpx_->key(zzsh));
	if(curr_col_it == lowidx_to_matidx_.end())
	{
	    std::cout << "Should not happen, backarrow. \n";
	    return; //case we didn't insert because IWD in max dim
	}

	Matrix_chain * curr_col    = &(*(curr_col_it->second));

	//Record all columns that get affected by the transpositions
	std::vector<Matrix_chain *> modified_columns;
	for(auto & hcell : *(curr_col->row_)) {
	    modified_columns.push_back(hcell.self_chain_);
	}
	//Sort by left-to-right order in the matrix_ (no order maintained in rows)
	std::stable_sort( modified_columns.begin(),modified_columns.end()
			  , [](Matrix_chain *mc1, Matrix_chain *mc2)
	{ return mc1->lowest_idx_ < mc2->lowest_idx_;} );

	// std::cout << "A\n";

	// std::cout << "      " << modified_columns.size() << "\n";

	//Modifies the pointer curr_col, not the other one.
	for(auto other_col_it = modified_columns.begin()+1;
	    other_col_it != modified_columns.end(); ++other_col_it)
	{ curr_col = arrow_transposition_case_study(curr_col, *other_col_it);
	    // display_mat();
	    // std::cout << "+++ done arr_transp.\n";
	}

	// std::cout << "B\n";


	//curr_col points to the column to remove by restriction of K to K-{\sigma}

	//to do cannot use cpx_->key(zzsh) for the birth here !
	if( curr_col->paired_col() == NULL ) { // in F
	    // // if(cpx_->filtration(zzsh) != curr_fil_) { //non-zero interval
	    //     persistence_diagram_.emplace_back( cpx_->dimension(zzsh)
	    //                                      , cpx_->filtration(zzsh)
	    //                                      , curr_fil_);
	    // // }
	    persistence_diagram_.emplace_back( cpx_->dimension(zzsh) - 1
					       , curr_col->birth()
					       , cpx_->key(zzsh));

	}
	else { //in H    -> paired c_g belongs to F now
	    curr_col->paired_col()->assign_paired_col(NULL);
	    curr_col->paired_col()->assign_birth(cpx_->key(zzsh));//num_arrow_); //closed interval WRONG
	    // curr_col->paired_col()->assign_fil(curr_fil_);
	}

	// std::cout << "C\n";


	if(curr_col->birth_ == -1) {  std::cout << "Error restrict a G column. \n"; } //cannot be in G

	//Erase the chain
	auto col_ptr = curr_col->column_;  auto row_ptr = curr_col->row_;
	for(auto c_it = col_ptr->begin(); c_it != col_ptr->end(); )
	{
	    auto tmp_it = c_it; ++c_it;
	    Cell * tmp_cell = &(*tmp_it);
	    tmp_it->base_hook_zzmat_h::unlink(); //rm from row
	    col_ptr->erase(tmp_it);
	    delete tmp_cell;
	}
	for(auto r_it = row_ptr->begin(); r_it != row_ptr->end(); )
	{
	    auto tmp_it = r_it; ++r_it;
	    Cell * tmp_cell = &(*tmp_it);
	    tmp_it->base_hook_zzmat_v::unlink(); //rm from col
	    row_ptr->erase(tmp_it);
	    delete tmp_cell;
	}

	// std::cout << curr_col->row_->size() << "  " << curr_col->column_->size() << "\n";


	// std::cout << "D\n";

	delete curr_col->row_;
	delete curr_col->column_;

	// std::cout << "DE\n";

	matrix_.erase(curr_col_it->second);
	lowidx_to_matidx_.erase(curr_col_it);


	// std::cout << "E\n";


	// std::cout << "F\n";
    }

    // void mini_function(int curr_col) {
    //     lowidx_to_matidx_[matrix_[curr_col].lowest_idx]   = curr_col;//update lowtomat
    //     lowidx_to_matidx_[matrix_[curr_col+1].lowest_idx] = curr_col+1;
    // }

    //   for(; curr_col < matrix_.size()-1; ++curr_col) //transpose curr_col with curr_col+1
    //   {
    //     flag = true;
    //     Simplex_key curr_low_idx = matrix_[curr_col].lowest_idx;
    //     for(auto & cell : matrix_[curr_col +1].column->col_)
    //     {
    //       if(cell.key_ == curr_low_idx) //M[i][i+1] != 0
    //       { arrow_transposition_case_study(curr_col);
    //         flag = false; break; }//end if M[i][i+1] !=0
    //     }
    //     //M[i][i+1] == 0
    //     if(flag) {
    //       std::swap(matrix_[curr_col], matrix_[curr_col+1]); //permute columns

    //       if(matrix_[curr_col].paired_col != -1)
    //         { matrix_[matrix_[curr_col].paired_col].paired_col = curr_col; }
    //       if(matrix_[curr_col+1].paired_col != -1)
    //         { matrix_[matrix_[curr_col+1].paired_col].paired_col = curr_col+1; }
    //     }
    //     mini_function(curr_col);
    //   }
    // }

    /* The two following methods
 * exchange members of the matrix_chains, except the column_ pointer. Modify
 * also the lowidx_to_matidx_ data structure, considering that the matrix chains
 * also exchange their lowest_idx_.
 *
 * Note that Cells in the matrix store a pointer to their matrix_chain_.
 * Consequently, exchanging columns would require to update all such pointers. That
 * is why we avoid doing it, and prefer exchanging all other attributes.
 */
    void exchange_lowest_indices_chains( Matrix_chain * curr_col
					 , Matrix_chain * other_col )
    { //exchange lowest_idx, update lowidx_to_matidx structure
	auto it_s = lowidx_to_matidx_.find(curr_col->lowest_idx_);
	auto it_t = lowidx_to_matidx_.find(other_col->lowest_idx_);
	std::swap(it_s->second, it_t->second);
	std::swap(curr_col->row_, other_col->row_);
	std::swap(curr_col->lowest_idx_, other_col->lowest_idx_);
    }
    void exchange_pairings_chains( Matrix_chain * curr_col
				   , Matrix_chain * other_col )
    { //exchange birth and pairing.
	std::swap(curr_col->birth_, other_col->birth_);
	std::swap(curr_col->paired_col_, other_col->paired_col_);
	// std::swap(curr_col->birth_fil_, curr_col->birth_fil_);
    }

    /*
 * Permute s and t, s goes up, whose insertions are adjacent. The bloc matrix gives:
 *                                 c_t                  c_t
 *                                  +                    +
 *      c_s c_t                    c_s c_s              c_s c_t
 *   s   1   1                  t   1   0            t   1   1
 *   t   0   1      --> either  s   0   1      or    s   0   1
 *
 */
    //return the new value of curr_col we continue with
    Matrix_chain * arrow_transposition_case_study( Matrix_chain * curr_col
						   , Matrix_chain * other_col )
    {
	// std::cout << "      arr_transp\n";
	switch( curr_col->birth() ) {
	case -2: {//c_s is in H
	    switch( other_col->birth() ) {
	    case -2: {//Case H x H, c_s+c_t paired w/ c_g+c_g', of death max<d{g,g'}=max
		// std::cout << "      H x H \n";
		auto curr_p_col  = curr_col->paired_col(); //c_s paired with c_g
		auto other_p_col = other_col->paired_col();//c_t paired with c_g'
		if( curr_p_col->lowest_idx_ < other_p_col->lowest_idx_) {//g<g', -->|c_s|
		    plus_equal_column( other_p_col, other_p_col->column_//c_g' <- c_g+c_g'
				       , curr_p_col->column_ );//of death g', low idx same, etc
		    plus_equal_column( other_col, other_col->column_
				       , curr_col->column_ ); //c_t <- c_t+c_s
		    return curr_col;//continue with c_s, paired with c_g of min death g
		}
		else {// g' < g, continue with --> |c_t|
		    plus_equal_column( curr_p_col, curr_p_col->column_
				       , other_p_col->column_ );//c_g <- c_g+c_g', death g, low idx same, etc
		    plus_equal_column( curr_col, curr_col->column_//c_s <- c_s+c_t, of low
				       , other_col->column_);//idx t now
		    //exchange lowest_idx, update lowidx_to_matidx structure
		    exchange_lowest_indices_chains(curr_col, other_col);
		    return other_col; //continue with c_t, paired with c_g' of min death g'
		}
	    }//end case HH
	    default: { //in H x (F U G) : c_s+c_t in H, paired with c_g
		// std::cout << "      H x (F U G) \n";
		// auto curr_p_col = curr_col->paired_col_; //c_s paired with c_g
		plus_equal_column(curr_col, curr_col->column_ //(still in H) <-> c_g
				  , other_col->column_);//c_s <- c_s+c_t
		//exchange lowest_idx, update lowidx_to_matidx structure
		exchange_lowest_indices_chains(curr_col, other_col);
		//exchange all members, except the column_ pointers, connected to cells
		// exchange_pairings_chains(curr_col, other_col);
		return other_col; //continue with c_t, still in F
	    }//end case H(F U G)
	    }//end switch
	}//end case H*
	case -1: {//c_s is in G
	    switch( other_col->birth() ) {
	    case -2: { //in G x H
		// std::cout << "      G x H \n";
		plus_equal_column( other_col, other_col->column_, curr_col->column_ );
		return curr_col;
	    }//end case GH
	    case -1: {//in G x G, c_s+c_t (in G) has death max<d{h,h'}=max{h,h'}
		// std::cout << "      G x G \n";
		auto curr_p_col = curr_col->paired_col_;  //c_s paired with c_h
		auto other_p_col = other_col->paired_col_;//c_t paired with c_h'
		if( curr_p_col->lowest_idx_ < other_p_col->lowest_idx_ ) {//h < h'
		    plus_equal_column( other_p_col, other_p_col->column_//c_h' <- c_h+c_h'
				       , curr_p_col->column_ );//of death h', low idx h'
		    plus_equal_column( other_col, other_col->column_//c_t <- c_s+c_t
				       , curr_col->column_ );//paired with c_h'
		    return curr_col;//continue with c_s, still paired with c_h
		}//endif
		else {//h' < h
		    plus_equal_column( curr_p_col, curr_p_col->column_//c_h <- c_h+c_h'
				       , other_p_col->column_ );//of death h, low idx h
		    plus_equal_column( curr_col, curr_col->column_//still paired with c_h
				       , other_col->column_);//c_s <- c_s+c_t, of low idx t
		    //exchange lowest_idx, update lowidx_to_matidx structure
		    exchange_lowest_indices_chains(curr_col, other_col);
		    return other_col; //continue with c_t, of min death h' and low idx s
		}//endelse
	    }//end case GG
	    default: { //in G x F, c_s+c_t (in F) has maximal birth
		// std::cout << "      G x F \n";
		plus_equal_column( other_col, other_col->column_
				   , curr_col->column_ );//c_t <- c_t+c_s, in F
		return curr_col;
	    }//end case GF
	    }//end switch
	}//end case G*
	default: {//c_s is in F
	    switch( other_col->birth() ) {
	    case -2: { // in F x H
		// std::cout << "      F x H \n";
		plus_equal_column( other_col, other_col->column_
				   , curr_col->column_ );//c_t <- c_s+c_t
		return curr_col;
	    }//end case FH
	    case -1: { //in F x G: c_s+c_t is in F, of max birth
		// std::cout << "      F x G \n";
		plus_equal_column(curr_col, curr_col->column_
				  , other_col->column_);//c_s <- c_s+c_t, still in F, low idx t
		//exchange lowest_idx, update lowidx_to_matidx structure
		exchange_lowest_indices_chains(curr_col, other_col);
		return other_col; //still in G
	    }//end case FG
	    default: { //in F x F: c_s+c_t has max<=b birth
		// std::cout << "      F x F \n";
		if(birth_ordering_.birth_order(curr_col->birth(), other_col->birth())) {
		    plus_equal_column( other_col, other_col->column_     //b_s <b b_t
				       , curr_col->column_ );//c_t <- c_s+c_t, of low idx t
		    return curr_col;
		}//endif
		else { //b_t <b b_s
		    plus_equal_column(curr_col, curr_col->column_
				      , other_col->column_);//c_s <- c_s+c_t, of low idx t

		    //exchange lowest_idx, update lowidx_to_matidx structure
		    exchange_lowest_indices_chains(curr_col, other_col);
		    return other_col;
		}//endelse
	    }//end case FF
	    }//end switch
	}//end case F*
	}//end switch
    }//end function



    //Class members
    Complex_ds                                           * cpx_; // complex
    std::map< Simplex_key //idx -> chain with lowest element at index idx in matrix_
    , typename std::list<Matrix_chain>::iterator > lowidx_to_matidx_;
    //arbitrary order for the matrix chains
    std::list< Matrix_chain >                              matrix_; // 0 ... m-1
    // birth_vector                                           birth_vector_; //<=b order
    birth_ordering                                         birth_ordering_;
    std::list< interval_t >                                persistence_diagram_;
    Simplex_key                                            num_arrow_; //current index
    //K_{i-1} -> K_i  prev_fil_ is the fil_val of K_{i-1}, curr_fil_ of K_i
    // Filtration_value                                       prev_fil_;
    // Filtration_value                                       curr_fil_;

    // filtration_values stores consecutive pairs (i,f) , (j,f') with f != f',
    // meaning that all inserted simplices with key in [i;j-1] have filtration value f
    std::vector< std::pair< Simplex_key, Filtration_value > > filtration_values_;


};

} //namespace zigzag_persistence

} //namespace Gudhi

#endif //_ZIGZAG_PERSISTENT_HOMOLOGY_

//----------------------------------------------------------------------------------









//FOR A POOL:
//In Zigzag_persistence creator:
//  , cell_pool_(new boost::object_pool< Cell > ()) //must be deleted
//As member of the class:
// boost::object_pool< Cell >                           * cell_pool_;
//In the code:
//        Cell * new_cell = new Cell(it2->key_, self1);    becomes,
//        Cell * new_cell = cell_pool_->construct(Cell(it2->key_, self1));
// and, for
//        Cell * tmp_ptr = &(*tmp_it); 
//        ...
//        delete tmp_ptr;                                  becomes,
//        cell_pool_->free(tmp_ptr);






//OLD arrow_transposition

//return the new value of curr_col
// matrix_chain * arrow_transposition_case_study( matrix_chain * curr_col
//                                              , matrix_chain * other_col )
// {
//   switch( curr_col->birth() ) 
//   {
//     case -2: 
//     { //in H
//       switch( other_col->birth() ) 
//       {
//         case -2: {//Case H x H, c_s+c_t paired w. c_g+c_g', of death max<d{g,g'}=max
//           auto curr_p_col  = curr_col->paired_col(); //c_s paired with c_g
//           auto other_p_col = other_col->paired_col();//c_t paired with c_g'
//           if( curr_p_col->lowest_idx_ < other_p_col->lowest_idx_) {//g<g', -->|c_s|
//             plus_equal_column( other_p_col, other_p_col->column_ 
//                                           , curr_p_col->column_ );//c_g' <- c_g+c_g'
//             plus_equal_column( other_col, other_col->column_ 
//                                         , curr_col->column_ ); //c_t <- c_t+c_s
//           }
//           else {// g' < g, continue with --> |c_t|
//             plus_equal_column( curr_p_col, curr_p_col->column_ 
//                              , other_p_col->column_ );//c_g <- c_g+c_g'
//             plus_equal_column( curr_col, curr_col->column_
//                                        , other_col->column_);//c_s <- c_s+c_t      
//             return other_col; //continue with c_t


//             // plus_equal_column( curr_p_col, curr_p_col->column_ 
//             //                              , other_p_col->column_ );//c_g <- c_g+c_g'
//             // col_swap( curr_col , other_col ); //exhange columns so as h' survives
//             // //change pairing
//             // std::swap( curr_p_col->paired_col_, other_p_col->paired_col_ );
//             // std::swap( curr_col->paired_col_, other_col->paired_col_ );
//             // plus_equal_column( other_col, other_col->column_ , curr_col->column_ ); //h' <- h+h'
//           }
//           break;
//         }




//         case -1: { //in H x G, (in H) c_s+c_t <-> c_g 
//           auto curr_p_col = curr_col->paired_col_; //c_s paired with c_g
//           plus_equal_column(curr_col, curr_col->column_ //(still in H) <-> c_g
//                                     , other_col->column_);//c_s <- c_s+c_t 
//           std::swap(curr_col, other_col); //continue with c_t, still in G      
//           break;


//           // auto curr_p_col = curr_col->paired_col_; 
//           // auto other_p_col = other_col->paired_col_;

//           // col_swap( curr_col , other_col );
//           // std::swap( curr_p_col->paired_col_, other_p_col->paired_col_ );
//           // std::swap( curr_col->paired_col_, other_col->paired_col_ );

//           // plus_equal_column( other_col, other_col->column_ , curr_col->column_ );

//           // other_col->assign_birth(-2); curr_col->assign_birth(-1);

//           // break;
//         }
//         default: 
//         { //in H x F
//           auto curr_p_col = curr_col->paired_col_; //c_s paired with c_g
//           plus_equal_column(curr_col, curr_col->column_ //(still in H) <-> c_g
//                                     , other_col->column_);//c_s <- c_s+c_t 
//           std::swap(curr_col, other_col); //continue with c_t, still in F      
//           break;

//           // col_swap( curr_col , other_col );
//           // plus_equal_column( other_col, other_col->column_ , curr_col->column_ );
//           // std::swap( curr_col->paired_col_, other_col->paired_col_ );
//           // other_col->paired_col_->paired_col_ = other_col;
//           // std::swap( curr_col->birth_ , other_col->birth_ );
//           // break;
//         }
//       }
//       break;
//     }
//     case -1: 
//     { //in G


//  //     std::cout << "Transposition with column in G \n \n \n ";

//       switch( other_col->birth() ) 
//       {
//         case -2: { //in G x H ok

//           // std::cout << "GH ";
//           // if(other_col->paired_col_ == curr_col) 
//           //   { std::cout << " Complain.......................................\n";}

//           plus_equal_column( other_col, other_col->column_ , curr_col->column_ );
//           break;
//         }
//         case -1: 
//         { //in G x G 

//           auto curr_p_col = curr_col->paired_col_; auto other_p_col = other_col->paired_col_;
//           if( curr_p_col->lowest_idx_ < other_p_col->lowest_idx_ ) 
//           {
//             plus_equal_column( other_p_col, other_p_col->column_ , curr_p_col->column_ ); //h' <- h+h'
//             plus_equal_column( other_col, other_col->column_ , curr_col->column_ ); //g' <- g+g'
//           }
//           else 
//           {
//             plus_equal_column( curr_p_col, curr_p_col->column_ , other_p_col->column_ ); 
//             col_swap( curr_col , other_col ); 
//             //change pairing
//             std::swap( curr_p_col->paired_col_, other_p_col->paired_col_ );
//             std::swap( curr_col->paired_col_, other_col->paired_col_ );

//             plus_equal_column( other_col, other_col->column_ , curr_col->column_ ); 
//           }
//           break;
//         }
//         default: 
//         { //in G x F
//           //---leave c1+c2 (belongs to F), continue with c1 in G
//           plus_equal_column( other_col, other_col->column_ , curr_col->column_ );
//           break;
//         }
//         break;
//       }
//     }   
//     default: 
//     { // in F
//       switch( other_col->birth() ) 
//       {
//         case -2: { // in F x H


//           plus_equal_column( other_col, other_col->column_ , curr_col->column_ );                  
//           break;
//         }
//         case -1: 
//         { //in F x G:
//           //---leave c1+c2 (belongs to F), continue with c2 in G
//           col_swap( curr_col , other_col ); 
//           //change pairing
//           std::swap( curr_col->paired_col_, other_col->paired_col_ );
//           curr_col->paired_col_->paired_col_ = curr_col;
//           std::swap( curr_col->birth_ , other_col->birth_ );
//           break;
//         }
//         default: 
//         { //in F x F: 
//           //---leave c1+c2 (has max birth), continue with c1 or c2 with min birth

//           if(birth_ordering_.birth_order(curr_col->birth(), other_col->birth()))
//        // if( birth_vector_[curr_col->birth()] < birth_vector_[other_col->birth()] )
//           {plus_equal_column( other_col, other_col->column_ , curr_col->column_ );}
//           else 
//           { 
//             col_swap( curr_col , other_col ); 
//             std::swap( curr_col->birth_ , other_col->birth_ );
//             plus_equal_column( other_col, other_col->column_ , curr_col->column_ );
//           }
//           break;
//         }
//         break;
//       }
//     }  
//   } //end switch
// }




































// //return the new value of curr_col
// matrix_chain * arrow_transposition_case_study( matrix_chain * curr_col, matrix_chain * other_col )
// {
//   switch( curr_col->birth_ ) 
//   {
//     case -2: 
//     { //in H
//       switch( other_col->birth_ ) 
//       {
//         case -2: 
//         { //in H x H
//           if( curr_col->paired_col_->lowest_idx_ < other_col->paired_col_->lowest_idx_) 
//           { //g < g'
//             plus_equal_column( other_col->paired_col_->column_
//                              , curr_col->paired_col_->column_); //g' <- g+g'
//             plus_equal_column( other_col->column_
//                              , curr_col_->column_ ); //h' <- h+h'


//             matrix_[matrix_[curr_col].paired_col].paired_col = curr_col+1;//new pairing gs
//             matrix_[matrix_[curr_col+1].paired_col].paired_col = curr_col;
//             std::swap(matrix_[curr_col], matrix_[curr_col+1]);          //permute columns
//             plus_equal_column(matrix_[curr_col].column, matrix_[curr_col+1].column); //h+h'
//           }
//           else {
//             plus_equal_column( matrix_[matrix_[curr_col].paired_col].column
//                              , matrix_[matrix_[curr_col+1].paired_col].column ); //g+g'
//             plus_equal_column( matrix_[curr_col].column
//                              , matrix_[curr_col+1].column );//h+h'
//             std::swap(matrix_[curr_col].lowest_idx, matrix_[curr_col+1].lowest_idx);
//             std::swap(matrix_[curr_col].row, matrix_[curr_col+1].row);
//           }
//           break;
//         }
//         case -1: { //in H x G
//           plus_equal_column( matrix_[curr_col].column
//                            , matrix_[curr_col+1].column );//g+h'
//           std::swap(matrix_[curr_col].lowest_idx, matrix_[curr_col+1].lowest_idx);
//           std::swap(matrix_[curr_col].row, matrix_[curr_col+1].row);
//           break;
//         }
//         default: { //in H x F
//           plus_equal_column( matrix_[curr_col].column
//                            , matrix_[curr_col+1].column );//f+h
//           std::swap(matrix_[curr_col].lowest_idx, matrix_[curr_col+1].lowest_idx);
//           std::swap(matrix_[curr_col].row, matrix_[curr_col+1].row);
//           break;
//         }
//       }
//       break;
//     }
//     case -1: 
//     { //in G
//       switch(matrix_[curr_col+1].birth) 
//       {
//         case -2: { //in G x H
//           matrix_[matrix_[curr_col].paired_col].paired_col = curr_col+1; //new pairing
//           matrix_[matrix_[curr_col+1].paired_col].paired_col = curr_col;
//           std::swap(matrix_[curr_col], matrix_[curr_col+1]);          //permute columns
//           plus_equal_column( matrix_[curr_col].column
//                            , matrix_[curr_col+1].column );//h'+g                  
//           break;
//         }
//         case -1: 
//         { //in G x G
//           if(matrix_[curr_col].paired_col < matrix_[curr_col+1].paired_col) {
//             plus_equal_column(  matrix_[matrix_[curr_col+1].paired_col].column
//                               , matrix_[matrix_[curr_col].paired_col].column ); //h+h'
//             matrix_[matrix_[curr_col].paired_col].paired_col = curr_col+1;//new pairing hs
//             matrix_[matrix_[curr_col+1].paired_col].paired_col = curr_col;
//             std::swap(matrix_[curr_col], matrix_[curr_col+1]);          //permute columns
//             plus_equal_column(matrix_[curr_col].column, matrix_[curr_col+1].column); //g+g'
//           }
//           else {
//             plus_equal_column( matrix_[matrix_[curr_col].paired_col].column
//                              , matrix_[matrix_[curr_col+1].paired_col].column ); //h+h'
//             plus_equal_column( matrix_[curr_col].column
//                              , matrix_[curr_col+1].column );//g+g'
//             std::swap(matrix_[curr_col].lowest_idx, matrix_[curr_col+1].lowest_idx);
//             std::swap(matrix_[curr_col].row, matrix_[curr_col+1].row);
//           }
//           break;
//         }
//         default: 
//         { //in G x F
//           matrix_[matrix_[curr_col].paired_col].paired_col = curr_col+1;//new pairing hs
//           std::swap(matrix_[curr_col], matrix_[curr_col+1]);          //permute columns
//           plus_equal_column( matrix_[curr_col].column
//                            , matrix_[curr_col+1].column );//f+g
//           break;
//         }
//         break;
//       }
//     }   
//     default: 
//     { // in F
//       switch(matrix_[curr_col+1].birth) 
//       {
//         case -2: { // in F x H
//           matrix_[matrix_[curr_col+1].paired_col].paired_col = curr_col; //new pairing
//           std::swap(matrix_[curr_col], matrix_[curr_col+1]);          //permute columns
//           plus_equal_column( matrix_[curr_col].column
//                            , matrix_[curr_col+1].column );//h+f                  
//           break;
//         }
//         case -1: { //in F x G
//           plus_equal_column( matrix_[curr_col].column
//                            , matrix_[curr_col+1].column );//f+g
//           std::swap(matrix_[curr_col].lowest_idx, matrix_[curr_col+1].lowest_idx);
//           std::swap(matrix_[curr_col].row, matrix_[curr_col+1].row);
//           break;
//         }
//         default: 
//         { //in F x F
//           if( birth_vector_[matrix_[curr_col].birth] < birth_vector_[matrix_[curr_col+1].birth] ) 
//           { std::swap(matrix_[curr_col], matrix_[curr_col+1]); }  //permute columns
//           else 
//             { std::swap(matrix_[curr_col].lowest_idx, matrix_[curr_col+1].lowest_idx); 
//               std::swap(matrix_[curr_col].row, matrix_[curr_col+1].row); 
//             }
//           plus_equal_column( matrix_[curr_col].column
//                            , matrix_[curr_col+1].column );//f+f'
//           break;
//         }
//         break;
//       }
//     }  
//   } //end switch
// }








//   std::map< Simplex_key, Column * >          mat_; // matrix

//   // to update properly
//   std::map< Simplex_key, int >               idx_to_colorder_; //idx_to_colorder_[simp_key] = order of col in mat
//   std::vector< Simplex_key >                 col_order_; //col_order_[ idx_to_colorder[key] ] = key
//   //------
//   std::map< Simplex_key, Simplex_key >       pairing_; // g -> h, f -> 0 and h -> -1
//   std::map< Simplex_key, Simplex_key >       low_to_col_; // idx -> index of reduced col with // lowest 1 at index idx
//   std::map< Simplex_key, Simplex_key >       col_to_low_; //inverse
//   std::map< Simplex_key, int >               colidx_to_birth_; //return the birth of a column from its idx
// //  std::map< Simplex_key, int >               colidx_to_death_;// same with death






// /*
// * check that everything is well updated ->colidx_to_birth, low_to_col etc  ?
// *
// */
// void forward_arrow( Simplex_handle zzsh )
// {
//   birth_vector_.add_birth_forward();

//   //store the boundary in a set
//   std::set< Simplex_key > col_bsh;
//   for( auto b_sh : cpx_->boundary_simplex_range(zzsh) ) 
//   {  col_bsh.insert( cpx_->key(b_sh) );  } //compute col for boundary of sigma, repr by a set

//   Simplex_key low_idx     = *(col_bsh.rbegin());   //idx of lowest element in col_bsh
//   Simplex_key col_low     =  low_to_col_[low_idx]; //idx of reduced col with low_idx as lowest index
//   Simplex_key paired_idx  =  pairing_[col_low];    //col with which col_low is paired
//   std::vector< Simplex_key > chains_in_H;          //for corresponding indices in H
//   std::vector< Simplex_key > chains_in_G;
//   std::pair< typename std::set< Simplex_key >::iterator, bool > res_insert;

//   while( true ) //reduce col_bsh with cycles of G
//   {
//     if( paired_idx > 0 ) //i.e. col_low \in G and paired_idx \in H
//     { 
//       chains_in_H.push_back(paired_idx); //remember the h with which g is paired
//       chains_in_G.push_back(col_low);    //remember the g
//       //add the column
//       for(auto &cell : mat_[col_low]->col_) {
//           res_insert = col_bsh.insert(cell.key_);
//           if( !res_insert.second ) { col_bsh.erase(res_insert.first); }
//       }

//       if(col_bsh.empty()) //col_bsh==0: [\partial sigma]=0 
//       {                    //sigma creates a cycle made with the col_h, h \in chains_in_H
//         for( auto idx_h : chains_in_H ) { //produce the sum of all col in chains_in_H
//           for(auto &cell : mat_[idx_h]->col_ ) {
//             res_insert = col_bsh.insert(cell.key_);
//             if( !res_insert.second ) { col_bsh.erase(res_insert.first); }
//           }
//         }
//         //Create a new column with the new cycle value
//         Column * new_col = new Column();
//         for( auto idx : col_bsh )  //add all indices in col_new with canonical order
//         {  
//           Cell * new_cell = new Cell(idx);
//           // todo      transverse_row
//           new_col->col_.push_back( *new_cell );  
//         }

//         Cell * new_cell = new Cell(cpx_->key(zzsh));
//         // todo      transverse_row

//         new_col->col_.push_back( *new_cell ); // and add sigma <- must have biggest idx of all

//         pairing_[cpx_->key(zzsh)] = 0;      //add index of sigma in F
//         mat_[cpx_->key(zzsh)]    = new_col; //add column to the matrix
//         colidx_to_birth_[cpx_->key(zzsh)] = cpx_->key(zzsh); //IWDP birth == key
//         col_to_low_[cpx_->key(zzsh)]      = cpx_->key(zzsh); //death == key
//         low_to_col_[cpx_->key(zzsh)]      = cpx_->key(zzsh);

//       //add row with cpx_->key(zzsh)
//       //......
//       /*******************************************************************/
//       /*******************************************************************/
//       /**************************** todo ********************************/
//       /*******************************************************************/
//       /*******************************************************************/

//         return;
//       }

//     //col_bsh != 0:
//       low_idx     = *(col_bsh.rbegin());   //idx of lowest element
//       col_low     =  low_to_col_[low_idx]; //reduced col with low_idx as lowest index
//       paired_idx  =  pairing_[col_low];    //col with which col_low is paired
//     }
//     else { //no more columns in G to reduce col_bsh, continue with F \cup G     [\partial sigma] != 0

//         std::vector< Simplex_key > chains_in_F;

//         while(true) {
//           if(paired_idx == 0) { chains_in_F.push_back(col_low); } //col_low \in F
//           else                { chains_in_H.push_back(paired_idx); } //paired_idx \in H

//           //add the col
//           for(auto &cell : mat_[col_low]->col_) {
//             res_insert = col_bsh.insert(cell.key_);
//             if( !res_insert.second ) { col_bsh.erase(res_insert.first); }
//           }

//           if(col_bsh.empty()) { break; }

//           low_idx     = *(col_bsh.rbegin());         //idx of lowest element
//           col_low     =  low_to_col_[low_idx];      //reduced col with low_idx as lowest index
//           paired_idx  =  pairing_[col_low];         //col with which col_low is paired
//         }

//     //Surjective weak diamond principle
//         // sort the fs by deaths (the closer from WD the bigger)
//         // the death is equal to the lowest index in the column
//         sort(chains_in_F.begin(),chains_in_F.end(), 
//             [this](Simplex_key k1, Simplex_key k2) {return col_to_low_[k1] < col_to_low_[k2];} );

//         //New value for col_fp: col_fp <- col_f1+...+col_fp \in G
//         Simplex_key fp = *(chains_in_F.rbegin()); //the one that gets cut
//         Column *col_fp = mat_[fp]; //corresponding column
//         chains_in_F.pop_back();    //remove fp from chains_in_F
//         for(auto c_idx : chains_in_F) {  plus_equal_column(col_fp,mat_[c_idx]); }

//         //New column col_m+1 that turns col_fp into a boundary
//         //col_m+1 <- sigma + col_h1 +...+ col_hi  with hi in chains_in_H. Here col_bsh is empty
//         for(auto c_idx : chains_in_H) {  
//             for(auto &cell : mat_[c_idx]->col_) {
//             res_insert = col_bsh.insert(cell.key_);
//             if( !res_insert.second ) { col_bsh.erase(res_insert.first); }
//           }
//         }
//         Column * new_col = new Column();
//         for( auto idx : col_bsh )  //add all indices in col_new 
//         {  
//           Cell * new_cell = new Cell(idx);
//           // todo      transverse_row
//           new_col->col_.push_back( *new_cell );  
//         }
//         Cell * new_cell = new Cell(cpx_->key(zzsh));
//         new_col->col_.push_back( *new_cell ); // and add sigma <- must have biggest idx of all
//                                               // todo new row
//         mat_[cpx_->key(zzsh)]     = new_col; //insert col for sigma in matrix
//         pairing_[cpx_->key(zzsh)] = -1;      //which belongs to H
//         pairing_[fp]              = cpx_->key(zzsh); //and is paired with fp now in G
//         col_to_low_[cpx_->key(zzsh)] = cpx_->key(zzsh);
//         low_to_col_[cpx_->key(zzsh)] = cpx_->key(zzsh);

//         //re-pairing of the fi
//         auto cmp_birth_vector = [this](Simplex_key k1, Simplex_key k2)->bool 
//                                     {return birth_vector_[k1] > birth_vector_[k2];};//ordered by reverse <=_b
//         std::map< Simplex_key, Simplex_key, decltype(cmp_birth_vector) > birth_to_idx(cmp_birth_vector); 

//         //for f1 to f_{p-1} sorted by <=_d
//         for(auto f_idx : chains_in_F) { birth_to_idx[ colidx_to_birth_[f_idx] ] = f_idx;} //inverse map
//         birth_to_idx[ colidx_to_birth_[fp] ] = fp; //contains p elements, the birth of fp must be available to others

//         auto bti_it = birth_to_idx.begin();  //points to max birth
//         Simplex_key max_birth = birth_to_idx.begin()->first;

//         Simplex_key curr_birth;
//         for(auto f_idx : chains_in_F) {
//           //colidx_to_birth_[f_idx] is the original birth of col_fidx
//           //find which reduced col has this birth
//           bti_it = birth_to_idx.find(colidx_to_birth_[f_idx]); 

//           if(bti_it->first == max_birth) { ++bti_it; } //if its max_birth, give next maximal birth
//           //while the col with curr_birth appears before in sorted chains_in_F
//           while( col_to_low_[bti_it->second] < col_to_low_[f_idx] ) 
//           {
//             plus_equal_column(mat_[f_idx],mat_[bti_it->second]);
//             ++bti_it;
//           }
//           bti_it->second = f_idx;
//         }
//     //update colidx_to_birth
//         bti_it = birth_to_idx.begin(); ++bti_it;
//         for(; bti_it != birth_to_idx.end(); ++bti_it) {
//           colidx_to_birth_[bti_it->second] = bti_it->first;
//         }
//         colidx_to_birth_.erase(fp); //fp not in F anymore

//     //update persistence diagram
//         persistence_diagram.emplace_back(max_birth,cpx_->key(zzsh)); //open interval
//         return;
//     }
//   }  
// }
