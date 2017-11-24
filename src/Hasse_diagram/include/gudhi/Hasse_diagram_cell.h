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

#include <iostream>
#include <vector>
#include <algorithm>
#include <string>
#include <sstream>
#include <type_traits>
#include <cstdlib>
#include <cstdio>




#ifndef HASSE_DIAGRAM_CELL_H
#define HASSE_DIAGRAM_CELL_H


namespace Gudhi {

namespace Hasse_diagram {


template < typename Cell_type > class Hasse_diagram;
template <typename Cell_type> class is_before_in_filtration;

template <typename Incidence_type_, typename Filtration_type_ , typename Additional_information_ = void>
class Hasse_diagram_cell
{
public:
	typedef Incidence_type_ Incidence_type;
	typedef Filtration_type_ Filtration_type;
	typedef Additional_information_ Additional_information;

    /**
     * Default constructor.
    **/
	Hasse_diagram_cell():dimension(0),position(0),deleted_(false){}

	/**
     * Constructor of a cell of dimension dim.
    **/
	Hasse_diagram_cell( unsigned dim ):dimension(dim),position(0),deleted_(false){}
	
	/**
     * Constructor of a cell of dimension dim.
    **/
	Hasse_diagram_cell( unsigned dim , Filtration_type filt_ ):dimension(dim),position(0),deleted_(false),filtration(filt_){}

	/**
     * Constructor of a cell of dimension dim with a given boundary.
    **/
	Hasse_diagram_cell( const std::vector< std::pair<Hasse_diagram_cell*,int> >& boundary_ , unsigned dim ):
	dimension(dim),boundary(boundary_),position(0),deleted_(false){}

	/**
     * Constructor of a cell of dimension dim with a given boundary and coboundary.
    **/
	Hasse_diagram_cell( const std::vector< std::pair<Hasse_diagram_cell*,int> >& boundary_ , const std::vector< std::pair<Hasse_diagram_cell*,int> >& coboundary_,
		 unsigned dim ):dimension(dim),boundary(boundary_),coBoundary(coboundary_),
		 position(0),deleted_(false){}

	/**
     * Constructor of a cell of dimension dim with a given boundary, coboundary and
     * additional information.
    **/
	Hasse_diagram_cell( const std::vector< std::pair<Hasse_diagram_cell*,int> >& boundary_ , const std::vector< std::pair<Hasse_diagram_cell*,int> >& coboundary_,
	Additional_information ai, unsigned dim ):
	dimension(dim),boundary(boundary_),coBoundary(coboundary_),additional_info(ai),
	position(0),deleted_(false){}

	/**
     * Construcor of a cell of dimension dim having given additional information.
    **/
	Hasse_diagram_cell(Additional_information ai, unsigned dim ):
	dimension(dim),additional_info(ai),position(0),deleted_(false){}

	/**
     * Procedure to get the boundary of a fiven cell. The output format
     * is a vector of pairs of pointers to boundary elements and incidence
     * coefficients.
    **/
	inline std::vector< std::pair<Hasse_diagram_cell*,int> >& get_boundary(){return this->boundary;}

	/**
     * Procedure to get the coboundary of a fiven cell. The output format
     * is a vector of pairs of pointers to coboundary elements and incidence
     * coefficients.
    **/
	inline std::vector< std::pair<Hasse_diagram_cell*,int> >& get_coBoundary(){return this->coBoundary;}

	/**
     * Procedure to get the dimension of a cell.
    **/
	inline unsigned& get_dimension(){return this->dimension;}

	/**
     * Procedure to get additional information about the cell.s
    **/
	inline Additional_information& get_additional_information(){return this->additional_info;}

	/**
	 * Procedure to retrive position of the cell in the structure. It is used in
	 * the implementation of Hasse diagram and set by it. Note that removal of
	 * cell and subsequent call of clean_up_the_structure will change those
	 * positions.
	**/
	inline size_t& get_position(){return this->position;}
	
	/**
	 * Accessing the filtration of the cell.
	**/
	inline Filtration_type& get_filtration(){return this->filtration;}

	/**
	 * A procedure used to check if the cell is deleted. It is used by the
	 * subsequent implementation of Hasse diagram that is absed on lazy
	 * delete.
	**/
	inline bool deleted(){return this->deleted_;}

	template < typename Cell_type >
	friend class Hasse_diagram;
	
	template < typename Cell_type >
	friend class is_before_in_filtration;
	
	
	template <typename Complex_type , typename Cell_type >  
	friend std::vector<Cell_type*> convert_to_vector_of_Cell_type( Complex_type& cmplx );

	/**
	 * Procedure to remove deleted boundary and coboundary elements from the
	 * vectors of boundary and coboundary elements of this cell.
	**/
	void remove_deleted_elements_from_boundary_and_coboundary()
	{
		std::vector< std::pair<Hasse_diagram_cell*,int> > new_boundary;
		new_boundary.reserve( this->boundary.size() );
		for ( size_t bd = 0 ; bd != this->boundary.size() ; ++bd )
		{
			if ( !this->boundary[bd].first->deleted() )
			{
				new_boundary.push_back( this->boundary[bd] );
			}
		}
		this->boundary.swap( new_boundary );

		std::vector< std::pair<Hasse_diagram_cell*,int> > new_coBoundary;
		new_coBoundary.reserve( this->coBoundary.size() );
		for ( size_t cbd = 0 ; cbd != this->coBoundary.size() ; ++cbd )
		{
			if ( !this->coBoundary[cbd].first->deleted() )
			{
				new_coBoundary.push_back( this->coBoundary[cbd] );
			}
		}
		this->coBoundary.swap( new_coBoundary );
	}

	/**
	 * Writing to a stream operator.
	**/
	friend std::ostream& operator<<( std::ostream& out, const Hasse_diagram_cell<Incidence_type,Filtration_type,Additional_information>& c )
	{
		 //cout << "position : " << c.position << ", dimension : " << c.dimension << ", filtration: " << c.filtration << ", size of boudary : " <<  c.boundary.size() << "\n";
		 out << c.position << " " << c.dimension << " " << c.filtration << std::endl;
		 for ( size_t bd = 0 ; bd != c.boundary.size() ; ++bd )
	     {
			 //do not write out the cells that has been deleted
			 if ( c.boundary[bd].first->deleted() )continue;
			 out <<  c.boundary[bd].first->position << " " << c.boundary[bd].second << " ";
		 }		
		 out << std::endl;
		return out;
	}
	
		
	/**
	 * Procedure that return vector of pointers to boundary elements of a given cell.
	**/  	
	inline std::vector< Hasse_diagram_cell* > get_list_of_boundary_elements()
	{
		std::vector< Hasse_diagram_cell* > result;	
		size_t size_of_boundary = this->boundary.size();
		result.reserve( size_of_boundary );	
		for ( size_t bd = 0 ; bd != size_of_boundary ; ++bd )
		{
			result.push_back( this->boundary[bd].first );
		}		
		return result;
	}
	
	/**
	 * Procedure that return vector of positios of boundary elements of a given cell.
	**/  	
	inline std::vector< unsigned > get_list_of_positions_of_boundary_elements()
	{
		std::vector< unsigned > result;	
		size_t size_of_boundary = this->boundary.size();
		result.reserve( size_of_boundary );	
		for ( size_t bd = 0 ; bd != size_of_boundary ; ++bd )
		{
			result.push_back( this->boundary[bd].first->position );
		}		
		return result;
	}
	
	/**
	 * Function that display a string being a signature of a structure. 
	 * Used mainly for debugging purposes. 
	**/ 
	std::string full_signature_of_the_structure()
	{
		std::string result;
		result += "dimension: ";
		result += std::to_string(this->dimension);
		result += " filtration: ";
		result += std::to_string(this->filtration);
		result += " position: ";
		result += std::to_string(this->position);
		result += " deleted_: ";
		result += std::to_string(this->deleted_);
		
		//if the Additional_information is not void, add them to
		//the signature as well.
		if ( std::is_same<Additional_information, void>::value )
		{			
			result += " Additional_information: ";
			result += std::to_string(this->additional_info);
		}
		result += " boundary ";
		for ( size_t bd = 0 ; bd != this->boundary.size() ; ++bd )
		{
			result += "( " + std::to_string(this->boundary[bd].first->position);
			result += " " + std::to_string(this->boundary[bd].second);
			result += ") ";			
		}
		
		result += " coBoundary ";
		for ( size_t cbd = 0 ; cbd != this->coBoundary.size() ; ++cbd )
		{
			result += "( " + std::to_string(this->coBoundary[cbd].first->position);
			result += " " + std::to_string(this->coBoundary[cbd].second);
			result += ") ";			
		}
		
		return result;
	}
		
	
protected:
	std::vector< std::pair<Hasse_diagram_cell*,int> > boundary;
	std::vector< std::pair<Hasse_diagram_cell*,int> > coBoundary;
	unsigned dimension;
	Additional_information additional_info;
	size_t position;
	bool deleted_;
	Filtration_type filtration;

	/**
	 * A procedure to delete a cell. It is a private function of the Hasse_diagram_cell
	 * class, since in the Hasse_diagram class I want to have a control
	 * of removal of cells. Therefore, to remove cell please use
	 * remove_cell in the Hasse_diagram structure.
	**/
	void delete_cell(){ this->deleted_ = true; }
};//Hasse_diagram_cell



}//namespace Hasse_diagram
}//namespace Gudhi

#endif //CELL_H
