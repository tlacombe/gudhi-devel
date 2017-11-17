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



#ifndef HASSE_DIAGRAM_H
#define HASSE_DIAGRAM_H


namespace Gudhi {

namespace Hasse_diagram {
	
template < typename Cell_type > class Hasse_diagram;

template <typename Incidence_type, typename Additional_information = void>
class Cell
{
public:
    /**
     * TODO
    **/  
	Cell():dimension(0),position(0),deleted_(false){}
	
	/**
     * TODO
    **/  
	Cell( unsigned dim ):dimension(dim),position(0),deleted_(false){}
	
	/**
     * TODO
    **/  
	Cell( const std::vector< std::pair<Cell*,int> >& boundary_ , unsigned dim ):
	dimension(dim),boundary(boundary_),position(0),deleted_(false){}
	
	/**
     * TODO
    **/  
	Cell( const std::vector< std::pair<Cell*,int> >& boundary_ , const std::vector< std::pair<Cell*,int> >& coboundary_,
		 unsigned dim ):dimension(dim),boundary(boundary_),coBoundary(coboundary_),
		 position(0),deleted_(false){}
	
	/**
     * TODO
    **/  
	Cell( const std::vector< std::pair<Cell*,int> >& boundary_ , const std::vector< std::pair<Cell*,int> >& coboundary_,
	Additional_information ai, unsigned dim ):
	dimension(dim),boundary(boundary_),coBoundary(coboundary_),additional_info(ai),
	position(0),deleted_(false){}
	
	/**
     * TODO
    **/  
	Cell(Additional_information ai, unsigned dim ):
	dimension(dim),additional_info(ai),position(0),deleted_(false){}
	
	/**
     * TODO
    **/  
	inline std::vector< std::pair<Cell*,int> >& get_boundary(){return this->boundary;}
	
	/**
     * TODO
    **/  
	inline std::vector< std::pair<Cell*,int> >& get_coBoundary(){return this->coBoundary;}
	
	/**
     * TODO
    **/  
	inline unsigned& get_dimension(){return this->dimension;}
	
	/**
     * TODO
    **/  
	inline Additional_information& get_additional_information(){return this->additional_info;}
	
	/**
	 * TODO
	**/ 
	inline size_t& get_position(){return this->position;}
	
	/**
	 * TODO
	**/ 	
	bool& deleted(){return this->deleted_;}
	
	template < typename Cell_type >
	friend class Hasse_diagram;
	
	/**
	 * Procedure to remove deleted boundary and coboundary elements from the 
	 * vectors of boundary and coboundary elements of this cell.
	**/ 
	void remove_deleted_elements_from_boundary_and_coboundary()
	{
		std::vector< std::pair<Cell*,int> > new_boundary;
		new_boundary.reserve( this->boundary.size() );
		for ( size_t bd = 0 ; bd != this->boundary.size() ; ++bd )
		{
			if ( !this->boundary[bd]->deleted() )
			{
				new_boundary.push_back( this->boundary[bd] );
			}
		}
		this->boundary.swap( new_boundary );
		
		std::vector< std::pair<Cell*,int> > new_coBoundary;
		new_coBoundary.reserve( this->coBoundary.size() );
		for ( size_t cbd = 0 ; cbd != this->coBoundary.size() ; ++cbd )
		{
			if ( !this->coBoundary[cbd]->deleted() )
			{
				new_coBoundary.push_back( this->coBoundary[cbd] );
			}
		}
		this->coBoundary.swap( new_coBoundary );
	}
private:	
	std::vector< std::pair<Cell*,int> > boundary;
	std::vector< std::pair<Cell*,int> > coBoundary;
	unsigned dimension;	
	Additional_information additional_info;	
	size_t position;
	bool deleted_; //for lazy delete	
};//Cell


template < typename Cell_type >
class Hasse_diagram
{
public:
	/**
	 * TODO
	**/ 
	Hasse_diagram(){};
	
	/**
	 * TODO
	**/ 	
    Hasse_diagram( const char* filename );
    
    /**
	 * TODO
	**/ 
    Hasse_diagram( const std::vector< Cell_type* >& cells_ ):cells(cells_),number_of_deleted_cells(0){};


	//those iterators should skip the deleted cells. 
	/**
	 * TODO
	**/ 
	typedef typename std::vector< size_t >::const_iterator All_cells_iterator;
    typedef typename std::vector< size_t > All_cells_iterator_range;

    All_cells_iterator_range all_cells_range() 
    {
        return this->cells;
    }
    
	/**
	 * After many operation of deleting cells, this->cells vector may became
	 * very fragmented. Also, the complexity of operation using all the iterators
	 * depends on the actual size of a structure (where the deleted elements are still
	 * stored. This procedure remove permanently all the deleted elements. Ideally, 
	 * it should be initialized when the percentage of deleted elements is larger than a 
	 * predefined constant. 
	**/     
    void clean_up_the_structure()
    {
		//count the number of not deleted cells:
		unsigned number_of_non_deleted_cells = 0;
		for ( size_t i = 0 ; i != this->cells.size() ; ++i )
		{
			if ( !this->cells[i]->deleted() )++number_of_non_deleted_cells;
		}
		//create a new vector to store the undeleted cells:
		std::vector< Cell_type* > new_cells;
		new_cells.reserve( number_of_non_deleted_cells );
		//fill the new vector in and adjust the new positions.
		//In the same time make sure that the boundary and coboundary vectors 
		//in every cell are valid.
		size_t counter = 0;
		for ( size_t i = 0 ; i != this->cells.size() ; ++i )
		{
			if ( !this->cells[i]->deleted() )
			{
				new_cells.push_back( this->cells[i] );
				this->cells[i]->position = counter;
				this->cells[i]->remove_deleted_elements_from_boundary_and_coboundary();
				++counter;
			}
			else
			{
				delete this->cells[i];
			}
		}
		this->cells.swap(new_cells);
	} 
private:	
	std::vector< Cell_type* > cells;
	
	size_t number_of_deleted_cells; //to check how fragmented the data structure is (as a result of removing cells). Can we even use it if deleting is happening on the level of cells??
	
	/**
	 * This procedure assumes that the boundaries are already set up for all
	 * the cells, and set up the coboundaries based on them.
	**/  
	void set_up_coboundaries();
};//Hasse_diagram


template < typename Cell_type >
Hasse_diagram<Cell_type>::Hasse_diagram( const char* filename )
{
	//We assume that the cells in the file are enumerated in increasing order. 
	//The idea of the i-th cell in the file is by default i (starting from zero). 
	//Moreover, the cells are provided from low dimensional to high dimensiona.
	//By doing so, we know that when constructing a given cell, all its boundary
	//has already been constructed. 
	//Here is the format of a file:
	//Number of cells
	//cell id cell dimension
	//ids of cell boundary elements follwoed by the incdence coefficient.
	//Note that coboundary vector will be computed based on boundary vector.
	bool dbg = true;
	
	this->number_of_deleted_cells = 0;
	
	std::ifstream in( filename );
	if ( !in.good() )
	{
		std::cout << "The file do not exist, program will now terminate.\n";
		throw "The file do not exist, program will now terminate.\n";
	}
	unsigned number_of_cells;
	in >> number_of_cells;
	this->cells.reserve( number_of_cells );	
	
	size_t size_of_last_boundary = 10;//to initially reserve a vector for bounary elements.
	for ( size_t i = 0 ; i != number_of_cells ; ++i )
	{
		Cell_type* new_cell = new Cell_type();
		new_cell->position = i;
		in >> new_cell->dimension;
		std::string line;
		std::getline(in, line);
		std::istringstream iss(line);
		unsigned cell_id;
		int incidence_coef;
		std::vector< std::pair< unsigned,int > > bdry;
		bdry.reserve( size_of_last_boundary );
		while ( iss.good() )
		{
			iss >> cell_id;
			if ( cell_id >= i )
			{
				std::cerr << "Wrong format of a file. THe cell number : " << i << " contain in a boundary a cell that has not been introduced yet.\n";
			}
			iss >> incidence_coef;			
			bdry.push_back( std::pair< unsigned,int >(cell_id,incidence_coef) );
		}
		size_of_last_boundary = bdry.size();
		new_cell->boundary.reserve( size_of_last_boundary );
		for ( size_t bd = 0 ; bd != size_of_last_boundary ; ++bd )
		{
			new_cell->boundary.push_back( this->cells[ bdry[bd].first ] , bdry[bd].second );
		}		
		this->cells.push_back( new_cell );
	}
	//now once the boundaries are set, we are to set up the coboundaries.
	
	this->set_up_coboundaries();
}

template < typename Cell_type >
void Hasse_diagram<Cell_type>::set_up_coboundaries()
{	
	//first we check the number of coboundary elements for each cell:
	size_t number_of_cells = this->cells.size();
	std::vector< unsigned > sizes_of_coboundary( number_of_cells , 0 );
	for ( size_t i = 0 ; i != number_of_cells ; ++i )
	{
		std::vector< std::pair<Cell_type*,int> > bdry = this->cells[i]->get_boundary();
		for ( size_t bd = 0 ; bd != bdry.size() ; ++bd )
		{
			sizes_of_coboundary[ bdry[bd].first->get_position() ]++;
		}
	}
	
	//now we reserve the space for all coboundaries
	for ( size_t i = 0 ; i != number_of_cells ; ++i )
	{
		this->cells[i]->coBoundary.reserve( sizes_of_coboundary[i] );
	}
	
	//and now we set up the coboundaries.
	for ( size_t i = 0 ; i != number_of_cells ; ++i )
	{
		for ( size_t bd = 0 ; bd != this->cells[i]->boundary.size() ; ++bd )
		{
			this->cells[i]->boundary[bd]->coBoundary.push_back( this->cells[i] );
		}
	}
}


}//namespace Hasse_diagram
}//namespace Gudhi

#endif //HASSE_DIAGRAM_H
