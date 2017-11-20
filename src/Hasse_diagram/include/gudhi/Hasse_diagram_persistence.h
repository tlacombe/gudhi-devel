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

#include <gudhi/Hasse_diagram_cell.h>
#include <gudhi/Hasse_diagram.h>

#ifndef HASSE_DIAGRAM_H
#define HASSE_DIAGRAM_H


namespace Gudhi {

namespace Hasse_diagram {
	
/**
 * This variable indicate if a warning should be given anytime a cell is 
 * deleted that have nondeleted cell in the coboundary.
**/ 

template < typename Cell_type >
class Hasse_diagram_persistence : public Hasse_diagram<Cell_type>
{
public:
	/**
	 * Default constructor.
	**/ 
	Hasse_diagram(){};
	
	/**
	 * Creating Hasse diagram from a file. The file format is the following:
	 * Number of cells
	 * cell dimension
	 * ids of cell boundary elements followed by the incidence coefficient.
	 * the two lines above are repeated for each cell. 
	 * It is assumed that the id of a cell is its position in the file. 	 
	**/ 	
    Hasse_diagram( const char* filename );
    
    /**
	 * Constructor to create a Hasse diagram from a vector of cells. It is assumed 
	 * that all the cells have boundaries set up. Setting up the coboundaries will
	 * be done in the constructor based on the information about boundaries. 
	**/ 
    Hasse_diagram( const std::vector< Cell_type* >& cells_ ):cells(cells_),number_of_deleted_cells(0)
    {
		this->set_up_coboundaries();
	};		

    
	//From here on we have implementation of methods that are required to use
	//this class with persistent homology engine.
	
	typedef typename Cell_type::Filtration_type Filtration_value;
    typedef unsigned Simplex_key;
    typedef Simplex_key Simplex_handle;
    
    size_t num_simplices() 
    {
		return this->cells.size();
    }
    
	
	Simplex_key key(Simplex_handle sh) 
	{
		return sh;
	}

	Simplex_key null_key() 
	{
		return std::numeric_limits<unsigned>::infinity();
	}

	Simplex_handle simplex(Simplex_key key) 
	{
		return key;
	}

	Simplex_handle null_simplex() 
	{
		return std::numeric_limits<unsigned>::infinity();
	}

	Filtration_value filtration(Simplex_handle sh) 
	{
		if (sh == null_simplex()) 
		{
		  return std::numeric_limits<Filtration_value>::infinity();
		}
		return this->cells[ sh ]->filtration;
	}

	int dimension(Simplex_handle sh) 
	{
		if (sh == null_simplex()) 
		{
		  return std::numeric_limits<Filtration_value>::infinity();
		}
		return this->cells[ sh ]->dimension;
	}

	int dimension() 
	{
		int top_dimension = 0;
		for ( size_t i = 0 ; i != this->cells.size() ; ++i )
		{
			if ( top_dimension < this->cells[i]->dimension )
			{
				top_dimension = this->cells[i]->dimension;
			}
		}
		return top_dimension;
	}

	std::pair<Simplex_handle, Simplex_handle> endpoints(Simplex_handle sh) 
	{
		return std::pair<Simplex_handle, Simplex_handle>(
							this->cells[sh].boundary[0].first->position, 
							this->cells[sh].boundary[1].first->position
														);
	}

	void assign_key(Simplex_handle sh, Simplex_key key) 
	{
		//TODO
		//complex_[sh].key_ = key;
	}

	//Boundary_simplex_range boundary_simplex_range(Simplex_handle sh) 
	//{
		//TODO
		//return Boundary_simplex_range(complex_[sh].boundary_.begin(), complex_[sh].boundary_.end());
	//}	
		
private:	
	std::vector< Cell_type* > cells;
	
	//to check how fragmented the data structure is (as a result of removing cells).
	size_t number_of_deleted_cells; 
	
	/**
	 * This procedure assumes that the boundaries are already set up for all
	 * the cells, and set up the coboundaries based on them.
	**/  
	void set_up_coboundaries();
	
	static double percentate_of_removed_cells_that_triggers_reorganization_of_structure;
};//Hasse_diagram

template <typename Cell_type>
double Hasse_diagram<Cell_type>::percentate_of_removed_cells_that_triggers_reorganization_of_structure = 0.8;


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
	//cell dimension
	//ids of cell boundary elements followed by the incidence coefficient.
	//Note that coboundary vector will be computed based on boundary vector.
	bool dbg = false;
	std::string line;
	
	this->number_of_deleted_cells = 0;
	
	std::ifstream in( filename );
	if ( !in.good() )
	{
		std::cout << "The file do not exist, program will now terminate.\n";
		throw "The file do not exist, program will now terminate.\n";
	}
	
	std::getline(in, line);
	while ( line[0] == '#' )
	{
		std::getline(in, line);
	}
	std::stringstream iss(line);
	
	
	unsigned number_of_cells;
	iss >> number_of_cells;
	this->cells.reserve( number_of_cells );	
	
	std::getline(in, line);
	
	if ( dbg )
	{
		std::cout << "Number of cells : " << number_of_cells << std::endl;		
	}
	
	size_t size_of_last_boundary = 10;//to initially reserve a vector for bounary elements.
	for ( size_t i = 0 ; i != number_of_cells ; ++i )
	{
		Cell_type* new_cell = new Cell_type();		
		while ( line[0] == '#' )
		{
			std::getline(in, line);	
		}
		
		iss.str("");
		iss.clear();
		iss << line;		
		iss >> new_cell->position >> new_cell->dimension;
		
		if ( dbg )std::cout << "Position and dimension of the cell : " << new_cell->position << " , " << new_cell->dimension << std::endl;
		
		if ( new_cell->position != i )
		{
			std::cerr << "Wrong numeration of cells in the file. Cell number : " << i << " is marked as : " << new_cell->position << " in the file." << std::endl;
			throw "Wrong numeration of cells in the file.";
		}
		if ( iss.good() )
		{
			//in this case we still have a filtration value to be read
			//from the file.
			iss >> new_cell->filtration;
			if ( dbg )std::cout << "Filtration of the cell : " << new_cell->filtration << std::endl;
		}
		else
		{
			new_cell->filtration = 0;
		}
		
		std::getline(in, line);
		while ( line[0] == '#' )
		{
			std::getline(in, line);
		}		

		iss.str("");
		iss.clear();
		iss << line;
		unsigned cell_id;
		int incidence_coef;
		std::vector< std::pair< unsigned,int > > bdry;
		bdry.reserve( size_of_last_boundary );
		while ( iss.good() )
		{			
			iss >> cell_id;
			if ( !iss.good() )continue;
						
			if ( cell_id >= i )
			{
				std::cerr << "Wrong format of a file. THe cell number : " << i << " contain in a boundary a cell that has not been introduced yet.\n";
			}
			iss >> incidence_coef;			
			if ( dbg )std::cout << "( " <<  cell_id << " , " << incidence_coef << " ), ";
			bdry.push_back( std::pair< unsigned,int >(cell_id,incidence_coef) );
		}
		size_of_last_boundary = bdry.size();
		new_cell->boundary.reserve( size_of_last_boundary );
		for ( size_t bd = 0 ; bd != size_of_last_boundary ; ++bd )
		{
			new_cell->boundary.push_back( std::make_pair(this->cells[ bdry[bd].first ] , bdry[bd].second) );
		}		
		this->cells.push_back( new_cell );
		if ( dbg )
		{
			std::cout << "Done with this cell. \n";
			getchar();
		}
		
		std::getline(in, line);		
		while ( line[0] == '#' )
		{
			std::getline(in, line);			
		}		
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
			this->cells[ this->cells[i]->boundary[bd].first->position ]
			->
			coBoundary.push_back
			                    ( 
			                    std::make_pair(this->cells[i], 
								this->cells[i]->boundary[bd].second)
			                    );
		}
	}
}

template < typename Cell_type >
void Hasse_diagram<Cell_type>::write_to_file( const char* filename )
{
	std::ofstream out( filename );
	//If there are any deleted cells, then we need to clean up the structure 
	//first before writing it to a file. The reason for that is because the 
	//file format assumes a continuous enumeration of cells (from zero to the
	//number od cells). This is not satisfied if we have in the structure any
	//deleted elements.	
	if ( this->number_of_deleted_cells != 0 )
	{
		this->clean_up_the_structure();
	}
	//now we can write to a file the number of (non deleted) cells in the structure. 
	out << this->cells.size() << std::endl;
	//and then the rest of the Hasse diagram.
	out << *this;		
	out.close();
}//template < typename Cell_type >


}//namespace Hasse_diagram
}//namespace Gudhi

#endif //HASSE_DIAGRAM_H
