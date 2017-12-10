 /*    This file is part of the Gudhi Library. The Gudhi library 
  *    (Geometric Understanding in Higher Dimensions) is a generic C++ 
  *    library for computational topology.
  *
  *    Author(s):       Pawel Dlotko
  *
  *    Copyright (C) 2017 Swansea University, UK
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

/** \brief The concept Cubical_complex_type describes the requirements 
  * for a type to use in Topological_inference_with_cubical_complexes 
  * class.
  */
struct Cubical_complex_type
{
/**
* A constructor creating a grid of sizes specified in resolution_of_a_grid_.
**/ 
Cubical_complex_type(const std::vector< unsigned >& resolution_of_a_grid_);

/**
* A constructor creating a grid of sizes specified in resolution_of_a_grid_.
* In addition the grid is periodic in the directions specified in directions_in_which_periodic_b_cond_are_to_be_imposed.
**/ 
Cubical_complex_type(const std::vector< unsigned >&  resolution_of_a_grid_,const std::vector<bool>& directions_in_which_periodic_b_cond_are_to_be_imposed);

/**
* After setting up the filtration values at all top dimensional cubes, this function is assumed
* to set it up at all lower dimensional cubes so that the resulting filtration is a lower star
* filtration.
**/ 
void impose_lower_star_filtration();

/**
* In some cases some work need to be done in the cubical complex data structure prior to 
* computing persistent homolgy of it. This function is supposed to do that work. From the user
* point of view, after the function is called on the data structure, any persistent homology 
* engine provided by Gudhi should correctly work on the given data structure. 
**/ 
void initialize_arrays_for_persistence_computation();

/**
* A procedure that store the top filtration values on top dimensional cubes in the
* cubical complex in the Perseus style.
**/ 
void store_in_perseus_format( const char* filename );

/**
* A type of iterator that iterates through top dimensional cubes of the complex.
**/ 
typename Top_dimensional_cells_iterator;

/**
* Begin iterator for Top_dimensional_cells_iterator.
**/ 
Top_dimensional_cells_iterator top_dimensional_cells_iterator_begin();

/**
* End iterator for Top_dimensional_cells_iterator.
**/ 
Top_dimensional_cells_iteratortop_dimensional_cells_iterator_end();


};
