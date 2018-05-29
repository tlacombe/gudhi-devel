/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Pawel Dlotko
 *
 *    Copyright (C) 2017  Swansea University
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

#ifndef PERSISTENCE_LANDSCAPE_ON_GRID_INTERFACE_H_
#define PERSISTENCE_LANDSCAPE_ON_GRID_INTERFACE_H_

#include <gudhi/Persistence_landscape_on_grid.h>

namespace Gudhi {
namespace Persistence_representations {


class Persistence_landscape_on_grid_interface : public Persistence_landscape_on_grid 
{
 public:
  Persistence_landscape_on_grid_interface():Persistence_landscape_on_grid(){}

  Persistence_landscape_on_grid_interface(const std::vector<std::pair<double, double> >& p, double grid_min_, double grid_max_,
  size_t number_of_points_):Persistence_landscape_on_grid(p, (grid_min_ == grid_max_ ? std::numeric_limits<double>::max() : grid_min_), 
                                                          (grid_min_ == grid_max_ ? std::numeric_limits<double>::max() : grid_max_),number_of_points_){}
                                

  Persistence_landscape_on_grid_interface(const std::vector<std::pair<double, double> >& p, double grid_min_, double grid_max_,
  size_t number_of_points_, unsigned number_of_levels_of_landscape):
  Persistence_landscape_on_grid(p, (grid_min_ == grid_max_ ? std::numeric_limits<double>::max() : grid_min_), 
  (grid_min_ == grid_max_ ? std::numeric_limits<double>::max() : grid_max_),number_of_points_, number_of_levels_of_landscape){}


  Persistence_landscape_on_grid_interface(const char* filename, double grid_min_, double grid_max_, size_t number_of_points_,
  unsigned number_of_levels_of_landscape, int dimension_ = -1):
  Persistence_landscape_on_grid(filename, (grid_min_ == grid_max_ ? std::numeric_limits<double>::max() : grid_min_), 
  (grid_min_ == grid_max_ ? std::numeric_limits<double>::max() : grid_max_), number_of_points_, number_of_levels_of_landscape, dimension_ ){}
  

  Persistence_landscape_on_grid_interface(const char* filename, double grid_min_, double grid_max_, size_t number_of_points_,
  int dimension_ = -1):Persistence_landscape_on_grid(filename,(grid_min_ == grid_max_ ? std::numeric_limits<double>::max() : grid_min_),
  (grid_min_ == grid_max_ ? std::numeric_limits<double>::max() : grid_max_),number_of_points_,dimension_ ){}


  Persistence_landscape_on_grid_interface(const char* filename, size_t number_of_points, unsigned number_of_levels_of_landscape, int dimension = -1):
  Persistence_landscape_on_grid(filename,number_of_points,number_of_levels_of_landscape,dimension){}


  Persistence_landscape_on_grid_interface(const char* filename, size_t number_of_points, int dimension = -1):
  Persistence_landscape_on_grid(filename,number_of_points,dimension){}  
  
  void new_compute_average(const std::vector<Persistence_landscape_on_grid_interface*>& to_average) 
  {
	  std::vector<Persistence_landscape_on_grid*> to_average_new;
	  to_average_new.reserve( to_average.size() );
	  for ( size_t i = 0 ; i != to_average.size() ; ++i )
	  {
		  to_average_new.push_back( (Persistence_landscape_on_grid*)to_average[i] );
	  }
	  this->compute_average(to_average_new);
  }
};   

}  // namespace Persistence_representations
}  // namespace Gudhi

#endif  // PERSISTENCE_LANDSCAPE_ON_GRID_INTERFACE_H_
	
