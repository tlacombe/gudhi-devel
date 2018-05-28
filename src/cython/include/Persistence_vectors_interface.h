/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Pawel Dlotko
 *
 *    Copyright (C) 2016  INRIA (France)
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

#ifndef PERSISTENCE_VECTORS_INTERFACE_H_
#define PERSISTENCE_VECTORS_INTERFACE_H_

// gudhi include
#include <gudhi/Persistence_vectors.h>

namespace Gudhi {
namespace Persistence_representations {


class Vector_distances_in_diagram_interface : public Vector_distances_in_diagram<Euclidean_distance> {
 public:
  Vector_distances_in_diagram_interface():Vector_distances_in_diagram(){}

  Vector_distances_in_diagram_interface(const std::vector<std::pair<double, double> >& intervals, size_t where_to_cut):
  Vector_distances_in_diagram(intervals,where_to_cut){}

  Vector_distances_in_diagram_interface(const char* filename, size_t where_to_cut,
                              unsigned dimension = std::numeric_limits<unsigned>::max()):
  Vector_distances_in_diagram_interface(filename,where_to_cut,dimension){}                             
  
  void compute_average_interface(const std::vector<Vector_distances_in_diagram_interface*>& to_average)
  {
	  std::vector<Vector_distances_in_diagram*> to_average_new;
	  to_average_new.reserve( to_average.size() );
	  for ( size_t i = 0 ; i != to_average.size() ; ++i )
	  {
		  to_average_new.push_back( (Vector_distances_in_diagram_interface*)to_average[i] );
	  }
	  this->compute_average(to_average_new);
  }  
};

}  // namespace Persistence_representations
}  // namespace Gudhi

#endif  // PERSISTENCE_VECTORS_INTERFACE_H_
