/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Siargey Kachanovich
 *
 *    Copyright (C) 2017  INRIA (France)
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

#ifndef SIB_VERTEX_PAIR_H_
#define SIB_VERTEX_PAIR_H_

#include <utility>

namespace Gudhi {

namespace witness_complex {

template < class Simplex_tree_,
           class Vertex_ >
class Sib_vertex_pair : public std::pair<typename Simplex_tree_::Siblings*, Vertex_> {
public:

  Sib_vertex_pair()
  {}
  
  Sib_vertex_pair(typename Simplex_tree_::Siblings* sib, Vertex_ v)
    : std::pair<typename Simplex_tree_::Siblings*, Vertex_>(sib,v)
  {
  }
  
  typename Simplex_tree_::Simplex_handle simplex_handle() const
  {
    return this->first->members().find(this->second);
  }
};

}
}
  
#endif
