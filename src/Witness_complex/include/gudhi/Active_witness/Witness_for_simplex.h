/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Siargey Kachanovich
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

#ifndef WITNESS_FOR_SIMPLEX_H_
#define WITNESS_FOR_SIMPLEX_H_

#include <vector>
#include <utility>

namespace Gudhi {

namespace witness_complex {

  /* \class Witness_for_simplex
   *  \brief Class that contains all the necessary information about witnessing a simplex by a given witness
   *  \details Contains the farthest simplex' landmark in the nearest landmark list, the pointer to the witness
   *           and the limit distance defined by the closest non-simplex landmark to the witness.
   */
template< typename Active_witness_iterator,
          typename Active_witness_list >
class Witness_for_simplex {
public:  
  typedef typename Active_witness_list::iterator AWL_iterator;

  Active_witness_iterator last_it_;
  AWL_iterator witness_;
  double limit_distance_;

  Witness_for_simplex()
  {
  }
  
  Witness_for_simplex(Active_witness_iterator last_it, AWL_iterator witness, double limit_distance)
    : last_it_(last_it), witness_(witness), limit_distance_(limit_distance)
  {
  }
 
};

}
}
  
#endif
