/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Clément Maria
 *
 *    Copyright (C) 2014  INRIA Sophia Antipolis-Méditerranée (France)
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

#ifndef SRC_SIMPLEX_TREE_INCLUDE_GUDHI_SIMPLEX_TREE_SIMPLEX_TREE_NODE_EXPLICIT_STORAGE_H_
#define SRC_SIMPLEX_TREE_INCLUDE_GUDHI_SIMPLEX_TREE_SIMPLEX_TREE_NODE_EXPLICIT_STORAGE_H_

#include <vector>

namespace Gudhi {

/* \addtogroup simplex_tree
 * Represents a node of a Simplex_tree.
 * @{
 */

/*
 * \brief Node of a simplex tree with filtration value
 * and simplex key.
 *
 * It stores explicitely its own filtration value and its own Simplex_key.
 */
template<class SimplexTree>
struct Simplex_tree_node_explicit_storage : SimplexTree::Simplex_data {
  typedef typename SimplexTree::Siblings Siblings;
  typedef typename SimplexTree::Simplex_data Simplex_data;

  // Default constructor.
  Simplex_tree_node_explicit_storage()
      : children_(NULL) {
	  this->set_key(-1);
	  this->set_filtration(0);
  }

  Simplex_tree_node_explicit_storage(Siblings * sib)
      : children_(sib) {
	  this->set_key(-1);
  }

  Simplex_tree_node_explicit_storage(Siblings * sib,
                                     Filtration_value filtration)
      : children_(sib) {
	  this->set_key(-1);
	  this->set_filtration(filtration);
  }

  /*
   * Assign a children to the node
   */
  void assign_children(Siblings * children) {
    children_ = children;
  }

  /* Careful -> children_ can be NULL*/
  Siblings * children() {
    return children_;
  }

 private:
  Siblings * children_;
};

/* @} */  // end addtogroup simplex_tree
}  // namespace Gudhi

#endif  // SRC_SIMPLEX_TREE_INCLUDE_GUDHI_SIMPLEX_TREE_SIMPLEX_TREE_NODE_EXPLICIT_STORAGE_H_
