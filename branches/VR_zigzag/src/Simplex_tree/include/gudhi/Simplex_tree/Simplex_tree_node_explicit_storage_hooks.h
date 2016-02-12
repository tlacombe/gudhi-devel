/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Clément Maria
 *
 *    Copyright (C) 2015  Clément Maria
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

#ifndef SRC_SIMPLEX_TREE_INCLUDE_GUDHI_SIMPLEX_TREE_SIMPLEX_TREE_NODE_EXPLICIT_STORAGE_HOOKS_H_
#define SRC_SIMPLEX_TREE_INCLUDE_GUDHI_SIMPLEX_TREE_SIMPLEX_TREE_NODE_EXPLICIT_STORAGE_HOOKS_H_

#include <vector>
#include <boost/intrusive/list.hpp>

namespace Gudhi {

/* \addtogroup simplex_tree
 * Represents a node of a Simplex_tree.
 * @{
 */

/* Hook for an intrusive list linking all node with label v in a subtree rooted 
 * at a vertex node u together (cofaces of {u,v}). */
// struct st_node_cof_hook_tag;
// typedef boost::intrusive::list_base_hook 
//             < boost::intrusive::tag < st_node_cof_hook_tag >  
//             , boost::intrusive::link_mode < boost::intrusive::auto_unlink > 
//             >      base_hook_st_node_cof;
/* Hook for all node with label u (cofaces of vertex u). */
struct st_node_lab_hook_tag;
typedef boost::intrusive::list_base_hook 
            < boost::intrusive::tag < st_node_lab_hook_tag >  
            , boost::intrusive::link_mode < boost::intrusive::auto_unlink > 
            >      base_hook_st_node_lab;

/*
 * \brief Node of a simplex tree with filtration value
 * and simplex key.
 *
 * It stores explicitely its own filtration value and its own Simplex_key.
 */
template<class SimplexTree >
class Simplex_tree_node_explicit_storage_hooks
: public base_hook_st_node_lab {
 public:
  typedef typename SimplexTree::Siblings         Siblings;
  typedef typename SimplexTree::Filtration_value Filtration_value;
  typedef typename SimplexTree::Simplex_key      Simplex_key;

  // Default constructor.
  // Simplex_tree_node_explicit_storage_hooks()
  //     : children_(NULL)
  //       // simplex_key_(-1),
  //     , filtration_(0) 
  //     , annotation_(NULL) {
  // }

  Simplex_tree_node_explicit_storage_hooks( Siblings * sib,
                                            Filtration_value filtration)
      : children_(sib)
        // simplex_key_(-1),
      , filtration_(filtration) 
      , annotation_(NULL) {}

  void assign_key(Simplex_key key) {
    simplex_key_ = key;
  }

  /*
   * Assign a children to the node
   */
  void assign_children(Siblings * children) {
    children_ = children;
  }
  /*
   *
   */
  void assign_filtration(double filtration_value) {
    filtration_ = filtration_value;
  }

  Filtration_value filtration() {
    return filtration_;
  }

  /* Careful -> children_ can be NULL*/
  Siblings * children() {
    return children_;
  }

  Simplex_key key() {
    return simplex_key_;
  }

  void * annotation() { return annotation_; }
  
  /** \brief Checks if two Nodes are equal. */
  bool operator==(const Simplex_tree_node_explicit_storage_hooks& other_node) {
    if ((children_ != other_node.children_) ||  // For children, just check the pointed value is the same or not
        (filtration_ != other_node.filtration_) ||
        (simplex_key_ != other_node.simplex_key_))
      return false;
    return true;
  }

  /** \brief Checks if two simplex trees are different. */
  bool operator!=(const Simplex_tree_node_explicit_storage_hooks& other_node) {
    return (!(*this == other_node));
  }

public: // private:
  Siblings * children_;

  // Data attached to simplex, explicit storage
  Simplex_key simplex_key_;
  Filtration_value filtration_;   // value in the filtration
  void * annotation_;
};

/* @} */  // end addtogroup simplex_tree
}  // namespace Gudhi

#endif  // SRC_SIMPLEX_TREE_INCLUDE_GUDHI_SIMPLEX_TREE_SIMPLEX_TREE_NODE_EXPLICIT_STORAGE_HOOKS_H_
