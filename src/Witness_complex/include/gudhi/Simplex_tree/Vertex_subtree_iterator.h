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

#ifndef SIMPLEX_TREE_VERTEX_SUBTREE_ITERATOR_H_
#define SIMPLEX_TREE_VERTEX_SUBTREE_ITERATOR_H_

#include <boost/iterator/iterator_facade.hpp>
#include <boost/version.hpp>
#if BOOST_VERSION >= 105600
# include <boost/container/static_vector.hpp>
#endif

#include <vector>

namespace Gudhi {

/* \brief Iterator over the simplices of the subtree of a given
 * vertex of the simplicial complex.
 *
 * Forward iterator, value_type is SimplexTree::Simplex_handle.*/
template<class SimplexTree>
class Simplex_tree_vertex_subtree_iterator : public boost::iterator_facade<
    Simplex_tree_vertex_subtree_iterator<SimplexTree>,
    typename SimplexTree::Simplex_handle const, boost::forward_traversal_tag> {
 public:
  typedef typename SimplexTree::Simplex_handle Simplex_handle;
  typedef typename SimplexTree::Siblings Siblings;
  typedef typename SimplexTree::Vertex_handle Vertex_handle;

// any end() iterator
  Simplex_tree_vertex_subtree_iterator()
      : sib_(nullptr),
        st_(nullptr),
        dim_skel_(0),
        curr_dim_(0) {
  }

  Simplex_tree_vertex_subtree_iterator(SimplexTree * st, Vertex_handle v, int dim_skel)
      : sib_(nullptr),
        st_(st),
        dim_skel_(dim_skel),
        curr_dim_(0) {
    if (st == nullptr || st->root() == nullptr || st->root()->members().empty()) {
      st_ = nullptr;
    } else {
      sh_ = st->root()->members().find(v);
      sib_ = st->root();
      if (sh_ == st->null_simplex() || sh_->second.children() == st->root())
        st_ = nullptr;
      else {
        while (st->has_children(sh_) && curr_dim_ < dim_skel_) {
          sib_ = sh_->second.children();
          sh_ = sib_->members().begin();
          ++curr_dim_;
        }
        while (curr_dim_ != dim_skel_ && st_ != nullptr)
          simple_increment();
      }
    }
  }
 private:
  friend class boost::iterator_core_access;

// valid when iterating along the SAME boundary.
  bool equal(Simplex_tree_vertex_subtree_iterator const& other) const {
    if (other.st_ == nullptr) {
      return (st_ == nullptr);
    }
    if (st_ == nullptr) {
      return false;
    }
    return (&(sh_->second) == &(other.sh_->second));
  }

  Simplex_handle const& dereference() const {
    return sh_;
  }

// Depth first traversal of the skeleton.
  void simple_increment() {
    ++sh_;
    if (sh_ == sib_->members().end()) {
      if (sib_->oncles() == nullptr || sib_->oncles() == st_->root()) { // If returned to root, then finish
        st_ = nullptr;
        return;
      }  // reach the end
      sh_ = sib_->oncles()->members().find(sib_->parent());
      sib_ = sib_->oncles();
      --curr_dim_;
      return;
    }
    while (st_->has_children(sh_) && curr_dim_ < dim_skel_) {
      sib_ = sh_->second.children();
      sh_ = sib_->members().begin();
      ++curr_dim_;
    }
  }
  
// Depth first traversal of the skeleton.
  void increment() {
    simple_increment();
    while (curr_dim_ != dim_skel_ && st_ != nullptr)
      simple_increment();
  }

  Simplex_handle sh_;
  Siblings * sib_;
  SimplexTree * st_;
  int dim_skel_;
  int curr_dim_;
};
  
/* @} */  // end addtogroup simplex_tree
}  // namespace Gudhi

#endif  // SIMPLEX_TREE_VERTEX_SUBTREE_ITERATOR_H_
