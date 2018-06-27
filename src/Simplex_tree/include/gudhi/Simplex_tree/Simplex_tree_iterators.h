/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Clément Maria
 *
 *    Copyright (C) 2014 Inria
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

#ifndef SIMPLEX_TREE_SIMPLEX_TREE_ITERATORS_H_
#define SIMPLEX_TREE_SIMPLEX_TREE_ITERATORS_H_

#include <gudhi/Debug_utils.h>

#include <boost/iterator/iterator_facade.hpp>
#include <boost/version.hpp>
#if BOOST_VERSION >= 105600
#include <boost/container/static_vector.hpp>
#endif

#include <vector>
#include <boost/iterator/filter_iterator.hpp>

namespace Gudhi {

/* \addtogroup simplex_tree
 * Iterators and range types for the Simplex_tree.
 * @{
 */

/* \brief Iterator over the vertices of a simplex
 * in a SimplexTree.
 *
 * Forward iterator, 'value_type' is SimplexTree::Vertex_handle.*/
template <class SimplexTree>
class Simplex_tree_simplex_vertex_iterator
    : public boost::iterator_facade<Simplex_tree_simplex_vertex_iterator<SimplexTree>,
                                    typename SimplexTree::Vertex_handle const, boost::forward_traversal_tag,
                                    typename SimplexTree::Vertex_handle const> {
 public:
  typedef typename SimplexTree::Simplex_handle Simplex_handle;
  typedef typename SimplexTree::Siblings Siblings;
  typedef typename SimplexTree::Vertex_handle Vertex_handle;

  explicit Simplex_tree_simplex_vertex_iterator(SimplexTree* st)
      :  // any end() iterator
        sib_(nullptr),
        v_(st->null_vertex()) {}

  Simplex_tree_simplex_vertex_iterator(SimplexTree* st, Simplex_handle sh)
      : sib_(st->self_siblings(sh)), v_(sh->first) {}

 private:
  friend class boost::iterator_core_access;

  bool equal(Simplex_tree_simplex_vertex_iterator const& other) const { return sib_ == other.sib_ && v_ == other.v_; }

  Vertex_handle const& dereference() const { return v_; }

  void increment() {
    v_ = sib_->parent();
    sib_ = sib_->oncles();
  }

  Siblings* sib_;
  Vertex_handle v_;
};

/*---------------------------------------------------------------------------*/
/* \brief Iterator over the simplices of the boundary of a
 *  simplex.
 *
 * Forward iterator, value_type is SimplexTree::Simplex_handle.*/
template <class SimplexTree>
class Simplex_tree_boundary_simplex_iterator
    : public boost::iterator_facade<Simplex_tree_boundary_simplex_iterator<SimplexTree>,
                                    typename SimplexTree::Simplex_handle const, boost::forward_traversal_tag> {
 public:
  typedef typename SimplexTree::Simplex_handle Simplex_handle;
  typedef typename SimplexTree::Vertex_handle Vertex_handle;
  typedef typename SimplexTree::Siblings Siblings;

  // any end() iterator
  explicit Simplex_tree_boundary_simplex_iterator(SimplexTree* st) : sib_(nullptr), sh_(st->null_simplex()), st_(st) {}

  template <class SimplexHandle>
  Simplex_tree_boundary_simplex_iterator(SimplexTree* st, SimplexHandle sh) : last_(sh->first), sib_(nullptr), st_(st) {
    // Only check once at the beginning instead of for every increment, as this is expensive.
    if (SimplexTree::Options::contiguous_vertices)
      GUDHI_CHECK(st_->contiguous_vertices(), "The set of vertices is not { 0, ..., n } without holes");
    Siblings* sib = st->self_siblings(sh);
    next_ = sib->parent();
    sib_ = sib->oncles();
    if (sib_ != nullptr) {
      if (SimplexTree::Options::contiguous_vertices && sib_->oncles() == nullptr){
        // Only relevant for edges
          sh_ = sib_->members_.begin();
          for (int i = 0; i < next_; i++){
              sh_++;
          }
      }
      else
        sh_ = sib_->find(next_);
    } else {
      sh_ = st->null_simplex();
    }  // vertex: == end()
  }

 private:
  friend class boost::iterator_core_access;
  // valid when iterating along the SAME boundary.
  bool equal(Simplex_tree_boundary_simplex_iterator const& other) const { return sh_ == other.sh_; }

  Simplex_handle const& dereference() const {
    assert(sh_ != st_->null_simplex());
    return sh_;
  }

  void increment() {
    if (sib_ == nullptr) {
      sh_ = st_->null_simplex();
      return;
    }

    Siblings* for_sib = sib_;
    Siblings* new_sib = sib_->oncles();
    auto rit = suffix_.rbegin();
    if (SimplexTree::Options::contiguous_vertices && new_sib == nullptr) {
      // We reached the root, use a short-cut to find a vertex.
      if (rit == suffix_.rend()) {
        // Segment, this vertex is the last boundary simplex
          sh_ = for_sib->members_.begin();
          for (int i = 0; i < last_; i++){
              sh_++;
          }
        //sh_ = for_sib->members_.begin() + last_;
        sib_ = nullptr;
        return;
      } else {
        // Dim >= 2, initial step of the descent
          sh_ = for_sib->members_.begin();
          for (int i = 0; i < *rit; i++){
              sh_++;
          }
        //sh_ = for_sib->members_.begin() + *rit;
        for_sib = sh_->second.children();
        ++rit;
      }
    }
    for (; rit != suffix_.rend(); ++rit) {
      sh_ = for_sib->find(*rit);
      for_sib = sh_->second.children();
    }
    sh_ = for_sib->find(last_);  // sh_ points to the right simplex now
    suffix_.push_back(next_);
    next_ = sib_->parent();
    sib_ = new_sib;
  }

  // Most of the storage should be moved to the range, iterators should be light.
  Vertex_handle last_;  // last vertex of the simplex
  Vertex_handle next_;  // next vertex to push in suffix_
#if BOOST_VERSION >= 105600
  // 40 seems a conservative bound on the dimension of a Simplex_tree for now,
  // as it would not fit on the biggest hard-drive.
  boost::container::static_vector<Vertex_handle, 40> suffix_;
  // static_vector still has some overhead compared to a trivial hand-made
  // version using std::aligned_storage, or compared to making suffix_ static.
#else
  std::vector<Vertex_handle> suffix_;
#endif
  Siblings* sib_;      // where the next search will start from
  Simplex_handle sh_;  // current Simplex_handle in the boundary
  SimplexTree* st_;    // simplex containing the simplicial complex
};
/*---------------------------------------------------------------------------*/
/* \brief Iterator over the simplices of a simplicial complex.
 *
 * Forward iterator, value_type is SimplexTree::Simplex_handle.*/
template <class SimplexTree>
class Simplex_tree_complex_simplex_iterator
    : public boost::iterator_facade<Simplex_tree_complex_simplex_iterator<SimplexTree>,
                                    typename SimplexTree::Simplex_handle const, boost::forward_traversal_tag> {
 public:
  typedef typename SimplexTree::Simplex_handle Simplex_handle;
  typedef typename SimplexTree::Siblings Siblings;
  typedef typename SimplexTree::Vertex_handle Vertex_handle;

  // any end() iterator
  Simplex_tree_complex_simplex_iterator() : sib_(nullptr), st_(nullptr) {}

  explicit Simplex_tree_complex_simplex_iterator(SimplexTree* st) : sib_(nullptr), st_(st) {
    if (st == nullptr || st->root() == nullptr || st->root()->members().empty()) {
      st_ = nullptr;
    } else {
      sh_ = st->root()->members().begin();
      sib_ = st->root();
      while (st->has_children(sh_)) {
        sib_ = sh_->second.children();
        sh_ = sib_->members().begin();
      }
    }
  }

 private:
  friend class boost::iterator_core_access;

  // valid when iterating along the SAME boundary.
  bool equal(Simplex_tree_complex_simplex_iterator const& other) const {
    if (other.st_ == nullptr) {
      return (st_ == nullptr);
    }
    if (st_ == nullptr) {
      return false;
    }
    return (&(sh_->second) == &(other.sh_->second));
  }

  Simplex_handle const& dereference() const { return sh_; }

  // Depth first traversal.
  void increment() {
    ++sh_;
    if (sh_ == sib_->members().end()) {
      if (sib_->oncles() == nullptr) {
        st_ = nullptr;
        return;
      }  // reach the end
      sh_ = sib_->oncles()->members().find(sib_->parent());
      sib_ = sib_->oncles();
      return;
    }
    while (st_->has_children(sh_)) {
      sib_ = sh_->second.children();
      sh_ = sib_->members().begin();
    }
  }

  Simplex_handle sh_;
  Siblings* sib_;
  SimplexTree* st_;
};

/* \brief Iterator over the simplices of the skeleton of a given
 * dimension of the simplicial complex.
 *
 * Forward iterator, value_type is SimplexTree::Simplex_handle.*/
template <class SimplexTree>
class Simplex_tree_skeleton_simplex_iterator
    : public boost::iterator_facade<Simplex_tree_skeleton_simplex_iterator<SimplexTree>,
                                    typename SimplexTree::Simplex_handle const, boost::forward_traversal_tag> {
 public:
  typedef typename SimplexTree::Simplex_handle Simplex_handle;
  typedef typename SimplexTree::Siblings Siblings;
  typedef typename SimplexTree::Vertex_handle Vertex_handle;

  // any end() iterator
  Simplex_tree_skeleton_simplex_iterator() : sib_(nullptr), st_(nullptr), dim_skel_(0), curr_dim_(0) {}

  Simplex_tree_skeleton_simplex_iterator(SimplexTree* st, int dim_skel)
      : sib_(nullptr), st_(st), dim_skel_(dim_skel), curr_dim_(0) {
    if (st == nullptr || st->root() == nullptr || st->root()->members().empty()) {
      st_ = nullptr;
    } else {
      sh_ = st->root()->members().begin();
      sib_ = st->root();
      while (st->has_children(sh_) && curr_dim_ < dim_skel_) {
        sib_ = sh_->second.children();
        sh_ = sib_->members().begin();
        ++curr_dim_;
      }
    }
  }

 private:
  friend class boost::iterator_core_access;

  // valid when iterating along the SAME boundary.
  bool equal(Simplex_tree_skeleton_simplex_iterator const& other) const {
    if (other.st_ == nullptr) {
      return (st_ == nullptr);
    }
    if (st_ == nullptr) {
      return false;
    }
    return (&(sh_->second) == &(other.sh_->second));
  }

  Simplex_handle const& dereference() const { return sh_; }

  // Depth first traversal of the skeleton.
  void increment() {
    ++sh_;
    if (sh_ == sib_->members().end()) {
      if (sib_->oncles() == nullptr) {
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

  Simplex_handle sh_;
  Siblings* sib_;
  SimplexTree* st_;
  int dim_skel_;
  int curr_dim_;
};

/* \brief Iterator over all the roots of subtrees containing cofaces
 * of a given simplex.
 *
 * Forward iterator, value_type is SimplexTree::Simplex_handle.*/
template <class SimplexTree>
class Simplex_tree_opt_cofaces_rooted_subtrees_simplex_iterator
    : public boost::iterator_facade<Simplex_tree_opt_cofaces_rooted_subtrees_simplex_iterator<SimplexTree>,
                                    typename SimplexTree::Simplex_handle const, boost::forward_traversal_tag> {
 public:
  typedef typename SimplexTree::Simplex_handle Simplex_handle;
  typedef typename SimplexTree::Siblings Siblings;
  typedef typename SimplexTree::Vertex_handle Vertex_handle;
  typedef typename SimplexTree::Node Node;

  //----------------------------------------------
  /* Predicate to check whether an input SimplexTree::Node represents a
   * coface of codimension codim_ of a simplex simp_,
   * stored as a std::vector of SimplexTree::Vertex_handle.
   * If the codimension codim is 0, returns true for all cofaces.
   *
   * Given a SimplexHandle in a simplex tree cpx_, traverse the tree upwards
   * to find the sequence of Vertex_handle of simp_. Does not test sh itself.
   * Used for filter_iterator in the optimized algorithm for
   * cofaces_simplex_range.
   */
  class is_coface {
   public:
    is_coface() : cpx_(NULL) {}
    is_coface(SimplexTree* cpx, std::vector<Vertex_handle> simp) : cpx_(cpx), simp_(simp) {}

    // Returns true iff traversing the Node upwards to the root reads a
    // coface of simp_ of codimension codim_
    bool operator()(typename SimplexTree::Hooks_simplex_base& curr_hooks) {
      int dim = 0;
      Node& curr_node = static_cast<Node&>(curr_hooks);
      auto vertex_it = simp_.begin();
      // first Node must always have label simp_.begin()
      auto curr_sib = cpx_->self_siblings(curr_node, *vertex_it);
      if (++vertex_it == simp_.end()) {
        return true;
      }
      while (curr_sib->oncles() != NULL) {  // todo is NULL valid?
        if (curr_sib->parent() == *vertex_it) {
          if (++vertex_it == simp_.end()) {  // we have a coface root
            return true;
          }
        }
        curr_sib = curr_sib->oncles();
        ++dim;
      }
      return false;
    }

   private:
    SimplexTree* cpx_;
    std::vector<Vertex_handle> simp_;  // vertices of simplex, reverse ordered
  };

  typedef boost::filter_iterator<is_coface, typename SimplexTree::List_max_vertex::iterator>
      Filtered_cofaces_simplex_iterator;
  // any end() iterator
  Simplex_tree_opt_cofaces_rooted_subtrees_simplex_iterator() : predicate_(), st_(nullptr) {}

  Simplex_tree_opt_cofaces_rooted_subtrees_simplex_iterator(SimplexTree* cpx, std::vector<Vertex_handle> simp)
      : predicate_(cpx, simp), st_(cpx) {
    max_v_ = *(simp.begin());
    auto list_ptr = st_->cofaces_data_structure_.access(max_v_);
    // assert(list_it != st_->cofaces_data_structure_.end());
    it_ = boost::make_filter_iterator(predicate_, list_ptr->begin(), list_ptr->end());
    end_ = boost::make_filter_iterator(predicate_, list_ptr->end(), list_ptr->end());
    Node& curr_node = static_cast<Node&>(*it_);
    auto curr_sib = st_->self_siblings(curr_node, max_v_);
    sh_ = curr_sib->find(max_v_);
  }

 private:
  friend class boost::iterator_core_access;

  // valid when iterating along the SAME list of max vertex.
  bool equal(Simplex_tree_opt_cofaces_rooted_subtrees_simplex_iterator const& other) const {
    if (other.st_ == nullptr) {
      return (st_ == nullptr);
    }
    if (st_ == nullptr) {
      return false;
    }
    return (&(*it_) == &(*(other.it_)));
  }

  Simplex_handle const& dereference() const { return sh_; }

  void increment() {
    if (++it_ == end_) {
      st_ = nullptr;
    }       //== end
    else {  // update sh_
      Node& curr_node = static_cast<Node&>(*it_);
      auto curr_sib = st_->self_siblings(curr_node, max_v_);
      sh_ = curr_sib->find(max_v_);
    }
  }

  is_coface predicate_;
  SimplexTree* st_;
  Filtered_cofaces_simplex_iterator it_;
  Filtered_cofaces_simplex_iterator end_;
  Vertex_handle max_v_;
  Simplex_handle sh_;  // Simplex_handle corresponding to Node pointed at by it_
};

/* \brief Iterator over all cofaces of dimension exactly max_dim_cofaces_ for
 * a simplex.
 *
 * uses hooks stored in the Node of the SimplexTree.*/
template <class SimplexTree>
class Simplex_tree_opt_cofaces_simplex_iterator
    : public boost::iterator_facade<Simplex_tree_opt_cofaces_simplex_iterator<SimplexTree>,
                                    typename SimplexTree::Simplex_handle const, boost::forward_traversal_tag> {
 public:
  typedef typename SimplexTree::Simplex_handle Simplex_handle;
  typedef typename SimplexTree::Siblings Siblings;
  typedef typename SimplexTree::Vertex_handle Vertex_handle;
  typedef typename SimplexTree::Node Node;

  // any end() iterator
  Simplex_tree_opt_cofaces_simplex_iterator() : st_(nullptr) {}

  Simplex_tree_opt_cofaces_simplex_iterator(SimplexTree* cpx, std::vector<Vertex_handle> simp, int codim)
      : st_(cpx),
        it_(cpx, simp),
        end_(),
        sh_(*it_),
        dim_root_(st_->dimension(sh_)),
        dim_sh_(dim_root_),
        sib_(st_->self_siblings(sh_)) {
    if (codim == 0) {
      max_dim_cofaces_ = cpx->dimension();
      all_cofaces = true;
    }       // all cofaces
    else {  // exact dimension
      max_dim_cofaces_ = codim + simp.size() - 1;
      all_cofaces = false;
      while (st_ != nullptr && dim_sh_ != max_dim_cofaces_) {
        increment_bounded_dimension(); //look for a coface of right codimension
      } 
    }
  }

  int dim_sh() { return dim_sh_; }

 private:
  friend class boost::iterator_core_access;

  // valid when iterating along the SAME list of max vertex.
  bool equal(Simplex_tree_opt_cofaces_simplex_iterator const& other) const {
    if (other.st_ == nullptr) {
      return (st_ == nullptr);
    }
    if (st_ == nullptr) {
      return false;
    }
    return (&(sh_->second) == &(other.sh_->second));
  }

  Simplex_handle const& dereference() const { return sh_; }

  void increment_bounded_dimension() {
    if (dim_root_ == dim_sh_) {  // sh_ points at a subtree root,
      // which is a valid coface of the right dimension.
      if (!st_->has_children(sh_) || dim_root_ == max_dim_cofaces_) {  // go to next subtree
        do {
          if (++it_ == end_) {
            st_ = nullptr;
            return;
          }
          sh_ = *it_;  // initialize fields for new subtree
          dim_root_ = st_->dimension(sh_);
          dim_sh_ = dim_root_;
          sib_ = st_->self_siblings(sh_);
        } while (dim_root_ > max_dim_cofaces_);
        return;
      } else {  // sh_ has children, of dimension <= dim_cofaces
        sib_ = sh_->second.children();
        sh_ = sib_->members().begin();
        ++dim_sh_;
        return;
      }
    } else { // we are inside a subtree
      if (!st_->has_children(sh_) || dim_sh_ == max_dim_cofaces_) { //++sh_
        ++sh_;
        while (sh_ == sib_->members().end()) {
          if (dim_sh_ == dim_root_ + 1) {  // go to next subtree with dim_root_ <= dim_cofaces
            do {
              if (++it_ == end_) {
                st_ = nullptr;
                return;
              }
              sh_ = *it_;  // initialize fields for new subtree
              dim_root_ = st_->dimension(sh_);
              dim_sh_ = dim_root_;
              sib_ = st_->self_siblings(sh_);
            } while (dim_root_ > max_dim_cofaces_);
            return;
          } else {  // continue going up
            Vertex_handle parent = sib_->parent();
            sib_ = sib_->oncles();
            sh_ = sib_->members().find(parent);
            ++sh_;
            --dim_sh_;
          }
        }
        return;
      } else {  // sh_ has children, of dimension <= dim_cofaces: go down DFS
        sib_ = sh_->second.children();
        sh_ = sib_->members().begin();
        ++dim_sh_;
        return;
      }
    }
  }

  void increment() {
    if (!all_cofaces) {
      do {
        increment_bounded_dimension();
      } while (st_ != nullptr && dim_sh_ != max_dim_cofaces_);
    } else {
      increment_bounded_dimension();
    }
  }

  SimplexTree* st_;
  Simplex_tree_opt_cofaces_rooted_subtrees_simplex_iterator<SimplexTree> it_;
  Simplex_tree_opt_cofaces_rooted_subtrees_simplex_iterator<SimplexTree> end_;
  Simplex_handle sh_;    // curr
  int dim_root_;         // dimension of simplex *it_
  int dim_sh_;           // dimension of simplex sh_
  int max_dim_cofaces_;  // bound on the dimension of the cofaces
  Siblings* sib_;
  bool all_cofaces;
};

/* @} */  // end addtogroup simplex_tree
}  // namespace Gudhi

#endif  // SIMPLEX_TREE_SIMPLEX_TREE_ITERATORS_H_
