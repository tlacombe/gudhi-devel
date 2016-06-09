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

#ifndef SIMPLEX_TREE_H_
#define SIMPLEX_TREE_H_

#include <gudhi/Simplex_tree/Simplex_tree_node_explicit_storage.h>
#include <gudhi/Simplex_tree/Simplex_tree_siblings.h>
#include <gudhi/Simplex_tree/Simplex_tree_iterators.h>
#include <gudhi/Simplex_tree/indexing_tag.h>

#include <boost/intrusive/list.hpp>
#include <boost/iterator/filter_iterator.hpp>

// #include <gudhi/Simplex_tree_zz/Simplex_tree_node_explicit_storage_hooks.h>
// #include <gudhi/Simplex_tree_zz/Simplex_tree_siblings_zigzag.h>
#include <gudhi/Simplex_tree/Simplex_tree_zigzag_classes.h>

#include <gudhi/reader_utils.h>
#include <gudhi/graph_simplicial_complex.h>
#include <gudhi/Debug_utils.h>

#include <boost/fusion/container/map.hpp>
#include <boost/fusion/include/map.hpp>
#include <boost/fusion/container/map/map_fwd.hpp>
#include <boost/fusion/include/map_fwd.hpp>

#include <boost/container/flat_map.hpp>
#include <boost/iterator/transform_iterator.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/range/adaptor/reversed.hpp>

#ifdef GUDHI_USE_TBB
#include <tbb/parallel_sort.h>
#endif

#include <utility>
#include <vector>
#include <functional>  // for greater<>
#include <stdexcept>
#include <limits>  // Inf
#include <initializer_list>
#include <algorithm>  // for std::max

namespace Gudhi {
/** \defgroup simplex_tree Filtered Complexes
 * \author    Cl&eacute;ment Maria
 *
 * A simplicial complex \f$\mathbf{K}\f$
 * on a set of vertices \f$V = \{1, \cdots ,|V|\}\f$ is a collection of simplices
 * \f$\{\sigma\}\f$,
 * \f$\sigma \subseteq V\f$ such that \f$\tau \subseteq \sigma \in \mathbf{K} \rightarrow \tau \in
 * \mathbf{K}\f$. The
 * dimension \f$n=|\sigma|-1\f$ of \f$\sigma\f$ is its number of elements minus \f$1\f$.
 *
 * A filtration of a simplicial complex is
 * a function \f$f:\mathbf{K} \rightarrow \mathbb{R}\f$ satisfying \f$f(\tau)\leq f(\sigma)\f$ whenever
 * \f$\tau \subseteq \sigma\f$. Ordering the simplices by increasing filtration values
 * (breaking ties so as a simplex appears after its subsimplices of same filtration value)
 * provides an indexing scheme.
 *

 <DT>Implementations:</DT>
 There are two implementation of complexes. The first on is the Simplex_tree data structure.
 The simplex tree is an efficient and flexible
 data structure for representing general (filtered) simplicial complexes. The data structure
 is described in \cite boissonnatmariasimplextreealgorithmica
 \image html "Simplex_tree_representation.png" "Simplex tree representation"

 The second one is the Hasse_complex. The Hasse complex is a data structure representing
 explicitly all co-dimension 1 incidence relations in a complex. It is consequently faster
 when accessing the boundary of a simplex, but is less compact and harder to construct from
 scratch.

 * \copyright GNU General Public License v3.
 * @{
 */

struct Simplex_tree_options_full_featured;

/**
 * \class Simplex_tree Simplex_tree.h gudhi/Simplex_tree.h
 * \brief Simplex Tree data structure for representing simplicial complexes.
 *
 * \details Every simplex \f$[v_0, \cdots ,v_d]\f$ admits a canonical orientation
 * induced by the order relation on vertices \f$ v_0 < \cdots < v_d \f$.
 *
 * Details may be found in \cite boissonnatmariasimplextreealgorithmica.
 *
 * \implements FilteredComplex
 *
 */

template<typename SimplexTreeOptions = Simplex_tree_options_full_featured>
class Simplex_tree {
 public:
  typedef SimplexTreeOptions Options;
  typedef typename Options::Indexing_tag Indexing_tag;
  /** \brief Type for the value of the filtration function.
   *
   * Must be comparable with <. */
  typedef typename Options::Filtration_value Filtration_value;
  /** \brief Key associated to each simplex.
   *
   * Must be a signed integer type. */
  typedef typename Options::Simplex_key Simplex_key;
  /** \brief Type for the vertex handle.
   *
   * Must be a signed integer type. It admits a total order <. */
  typedef typename Options::Vertex_handle Vertex_handle;

  /* Type of node in the simplex tree. */
  typedef Simplex_tree_node_explicit_storage<Simplex_tree> Node;

  /* Type of dictionary Vertex_handle -> Node for traversing the simplex tree. */
  // Note: this wastes space when Vertex_handle is 32 bits and Node is aligned on 64 bits. It would be better to use a
  // flat_set (with our own comparator) where we can control the layout of the struct (put Vertex_handle and
  // Simplex_key next to each other).
  typedef typename boost::container::flat_map<Vertex_handle, Node> Dictionary;
  /* Need a fast deletion operation for zigzag persistence.*/
  // typedef boost::containe::map<Vertex_handle, Node> Dictionary;

  /* \brief Set of nodes sharing a same parent in the simplex tree. */
  typedef Simplex_tree_siblings<Simplex_tree, Dictionary> Siblings;
  // typedef Simplex_tree_siblings_zigzag<Simplex_tree, Dictionary> Siblings; //to do is this necessary?



  struct Key_simplex_base_real {
    Key_simplex_base_real() : key_(-1) {}
    void assign_key(Simplex_key k) { key_ = k; }
    Simplex_key key() const { return key_; }
   private:
    Simplex_key key_;
  };
  struct Key_simplex_base_dummy {
    Key_simplex_base_dummy() {}
    void assign_key(Simplex_key) { }
    Simplex_key key() const { assert(false); return -1; }
  };
  typedef typename std::conditional<Options::store_key, Key_simplex_base_real, Key_simplex_base_dummy>::type
      Key_simplex_base;

  struct Filtration_simplex_base_real {
    Filtration_simplex_base_real() : filt_(0) {}
    void assign_filtration(Filtration_value f) { filt_ = f; }
    Filtration_value filtration() const { return filt_; }
   private:
    Filtration_value filt_;
  };
  struct Filtration_simplex_base_dummy {
    Filtration_simplex_base_dummy() {}
    void assign_filtration(Filtration_value f) { assert(f == 0); }
    Filtration_value filtration() const { return 0; }
  };
  typedef typename std::conditional<Options::store_filtration, Filtration_simplex_base_real,
    Filtration_simplex_base_dummy>::type Filtration_simplex_base;

/* 
 * Stores members hooks in Node class in order to maintain them
 * in a map for fast coface location.
 */

  typedef boost::intrusive::list_member_hook< //allows .unlink()
              boost::intrusive::link_mode< boost::intrusive::auto_unlink > 
                                            >                    Member_hook_t;
  struct Hooks_simplex_base_dummy {
  };
  struct Hooks_simplex_base_cofaces {
    Hooks_simplex_base_cofaces() {}
    //the copy constructor, inherited by the Node class, exchanges hooks, 
    //and make the ones of this invalid.
    //this is used when stored in a map like DS, using copies when 
    //performing insertions and rebalancing of the rbtree
    Hooks_simplex_base_cofaces(const Hooks_simplex_base_cofaces & other)
    { list_max_vertex_hook_.swap_nodes(other.list_max_vertex_hook_); }

    void unlink_hooks() { list_max_vertex_hook_.unlink(); }

  public:
    mutable Member_hook_t          list_max_vertex_hook_;
  };

  typedef boost::intrusive::member_hook< 
            Hooks_simplex_base_cofaces
          , Member_hook_t
          , &Hooks_simplex_base_cofaces::list_max_vertex_hook_ > 
                                                            List_member_hook_t;

  typedef boost::intrusive::list< Hooks_simplex_base_cofaces
                                , List_member_hook_t 
                                , boost::intrusive::constant_time_size<false> 
                                >                              List_max_vertex;

  typedef typename std::conditional< Options::link_simplices_through_max_vertex
                                   , Hooks_simplex_base_cofaces
                                   , Hooks_simplex_base_dummy >::type 
                                                            Hooks_simplex_base;

  // struct Annotation_simplex_base_real {}                                                            

 public:
  /** \brief Handle type to a simplex contained in the simplicial complex represented
   * by the simplex tree. */
  typedef typename Dictionary::iterator Simplex_handle;

 private:
  typedef typename Dictionary::iterator Dictionary_it;
  typedef typename Dictionary_it::value_type Dit_value_t;

  struct return_first {
    Vertex_handle operator()(const Dit_value_t& p_sh) const {
      return p_sh.first;
    }
  };

 public:
  /** \name Range and iterator types
   *
   * The naming convention is Container_content_(iterator/range). A Container_content_range is
   * essentially an object on which the methods begin() and end() can be called. They both return
   * an object of type Container_content_iterator, and allow the traversal of the range
   * [ begin();end() ).
   * @{ */

  /** \brief Iterator over the vertices of the simplicial complex.
   *
   * 'value_type' is Vertex_handle. */
  typedef boost::transform_iterator<return_first, Dictionary_it> Complex_vertex_iterator;
  /** \brief Range over the vertices of the simplicial complex. */
  typedef boost::iterator_range<Complex_vertex_iterator> Complex_vertex_range;
  /** \brief Iterator over the vertices of a simplex.
   *
   * 'value_type' is Vertex_handle. */
  typedef Simplex_tree_simplex_vertex_iterator<Simplex_tree> Simplex_vertex_iterator;
  /** \brief Range over the vertices of a simplex. */
  typedef boost::iterator_range<Simplex_vertex_iterator> Simplex_vertex_range;
  /** \brief Range over the cofaces of a simplex. */
  typedef std::vector<Simplex_handle> Cofaces_simplex_range;
  /** \brief Iterator over the simplices of the boundary of a simplex.
   *
   * 'value_type' is Simplex_handle. */
  typedef Simplex_tree_boundary_simplex_iterator<Simplex_tree> Boundary_simplex_iterator;
  /** \brief Range over the simplices of the boundary of a simplex. */
  typedef boost::iterator_range<Boundary_simplex_iterator> Boundary_simplex_range;
  /** \brief Iterator over the simplices of the simplicial complex.
   *
   * 'value_type' is Simplex_handle. */
  typedef Simplex_tree_complex_simplex_iterator<Simplex_tree> Complex_simplex_iterator;
  /** \brief Range over the simplices of the simplicial complex. */
  typedef boost::iterator_range<Complex_simplex_iterator> Complex_simplex_range;
  /** \brief Iterator over the simplices of the skeleton of the simplicial complex, for a given 
   * dimension.
   *
   * 'value_type' is Simplex_handle. */
  typedef Simplex_tree_skeleton_simplex_iterator<Simplex_tree> Skeleton_simplex_iterator;
  /** \brief Range over the simplices of the skeleton of the simplicial complex, for a given 
   * dimension. */
  typedef boost::iterator_range<Skeleton_simplex_iterator> Skeleton_simplex_range;
  /** \brief Range over the simplices of the simplicial complex, ordered by the filtration. */
  typedef std::vector<Simplex_handle> Filtration_simplex_range;
  /** \brief Iterator over the simplices of the simplicial complex, ordered by the filtration.
   *
   * 'value_type' is Simplex_handle. */
  typedef typename Filtration_simplex_range::const_iterator Filtration_simplex_iterator;

  /* @} */  // end name range and iterator types
  /** \name Range and iterator methods
   * @{ */

  /** \brief Returns a range over the vertices of the simplicial complex. 
   * The order is increasing according to < on Vertex_handles.*/
  Complex_vertex_range complex_vertex_range() {
    return Complex_vertex_range(
                                boost::make_transform_iterator(root_.members_.begin(), return_first()),
                                boost::make_transform_iterator(root_.members_.end(), return_first()));
  }

  /** \brief Returns a range over the simplices of the simplicial complex.
   *
   * In the Simplex_tree, the tree is traverse in a depth-first fashion.
   * Consequently, simplices are ordered according to lexicographic order on the list of
   * Vertex_handles of a simplex, read in increasing < order for Vertex_handles. */
  Complex_simplex_range complex_simplex_range() {
    return Complex_simplex_range(Complex_simplex_iterator(this),
                                 Complex_simplex_iterator());
  }

  /** \brief Returns a range over the simplices of the dim-skeleton of the simplicial complex.
   *
   * The \f$d\f$-skeleton of a simplicial complex \f$\mathbf{K}\f$ is the simplicial complex containing the
   * simplices of \f$\mathbf{K}\f$ of dimension at most \f$d\f$.
   *
   * @param[in] dim The maximal dimension of the simplices in the skeleton.
   *
   * The simplices are ordered according to lexicographic order on the list of
   * Vertex_handles of a simplex, read in increasing < order for Vertex_handles. */
  Skeleton_simplex_range skeleton_simplex_range(int dim) {
    return Skeleton_simplex_range(Skeleton_simplex_iterator(this, dim),
                                  Skeleton_simplex_iterator());
  }

  /** \brief Returns a range over the simplices of the simplicial complex,
   * in the order of the filtration.
   *
   * The filtration is a monotonic function \f$ f: \mathbf{K} \rightarrow \mathbb{R} \f$, i.e. if two simplices
   * \f$\tau\f$ and \f$\sigma\f$ satisfy \f$\tau \subseteq \sigma\f$ then
   * \f$f(\tau) \leq f(\sigma)\f$.
   *
   * The method returns simplices ordered according to increasing filtration values. Ties are
   * resolved by considering inclusion relation (subsimplices appear before their cofaces). If two
   * simplices have same filtration value but are not comparable w.r.t. inclusion, lexicographic
   * order is used.
   *
   * The filtration must be valid. If the filtration has not been initialized yet, the
   * method initializes it (i.e. order the simplices). If the complex has changed since the last time the filtration
   * was initialized, please call `initialize_filtration()` to recompute it. */
  Filtration_simplex_range const& filtration_simplex_range(Indexing_tag = Indexing_tag()) {
    if (filtration_vect_.empty()) {
      initialize_filtration();
    }
    return filtration_vect_;
  }

  /** \brief Returns a range over the vertices of a simplex.
   *
   * The order in which the vertices are visited is the decreasing order for < on Vertex_handles,
   * which is consequenlty
   * equal to \f$(-1)^{\text{dim} \sigma}\f$ the canonical orientation on the simplex.
   */
  Simplex_vertex_range simplex_vertex_range(Simplex_handle sh) {
    assert(sh != null_simplex());  // Empty simplex
    return Simplex_vertex_range(Simplex_vertex_iterator(this, sh),
                                Simplex_vertex_iterator(this));
  }

  /** \brief Returns a range over the simplices of the boundary of a simplex.
   *
   * The boundary of a simplex is the set of codimension \f$1\f$ subsimplices of the simplex.
   * If the simplex is \f$[v_0, \cdots ,v_d]\f$, with canonical orientation
   * induced by \f$ v_0 < \cdots < v_d \f$, the iterator enumerates the
   * simplices of the boundary in the order:
   * \f$[v_0,\cdots,\widehat{v_i},\cdots,v_d]\f$ for \f$i\f$ from \f$0\f$ to \f$d\f$,
   * where \f$\widehat{v_i}\f$ means that the vertex \f$v_i\f$ is omitted.
   *
   * We note that the alternate sum of the simplices given by the iterator
   * gives \f$(-1)^{\text{dim} \sigma}\f$ the chains corresponding to the boundary
   * of the simplex.
   *
   * @param[in] sh Simplex for which the boundary is computed. */
  template<class SimplexHandle>
  Boundary_simplex_range boundary_simplex_range(SimplexHandle sh) {
    return Boundary_simplex_range(Boundary_simplex_iterator(this, sh),
                                  Boundary_simplex_iterator(this));
  }

  /** @} */  // end range and iterator methods
  /** \name Constructor/Destructor
   * @{ */

  /** \brief Constructs an empty simplex tree. */
  Simplex_tree()
      : null_vertex_(-1),
      threshold_(0),
      root_(nullptr, null_vertex_),
      filtration_vect_(),
      dimension_(-1) { }

  /** \brief User-defined copy constructor reproduces the whole tree structure. */
  Simplex_tree(const Simplex_tree& simplex_source)
      : null_vertex_(simplex_source.null_vertex_),
      threshold_(simplex_source.threshold_),
      root_(nullptr, null_vertex_ , simplex_source.root_.members_),
      filtration_vect_(),
      dimension_(simplex_source.dimension_) {
    auto root_source = simplex_source.root_;
    rec_copy(&root_, &root_source);
  }

  /** \brief depth first search, inserts simplices when reaching a leaf. */
  void rec_copy(Siblings *sib, Siblings *sib_source) {
    for (auto sh = sib->members().begin(), sh_source = sib_source->members().begin();
         sh != sib->members().end(); ++sh, ++sh_source) {
      if (has_children(sh_source)) {
        Siblings * newsib = new Siblings(sib, sh_source->first);
        newsib->members_.reserve(sh_source->second.children()->members().size());
        for (auto & child : sh_source->second.children()->members())
          newsib->members_.emplace_hint(newsib->members_.end(), child.first, Node(newsib, child.second.filtration()));
        rec_copy(newsib, sh_source->second.children());
        sh->second.assign_children(newsib);
      }
    }
  }

  /** \brief User-defined move constructor moves the whole tree structure. */
  Simplex_tree(Simplex_tree && old)
      : null_vertex_(std::move(old.null_vertex_)),
      threshold_(std::move(old.threshold_)),
      root_(std::move(old.root_)),
      filtration_vect_(std::move(old.filtration_vect_)),
      dimension_(std::move(old.dimension_)) {
    old.dimension_ = -1;
    old.threshold_ = 0;
    old.root_ = Siblings(nullptr, null_vertex_);
  }

  /** \brief Destructor; deallocates the whole tree structure. */
  ~Simplex_tree() {
    for (auto sh = root_.members().begin(); sh != root_.members().end(); ++sh) {
      if (has_children(sh)) {
        rec_delete(sh->second.children());
      }
    }
  }
  /** @} */  // end constructor/destructor
 private:
  // Recursive deletion
  void rec_delete(Siblings * sib) {
    for (auto sh = sib->members().begin(); sh != sib->members().end(); ++sh) {
      if (has_children(sh)) {
        rec_delete(sh->second.children());
      }
    }
    delete sib;
  }

 public:
  /** \brief Checks if two simplex trees are equal. */
  bool operator==(Simplex_tree& st2) {
    if ((null_vertex_ != st2.null_vertex_) ||
        (threshold_ != st2.threshold_) ||
        (dimension_ != st2.dimension_))
      return false;
    return rec_equal(&root_, &st2.root_);
  }

  /** \brief Checks if two simplex trees are different. */
  bool operator!=(Simplex_tree& st2) {
    return (!(*this == st2));
  }

 private:
  /** rec_equal: Checks recursively whether or not two simplex trees are equal, using depth first search. */
  bool rec_equal(Siblings* s1, Siblings* s2) {
    if (s1->members().size() != s2->members().size())
      return false;
    for (auto sh1 = s1->members().begin(), sh2 = s2->members().begin();
         (sh1 != s1->members().end() && sh2 != s2->members().end()); ++sh1, ++sh2) {
      if (sh1->first != sh2->first || sh1->second.filtration() != sh2->second.filtration())
        return false;
      if (has_children(sh1) != has_children(sh2))
        return false;
      // Recursivity on children only if both have children
      else if (has_children(sh1))
        if (!rec_equal(sh1->second.children(), sh2->second.children()))
          return false;
    }
    return true;
  }

 public:
  /** \brief Returns the key associated to a simplex.
   *
   * The filtration must be initialized.
   * \pre SimplexTreeOptions::store_key
   */
  static Simplex_key key(Simplex_handle sh) {
    return sh->second.key();
  }

  /** \brief Returns the simplex associated to a key.
   *
   * The filtration must be initialized.
   * \pre SimplexTreeOptions::store_key
   */
  Simplex_handle simplex(Simplex_key key) const {
    return filtration_vect_[key];
  }

  /** \brief Returns the filtration value of a simplex.
   *
   * Called on the null_simplex, returns INFINITY.
   * If SimplexTreeOptions::store_filtration is false, returns 0.
   */
  static Filtration_value filtration(Simplex_handle sh) {
    if (sh != null_simplex()) {
      return sh->second.filtration();
    } else {
      return INFINITY;
    }
  }

  /** \brief Sets the filtration value of a simplex.
   * \exception std::invalid_argument In debug mode, if sh is a null_simplex.
   */
  void assign_filtration(Simplex_handle sh, Filtration_value fv) {
    GUDHI_CHECK(sh != null_simplex(),
                std::invalid_argument("Simplex_tree::assign_filtration - cannot assign filtration on null_simplex"));
    sh->second.assign_filtration(fv);
  }

  /** \brief Returns an upper bound of the filtration values of the simplices. */
  Filtration_value filtration() const {
    return threshold_;
  }

  /** \brief Returns a Simplex_handle different from all Simplex_handles
   * associated to the simplices in the simplicial complex.
   *
   * One can call filtration(null_simplex()). */
  static Simplex_handle null_simplex() {
    return Dictionary_it(nullptr);
  }

  /** \brief Returns a key different for all keys associated to the
   * simplices of the simplicial complex. */
  static Simplex_key null_key() {
    return -1;
  }

  /** \brief Returns a Vertex_handle different from all Vertex_handles associated
   * to the vertices of the simplicial complex. */
  Vertex_handle null_vertex() const {
    return null_vertex_;
  }

  /** \brief Returns the number of vertices in the complex. */
  size_t num_vertices() const {
    return root_.members_.size();
  }

 public:
  /** \brief returns the number of simplices in the simplex_tree. */
  size_t num_simplices() {
    return num_simplices(&root_);
  }

 private:
  /** \brief returns the number of simplices in the simplex_tree. */
  size_t num_simplices(Siblings * sib) {
    auto sib_begin = sib->members().begin();
    auto sib_end = sib->members().end();
    size_t simplices_number = sib_end - sib_begin;
    for (auto sh = sib_begin; sh != sib_end; ++sh) {
      if (has_children(sh)) {
        simplices_number += num_simplices(sh->second.children());
      }
    }
    return simplices_number;
  }

 public:
  /** \brief Returns the dimension of a simplex.
   *
   * Must be different from null_simplex().*/
  int dimension(Simplex_handle sh) {
    Siblings * curr_sib = self_siblings(sh);
    int dim = 0;
    while (curr_sib != nullptr) {
      ++dim;
      curr_sib = curr_sib->oncles();
    }
    return dim - 1;
  }

  /** \brief Returns an upper bound on the dimension of the simplicial complex. */
  int dimension() const {
    return dimension_;
  }

  /** \brief Returns true if the node in the simplex tree pointed by
   * sh has children.*/
  template<class SimplexHandle>
  bool has_children(SimplexHandle sh) const {
    return (sh->second.children()->parent() == sh->first);
  }

    /** \brief Given a range of Vertex_handles, returns the Simplex_handle
   * of the simplex in the simplicial complex containing the corresponding
   * vertices. Return null_simplex() if the simplex is not in the complex.
   *
   * The type InputVertexRange must be a range of <CODE>Vertex_handle</CODE>
   * on which we can call std::begin() function
   */
  template<class InputVertexRange = std::initializer_list<Vertex_handle>>
  Simplex_handle find(const InputVertexRange & s) {
    auto first = std::begin(s);
    auto last = std::end(s);

    if (first == last)
      return null_simplex();  // ----->>

    // Copy before sorting
    std::vector<Vertex_handle> copy(first, last);
    std::sort(std::begin(copy), std::end(copy));
    return find_simplex(copy);
  }

 private:
  /** Find function, with a sorted range of vertices. */
  Simplex_handle find_simplex(const std::vector<Vertex_handle> & simplex) {
    Siblings * tmp_sib = &root_;
    Dictionary_it tmp_dit;
    Vertex_handle last = simplex.back();
    for (auto v : simplex) {
      tmp_dit = tmp_sib->members_.find(v);
      if (tmp_dit == tmp_sib->members_.end()) {
        return null_simplex();
      }
      if (!has_children(tmp_dit) && v != last) {
        return null_simplex();
      }
      tmp_sib = tmp_dit->second.children();
    }
    return tmp_dit;
  }

  /** \brief Returns the Simplex_handle corresponding to the 0-simplex
   * representing the vertex with Vertex_handle v. */
  Simplex_handle find_vertex(Vertex_handle v) {
    if (Options::contiguous_vertices) {
      assert(contiguous_vertices());
      return root_.members_.begin() + v;
    } else {
      return root_.members_.find(v);
    }
  }

 public:
  /** \private \brief Test if the vertices have contiguous numbering: 0, 1, etc.  */
  bool contiguous_vertices() const {
    if (root_.members_.empty()) return true;
    if (root_.members_.begin()->first != 0) return false;
    if (std::prev(root_.members_.end())->first != static_cast<Vertex_handle>(root_.members_.size() - 1)) return false;
    return true;
  }

 private:
  /** \brief Inserts a simplex represented by a vector of vertex.
   * @param[in]  simplex    vector of Vertex_handles, representing the vertices of the new simplex. The vector must be
   * sorted by increasing vertex handle order.
   * @param[in]  filtration the filtration value assigned to the new simplex.
   * @return If the new simplex is inserted successfully (i.e. it was not in the
   * simplicial complex yet) the bool is set to true and the Simplex_handle is the handle assigned
   * to the new simplex.
   * If the insertion fails (the simplex is already there), the bool is set to false. If the insertion
   * fails and the simplex already in the complex has a filtration value strictly bigger than 'filtration',
   * we assign this simplex with the new value 'filtration', and set the Simplex_handle field of the
   * output pair to the Simplex_handle of the simplex. Otherwise, we set the Simplex_handle part to
   * null_simplex.
   * 
  */
  std::pair<Simplex_handle, bool> insert_vertex_vector(const std::vector<Vertex_handle>& simplex,
                                                     Filtration_value filtration) {
    Siblings * curr_sib = &root_;
    std::pair<Simplex_handle, bool> res_insert;
    auto vi = simplex.begin();
    for (; vi != simplex.end() - 1; ++vi) {
      res_insert = curr_sib->members_.emplace(*vi, Node(curr_sib, filtration));
      if (!(has_children(res_insert.first))) {
        res_insert.first->second.assign_children(new Siblings(curr_sib, *vi));
      }
      curr_sib = res_insert.first->second.children();
    }
    res_insert = curr_sib->members_.emplace(*vi, Node(curr_sib, filtration));
    if (!res_insert.second) {
      // if already in the complex
      if (res_insert.first->second.filtration() > filtration) {
        // if filtration value modified
        res_insert.first->second.assign_filtration(filtration);
        return res_insert;
      }
      // if filtration value unchanged
      return std::pair<Simplex_handle, bool>(null_simplex(), false);
    }
    // otherwise the insertion has succeeded
    return res_insert;
  }

 public:
  /** \brief Insert a simplex, represented by a range of Vertex_handles, in the simplicial complex.
   *
   * @param[in]  simplex    range of Vertex_handles, representing the vertices of the new simplex
   * @param[in]  filtration the filtration value assigned to the new simplex.
   * @return If the new simplex is inserted successfully (i.e. it was not in the
   * simplicial complex yet) the bool is set to true and the Simplex_handle is the handle assigned
   * to the new simplex.
   * If the insertion fails (the simplex is already there), the bool is set to false. If the insertion
   * fails and the simplex already in the complex has a filtration value strictly bigger than 'filtration',
   * we assign this simplex with the new value 'filtration', and set the Simplex_handle field of the
   * output pair to the Simplex_handle of the simplex. Otherwise, we set the Simplex_handle part to
   * null_simplex.
   *
   * All subsimplices do not necessary need to be already in the simplex tree to proceed to an
   * insertion. However, the property of being a simplicial complex will be violated. This allows
   * us to insert a stream of simplices contained in a simplicial complex without considering any
   * order on them.
   *
   * The filtration value
   * assigned to the new simplex must preserve the monotonicity of the filtration.
   *
   * The type InputVertexRange must be a range for which .begin() and
   * .end() return input iterators, with 'value_type' Vertex_handle. */
  template<class InputVertexRange = std::initializer_list<Vertex_handle>>
  std::pair<Simplex_handle, bool> insert_simplex(const InputVertexRange & simplex,
                                                 Filtration_value filtration = 0) {
    auto first = std::begin(simplex);
    auto last = std::end(simplex);

    if (first == last)
      return std::pair<Simplex_handle, bool>(null_simplex(), true);  // ----->>

    // Copy before sorting
    std::vector<Vertex_handle> copy(first, last);
    std::sort(std::begin(copy), std::end(copy));
    return insert_vertex_vector(copy, filtration);
  }

    /** \brief Insert a N-simplex and all his subfaces, from a N-simplex represented by a range of
   * Vertex_handles, in the simplicial complex.
   *
   * @param[in]  Nsimplex   range of Vertex_handles, representing the vertices of the new N-simplex
   * @param[in]  filtration the filtration value assigned to the new N-simplex.
   * @return If the new simplex is inserted successfully (i.e. it was not in the
   * simplicial complex yet) the bool is set to true and the Simplex_handle is the handle assigned
   * to the new simplex.
   * If the insertion fails (the simplex is already there), the bool is set to false. If the insertion
   * fails and the simplex already in the complex has a filtration value strictly bigger than 'filtration',
   * we assign this simplex with the new value 'filtration', and set the Simplex_handle field of the
   * output pair to the Simplex_handle of the simplex. Otherwise, we set the Simplex_handle part to
   * null_simplex.
   */
  template<class InputVertexRange = std::initializer_list<Vertex_handle>>
  std::pair<Simplex_handle, bool> insert_simplex_and_subfaces(const InputVertexRange& Nsimplex,
                                   Filtration_value filtration = 0) {
    auto first = std::begin(Nsimplex);
    auto last = std::end(Nsimplex);

    if (first == last)
      return std::pair<Simplex_handle, bool>(null_simplex(), true);  // ----->>

    // Copy before sorting
    std::vector<Vertex_handle> copy(first, last);
    std::sort(std::begin(copy), std::end(copy));

    std::vector<std::vector<Vertex_handle>> to_be_inserted;
    std::vector<std::vector<Vertex_handle>> to_be_propagated;
    return rec_insert_simplex_and_subfaces(copy, to_be_inserted, to_be_propagated, filtration);
  }

 private:
  std::pair<Simplex_handle, bool> rec_insert_simplex_and_subfaces(std::vector<Vertex_handle>& the_simplex,
                                                             std::vector<std::vector<Vertex_handle>>& to_be_inserted,
                                                             std::vector<std::vector<Vertex_handle>>& to_be_propagated,
                                                             Filtration_value filtration = 0.0) {
    std::pair<Simplex_handle, bool> insert_result;
    if (the_simplex.size() > 1) {
      // Get and remove last vertex
      Vertex_handle last_vertex = the_simplex.back();
      the_simplex.pop_back();
      // Recursive call after last vertex removal
      insert_result = rec_insert_simplex_and_subfaces(the_simplex, to_be_inserted, to_be_propagated, filtration);

      // Concatenation of to_be_inserted and to_be_propagated
      to_be_inserted.insert(to_be_inserted.begin(), to_be_propagated.begin(), to_be_propagated.end());
      to_be_propagated = to_be_inserted;

      // to_be_inserted treatment
      for (auto& simplex_tbi : to_be_inserted) {
        simplex_tbi.push_back(last_vertex);
      }
      std::vector<Vertex_handle> last_simplex(1, last_vertex);
      to_be_inserted.insert(to_be_inserted.begin(), last_simplex);
      // i.e. (0,1,2) =>
      // [to_be_inserted | to_be_propagated] = [(1) (0,1) | (0)]
      // [to_be_inserted | to_be_propagated] = [(2) (0,2) (1,2) (0,1,2) | (0) (1) (0,1)]
      // N.B. : it is important the last inserted to be the highest in dimension
      // in order to return the "last" insert_simplex result

      // insert all to_be_inserted
      for (auto& simplex_tbi : to_be_inserted) {
        insert_result = insert_vertex_vector(simplex_tbi, filtration);
      }
    } else if (the_simplex.size() == 1) {
      // When reaching the end of recursivity, vector of simplices shall be empty and filled on back recursive
      if ((to_be_inserted.size() != 0) || (to_be_propagated.size() != 0)) {
        std::cerr << "Simplex_tree::rec_insert_simplex_and_subfaces - Error vector not empty\n";
        exit(-1);
      }
      std::vector<Vertex_handle> first_simplex(1, the_simplex.back());
      // i.e. (0,1,2) => [to_be_inserted | to_be_propagated] = [(0) | ]
      to_be_inserted.push_back(first_simplex);

      insert_result = insert_vertex_vector(first_simplex, filtration);
    } else {
        std::cerr << "Simplex_tree::rec_insert_simplex_and_subfaces - Recursivity error\n";
        exit(-1);
    }
    return insert_result;
  }

 public:
  /** \brief Assign a value 'key' to the key of the simplex
   * represented by the Simplex_handle 'sh'. */
  void assign_key(Simplex_handle sh, Simplex_key key) {
    sh->second.assign_key(key);
  }

  /** Returns the two Simplex_handle corresponding to the endpoints of
   * and edge. sh must point to a 1-dimensional simplex. This is an
   * optimized version of the boundary computation. */
  std::pair<Simplex_handle, Simplex_handle> endpoints(Simplex_handle sh) {
    assert(dimension(sh) == 1);
    return { find_vertex(sh->first), find_vertex(self_siblings(sh)->parent()) };
  }

  /** Returns the Siblings containing a simplex.*/
  template<class SimplexHandle>
  Siblings* self_siblings(SimplexHandle sh) {
    if (sh->second.children()->parent() == sh->first)
      return sh->second.children()->oncles();
    else
      return sh->second.children();
  }

 public:
  /** Returns a pointer to the root nodes of the simplex tree. */
  Siblings * root() {
    return &root_;
  }

  /** Set an upper bound for the filtration values. */
  void set_filtration(Filtration_value fil) {
    threshold_ = fil;
  }

  /** Set a dimension for the simplicial complex. */
  void set_dimension(int dimension) {
    dimension_ = dimension;
  }

 public:
  /** \brief Initializes the filtrations, i.e. sort the
   * simplices according to their order in the filtration and initializes all Simplex_keys.
   *
   * After calling this method, filtration_simplex_range() becomes valid, and each simplex is
   * assigned a Simplex_key corresponding to its order in the filtration (from 0 to m-1 for a
   * simplicial complex with m simplices).
   *
   * Will be automatically called when calling filtration_simplex_range()
   * if the filtration has never been initialized yet. */
  void initialize_filtration() {
    filtration_vect_.clear();
    filtration_vect_.reserve(num_simplices());
    for (Simplex_handle sh : complex_simplex_range())
      filtration_vect_.push_back(sh);

    /* We use stable_sort here because with libstdc++ it is faster than sort.
     * is_before_in_filtration is now a total order, but we used to call
     * stable_sort for the following heuristic:
     * The use of a depth-first traversal of the simplex tree, provided by
     * complex_simplex_range(), combined with a stable sort is meant to
     * optimize the order of simplices with same filtration value. The
     * heuristic consists in inserting the cofaces of a simplex as soon as
     * possible.
     */
#ifdef GUDHI_USE_TBB
    tbb::parallel_sort(filtration_vect_.begin(), filtration_vect_.end(), is_before_in_filtration(this));
#else
    std::stable_sort(filtration_vect_.begin(), filtration_vect_.end(), is_before_in_filtration(this));
#endif
  }

 private:
  /** Recursive search of cofaces
   * This function uses DFS
   *\param vertices contains a list of vertices, which represent the vertices of the simplex not found yet.
   *\param curr_nbVertices represents the number of vertices of the simplex we reached by going through the tree.
   *\param cofaces contains a list of Simplex_handle, representing all the cofaces asked.
   *\param star true if we need the star of the simplex
   *\param nbVertices number of vertices of the cofaces we search
   * Prefix actions : When the bottom vertex matches with the current vertex in the tree, we remove the bottom vertex from vertices.
   * Infix actions : Then we call or not the recursion.
   * Postfix actions : Finally, we add back the removed vertex into vertices, and remove this vertex from curr_nbVertices so that we didn't change the parameters.
   * If the vertices list is empty, we need to check if curr_nbVertices matches with the dimension of the cofaces asked.
   */
  void rec_coface(std::vector<Vertex_handle> &vertices, Siblings *curr_sib, int curr_nbVertices,
                  std::vector<Simplex_handle>& cofaces, bool star, int nbVertices) {
    if (!(star || curr_nbVertices <= nbVertices))  // dimension of actual simplex <= nbVertices
      return;
    for (Simplex_handle simplex = curr_sib->members().begin(); simplex != curr_sib->members().end(); ++simplex) {
      if (vertices.empty()) {
        // If we reached the end of the vertices, and the simplex has more vertices than the given simplex
        // => we found a coface

        // Add a coface if we wan't the star or if the number of vertices of the current simplex matches with nbVertices
        bool addCoface = (star || curr_nbVertices == nbVertices);
        if (addCoface)
          cofaces.push_back(simplex);
        if ((!addCoface || star) && has_children(simplex))  // Rec call
          rec_coface(vertices, simplex->second.children(), curr_nbVertices + 1, cofaces, star, nbVertices);
      } else {
        if (simplex->first == vertices.back()) {
          // If curr_sib matches with the top vertex
          bool equalDim = (star || curr_nbVertices == nbVertices);  // dimension of actual simplex == nbVertices
          bool addCoface = vertices.size() == 1 && equalDim;
          if (addCoface)
            cofaces.push_back(simplex);
          if ((!addCoface || star) && has_children(simplex)) {
            // Rec call
            Vertex_handle tmp = vertices.back();
            vertices.pop_back();
            rec_coface(vertices, simplex->second.children(), curr_nbVertices + 1, cofaces, star, nbVertices);
            vertices.push_back(tmp);
          }
        } else if (simplex->first > vertices.back()) {
          return;
        } else {
          // (simplex->first < vertices.back()
          if (has_children(simplex))
            rec_coface(vertices, simplex->second.children(), curr_nbVertices + 1, cofaces, star, nbVertices);
        }
      }
    }
  }

 public:
  /** \brief Compute the star of a n simplex
   * \param simplex represent the simplex of which we search the star
   * \return Vector of Simplex_handle, empty vector if no cofaces found.
   */

  Cofaces_simplex_range star_simplex_range(const Simplex_handle simplex) {
    return cofaces_simplex_range(simplex, 0);
  }

  /** \brief Compute the cofaces of a n simplex
   * \param simplex represent the n-simplex of which we search the n+codimension cofaces
   * \param codimension The function returns the n+codimension-cofaces of the n-simplex. If codimension = 0, 
   * return all cofaces (equivalent of star function)
   * \return Vector of Simplex_handle, empty vector if no cofaces found.
   */

  Cofaces_simplex_range cofaces_simplex_range(const Simplex_handle simplex, int codimension) {
    Cofaces_simplex_range cofaces;
    // codimension must be positive or null integer
    assert(codimension >= 0);
    Simplex_vertex_range rg = simplex_vertex_range(simplex);
    std::vector<Vertex_handle> copy(rg.begin(), rg.end());
    if (codimension + static_cast<int>(copy.size()) > dimension_ + 1 ||
        (codimension == 0 && static_cast<int>(copy.size()) > dimension_))  // n+codimension greater than dimension_
      return cofaces;
    // must be sorted in decreasing order
    assert(std::is_sorted(copy.begin(), copy.end(), std::greater<Vertex_handle>()));
    bool star = codimension == 0;
    rec_coface(copy, &root_, 1, cofaces, star, codimension + static_cast<int>(copy.size()));
    return cofaces;
  }

 private:
  /** \brief Returns true iff the list of vertices of sh1
   * is smaller than the list of vertices of sh2 w.r.t.
   * lexicographic order on the lists read in reverse.
   *
   * It defines a StrictWeakOrdering on simplices. The Simplex_vertex_iterators
   * must traverse the Vertex_handle in decreasing order. Reverse lexicographic order satisfy
   * the property that a subsimplex of a simplex is always strictly smaller with this order. */
  bool reverse_lexicographic_order(Simplex_handle sh1, Simplex_handle sh2) {
    Simplex_vertex_range rg1 = simplex_vertex_range(sh1);
    Simplex_vertex_range rg2 = simplex_vertex_range(sh2);
    Simplex_vertex_iterator it1 = rg1.begin();
    Simplex_vertex_iterator it2 = rg2.begin();
    while (it1 != rg1.end() && it2 != rg2.end()) {
      if (*it1 == *it2) {
        ++it1;
        ++it2;
      } else {
        return *it1 < *it2;
      }
    }
    return ((it1 == rg1.end()) && (it2 != rg2.end()));
  }

  /** \brief StrictWeakOrdering, for the simplices, defined by the filtration.
   *
   * It corresponds to the partial order
   * induced by the filtration values, with ties resolved using reverse lexicographic order.
   * Reverse lexicographic order has the property to always consider the subsimplex of a simplex
   * to be smaller. The filtration function must be monotonic. */
  struct is_before_in_filtration {
    explicit is_before_in_filtration(Simplex_tree * st)
        : st_(st) { }

    bool operator()(const Simplex_handle sh1, const Simplex_handle sh2) const {
      // Not using st_->filtration(sh1) because it uselessly tests for null_simplex.
      if (sh1->second.filtration() != sh2->second.filtration()) {
        return sh1->second.filtration() < sh2->second.filtration();
      }
      // is sh1 a proper subface of sh2
      return st_->reverse_lexicographic_order(sh1, sh2);
    }

    Simplex_tree * st_;
  };

 public:
  /** \brief Inserts a 1-skeleton in an empty Simplex_tree.
   *
   * The Simplex_tree must contain no simplex when the method is
   * called.
   *
   * Inserts all vertices and edges given by a OneSkeletonGraph.
   * OneSkeletonGraph must be a model of boost::AdjacencyGraph,
   * boost::EdgeListGraph and boost::PropertyGraph.
   *
   * The vertex filtration value is accessible through the property tag
   * vertex_filtration_t.
   * The edge filtration value is accessible through the property tag
   * edge_filtration_t.
   *
   * boost::graph_traits<OneSkeletonGraph>::vertex_descriptor
   *                                    must be Vertex_handle.
   * boost::graph_traits<OneSkeletonGraph>::directed_category
   *                                    must be undirected_tag. */
  template<class OneSkeletonGraph>
  void insert_graph(const OneSkeletonGraph& skel_graph) {
    // the simplex tree must be empty
    assert(num_simplices() == 0);

    if (boost::num_vertices(skel_graph) == 0) {
      return;
    }
    if (num_edges(skel_graph) == 0) {
      dimension_ = 0;
    } else {
      dimension_ = 1;
    }

    root_.members_.reserve(boost::num_vertices(skel_graph));

    typename boost::graph_traits<OneSkeletonGraph>::vertex_iterator v_it,
        v_it_end;
    for (std::tie(v_it, v_it_end) = boost::vertices(skel_graph); v_it != v_it_end;
         ++v_it) {
      root_.members_.emplace_hint(
                                  root_.members_.end(), *v_it,
                                  Node(&root_, boost::get(vertex_filtration_t(), skel_graph, *v_it)));
    }
    typename boost::graph_traits<OneSkeletonGraph>::edge_iterator e_it,
        e_it_end;
    for (std::tie(e_it, e_it_end) = boost::edges(skel_graph); e_it != e_it_end;
         ++e_it) {
      auto u = source(*e_it, skel_graph);
      auto v = target(*e_it, skel_graph);
      if (u < v) {
        // count edges only once { std::swap(u,v); } // u < v
        auto sh = find_vertex(u);
        if (!has_children(sh)) {
          sh->second.assign_children(new Siblings(&root_, sh->first));
        }

        sh->second.children()->members().emplace(
                                                 v,
                                                 Node(sh->second.children(),
                                                      boost::get(edge_filtration_t(), skel_graph, *e_it)));
      }
    }
  }

  /** \brief Expands the Simplex_tree containing only its one skeleton
   * until dimension max_dim.
   *
   * The expanded simplicial complex until dimension \f$d\f$
   * attached to a graph \f$G\f$ is the maximal simplicial complex of
   * dimension at most \f$d\f$ admitting the graph \f$G\f$ as \f$1\f$-skeleton.
   * The filtration value assigned to a simplex is the maximal filtration
   * value of one of its edges.
   *
   * The Simplex_tree must contain no simplex of dimension bigger than
   * 1 when calling the method. */
  void expansion(int max_dim) {
    dimension_ = max_dim;
    for (Dictionary_it root_it = root_.members_.begin();
         root_it != root_.members_.end(); ++root_it) {
      if (has_children(root_it)) {
        siblings_expansion(root_it->second.children(), max_dim - 1);
      }
    }
    dimension_ = max_dim - dimension_;
  }

 private:
  /** \brief Recursive expansion of the simplex tree.*/
  void siblings_expansion(Siblings * siblings,  // must contain elements
                          int k) {
    if (dimension_ > k) {
      dimension_ = k;
    }
    if (k == 0)
      return;
    Dictionary_it next = siblings->members().begin();
    ++next;

    static std::vector<std::pair<Vertex_handle, Node> > inter;  // static, not thread-safe.
    for (Dictionary_it s_h = siblings->members().begin();
         s_h != siblings->members().end(); ++s_h, ++next) {
      Simplex_handle root_sh = find_vertex(s_h->first);
      if (has_children(root_sh)) {
        intersection(
                     inter,  // output intersection
                     next,  // begin
                     siblings->members().end(),  // end
                     root_sh->second.children()->members().begin(),
                     root_sh->second.children()->members().end(),
                     s_h->second.filtration());
        if (inter.size() != 0) {
          Siblings * new_sib = new Siblings(siblings,  // oncles
                                            s_h->first,  // parent
                                            inter);  // boost::container::ordered_unique_range_t
          inter.clear();
          s_h->second.assign_children(new_sib);
          siblings_expansion(new_sib, k - 1);
        } else {
          // ensure the children property
          s_h->second.assign_children(siblings);
          inter.clear();
        }
      }
    }
  }

  /** \brief Intersects Dictionary 1 [begin1;end1) with Dictionary 2 [begin2,end2)
   * and assigns the maximal possible Filtration_value to the Nodes. */
  static void intersection(std::vector<std::pair<Vertex_handle, Node> >& intersection,
                           Dictionary_it begin1, Dictionary_it end1,
                           Dictionary_it begin2, Dictionary_it end2,
                           Filtration_value filtration_) {
    if (begin1 == end1 || begin2 == end2)
      return;  // ----->>
    while (true) {
      if (begin1->first == begin2->first) {
        Filtration_value filt = (std::max)({begin1->second.filtration(), begin2->second.filtration(), filtration_});
        intersection.emplace_back(begin1->first, Node(nullptr, filt));
        if (++begin1 == end1 || ++begin2 == end2)
          return;  // ----->>
      } else if (begin1->first < begin2->first) {
        if (++begin1 == end1)
          return;
      } else /* begin1->first > begin2->first */ {
        if (++begin2 == end2)
          return;  // ----->>
      }
    }
  }

 public:
  /** \brief Write the hasse diagram of the simplicial complex in os.
   *
   * Each row in the file correspond to a simplex. A line is written:
   * dim idx_1 ... idx_k fil   where dim is the dimension of the simplex,
   * idx_1 ... idx_k are the row index (starting from 0) of the simplices of the boundary
   * of the simplex, and fil is its filtration value. */
  void print_hasse(std::ostream& os) {
    os << num_simplices() << " " << std::endl;
    for (auto sh : filtration_simplex_range()) {
      os << dimension(sh) << " ";
      for (auto b_sh : boundary_simplex_range(sh)) {
        os << key(b_sh) << " ";
      }
      os << filtration(sh) << " \n";
    }
  }

 public:
  /** \brief Browse the simplex tree to ensure the filtration is not decreasing.
   * The simplex tree is browsed starting from the root until the leaf, and the filtration values are set with their
   * parent value (increased), in case the values are decreasing.
   * @return The filtration modification information.
   * \post Some simplex tree functions require the filtration to be valid. `make_filtration_non_decreasing()`
   * function is not launching `initialize_filtration()` but returns the filtration modification information. If the
   * complex has changed , please call `initialize_filtration()` to recompute it.
   */
  bool make_filtration_non_decreasing() {
    bool modified = false;
    // Loop must be from the end to the beginning, as higher dimension simplex are always on the left part of the tree
    for (auto& simplex : boost::adaptors::reverse(root_.members())) {
      if (has_children(&simplex)) {
        modified |= rec_make_filtration_non_decreasing(simplex.second.children());
      }
    }
    return modified;
  }

 private:
  /** \brief Recursively Browse the simplex tree to ensure the filtration is not decreasing.
   * @param[in] sib Siblings to be parsed.
   * @return The filtration modification information in order to trigger initialize_filtration.
   */
  bool rec_make_filtration_non_decreasing(Siblings * sib) {
    bool modified = false;

    // Loop must be from the end to the beginning, as higher dimension simplex are always on the left part of the tree
    for (auto& simplex : boost::adaptors::reverse(sib->members())) {
      // Find the maximum filtration value in the border
      Boundary_simplex_range boundary = boundary_simplex_range(&simplex);
      Boundary_simplex_iterator max_border = std::max_element(std::begin(boundary), std::end(boundary),
                                                              [](Simplex_handle sh1, Simplex_handle sh2) {
                                                                return filtration(sh1) < filtration(sh2);
                                                              });

      Filtration_value max_filt_border_value = filtration(*max_border);
      if (simplex.second.filtration() < max_filt_border_value) {
        // Store the filtration modification information
        modified = true;
        simplex.second.assign_filtration(max_filt_border_value);
      }
      if (has_children(&simplex)) {
        modified |= rec_make_filtration_non_decreasing(simplex.second.children());
      }
    }
    // Make the modified information to be traced by upper call
    return modified;
  }

 public:
  /** \brief Prune above filtration value given as parameter.
   * @param[in] filtration Maximum threshold value.
   * @return The filtration modification information.
   * \post Some simplex tree functions require the filtration to be valid. `prune_above_filtration()`
   * function is not launching `initialize_filtration()` but returns the filtration modification information. If the
   * complex has changed , please call `initialize_filtration()` to recompute it.
   */
  bool prune_above_filtration(Filtration_value filtration) {
    return rec_prune_above_filtration(root(), filtration);
  }

 private:
  bool rec_prune_above_filtration(Siblings* sib, Filtration_value filt) {
    auto&& list = sib->members();
    auto last = std::remove_if(list.begin(), list.end(), [=](Dit_value_t& simplex) {
        if (simplex.second.filtration() <= filt) return false;
        if (has_children(&simplex)) rec_delete(simplex.second.children());
        return true;
      });

    bool modified = (last != list.end());
    if (last == list.begin() && sib != root()) {
      // Removing the whole siblings, parent becomes a leaf.
      sib->oncles()->members()[sib->parent()].assign_children(sib->oncles());
      delete sib;
      return true;
    } else {
      // Keeping some elements of siblings. Remove the others, and recurse in the remaining ones.
      list.erase(last, list.end());
      for (auto&& simplex : list)
        if (has_children(&simplex))
          modified |= rec_prune_above_filtration(simplex.second.children(), filt);
    }
    return modified;
  }

 public:
  /** \brief Remove a maximal simplex.
   * @param[in] sh Simplex handle on the maximal simplex to remove.
   * \pre Please check the simplex has no coface before removing it.
   * \exception std::invalid_argument In debug mode, if sh has children.
   * \post Be aware that removing is shifting data in a flat_map (initialize_filtration to be done).
   */
  void remove_maximal_simplex(Simplex_handle sh) {
    // Guarantee the simplex has no children
    GUDHI_CHECK(!has_children(sh),
                std::invalid_argument("Simplex_tree::remove_maximal_simplex - argument has children"));

    // Simplex is a leaf, it means the child is the Siblings owning the leaf
    Siblings* child = sh->second.children();

    if ((child->size() > 1) || (child == root())) {
      // Not alone, just remove it from members
      // Special case when child is the root of the simplex tree, just remove it from members
      child->erase(sh);
    } else {
      // Sibling is emptied : must be deleted, and its parent must point on his own Sibling
      child->oncles()->members().at(child->parent()).assign_children(child->oncles());
      delete child;
    }
  }


/***********************************************/
/********* Zigzag flag complex methods *********/
/***********************************************/
/* 
 * These methods are dedicated to a stream-like construction of zigzag filtrations 
 * based on flag complexes. The complexes represented are always flag complexes. 
 * Consequently, the full zigzag filtration is encoded in the sequence of insertions 
 * and deletions of edges. 
 *
 * The major difference with the standard simplex_insertion(...) method of the 
 * Simplex_tree is that we maintain in parallel a data structure for locating 
 * cofaces and removing easily an edge (and its cofaces) from the complex.
 *
 * Vertex_handle must be contigue (or at least not widely spread...)
 * Must know an upper bound on the Vertex_handle to assign the coface table 
 */
public:
/** Type of edges for representing implicetely the flag zigzag filtration.*/
  typedef Zigzag_edge< Simplex_tree >                                Edge_type;
/** Forward iterator on the simplices (insertion and deletion) of the zigzag 
  * filtration.
  *
  * 'value_type' is Simplex_handle.
  */ 
  typedef Flagzigzagfiltration_simplex_iterator< Simplex_tree > 
                                             Zigzagfiltration_simplex_iterator;
/** Range for the flag zigzag filtration.*/
  typedef boost::iterator_range< Zigzagfiltration_simplex_iterator > 
                                                Zigzagfiltration_simplex_range;
 




//////////////////////to remove
private:
  /* Returns true iff the node in the simplex tree pointed by
   * sh has children. For technical reasons, implemented with a Node and a 
   * Vertex_handle as inputs, instead of a Simplex_handle. */
  bool has_children(Node & node, Vertex_handle u) const {
    return (node.children()->parent() == u);
  }
  Siblings * self_siblings(Node & node, Vertex_handle u) {
    if(node.children()->parent() == u) { return node.children()->oncles(); } 
    else                               { return node.children(); }
  }
 /* Returns the dimension of a simplex.
  * For technical reasons, implemented with a Node and a 
  * Vertex_handle  as inputs, instead of a Simplex_handle. */
  int dimension(Node & node, Vertex_handle u) {
    Siblings * curr_sib = self_siblings(node, u);
    int dim = 0;
    while (curr_sib != NULL) {
      ++dim;
      curr_sib = curr_sib->oncles();
    }
    return dim - 1;
  }
////////////////////////////////////////



public:
  /* Vertex handle siblings data structure. For every vertex u in the complex,
   * apart from the vertices (dim 0), maintains 
   * the list of nodes with label u in vh_siblings_[u]. Used for cofaces location, and/or the local 
   * expansion of the complex when adding the new edge {u,v} in order to find 
   * siblings with both u and v. */
  // typedef boost::intrusive::list < Node
  //                                , boost::intrusive::base_hook< base_hook_st_node_lab >
  //                                , boost::intrusive::constant_time_size<false>
  //                                >                                  Siblings_vector;

  // typedef boost::intrusive::list < Node
  //   , boost::intrusive::member_hook< Node, member_hook_st_node_lab, &Node::hook_ >
  //   , boost::intrusive::constant_time_size<false>
  //   >                                    Siblings_vector;

   /*at index u, list of Siblings containing a node with label u.*/
  typedef std::set< Siblings * >                               Siblings_vector;
  typedef std::vector< Siblings_vector >                        Vh_siblings_ds;
  Vh_siblings_ds                                                  vh_siblings_;
//maintain the list of edges
  std::vector< Edge_type >                    * zigzag_edge_filtration_; 
//bound on the dimension of the complex.
  int                                                          dim_max_; 

 /*
  * A ZigzagEdge type must contain operations:
  * Vertex_handle u() return the endpoint with smaller label,
  * Vertex_handle v() return the endpoint with bigger label,
  * Filtration_value fil() return the filtration value in the zigzag filtration,
  * bool type() return true if the edge is inserted, false if it is removed.
  *
  * Must be Zigzag_edge.
  */
  template<class ZigzagEdge>
  Zigzagfiltration_simplex_range 
  zigzagfiltration_simplex_range( int num_vertices //total number of vertices
                                , std::vector< ZigzagEdge > & zz_edge_fil
                                , int dim_max)
  { //initialize the cofaces data structure, vertices must be labeled from 0 to num_vertices-1
    vh_siblings_.clear(); //vh_siblings_.assign(num_vertices, Siblings_vector()); 
    vh_siblings_.resize(num_vertices);
    for(auto & ref : vh_siblings_) { ref = Siblings_vector(); }
    zigzag_edge_filtration_ = &zz_edge_fil;
    dim_max_ = dim_max;
    return Zigzagfiltration_simplex_range( Zigzagfiltration_simplex_iterator(this)
                                         , Zigzagfiltration_simplex_iterator() );
  }

  Zigzagfiltration_simplex_range 
  zigzagfiltration_simplex_range()
  { 
    return Zigzagfiltration_simplex_range( Zigzagfiltration_simplex_iterator()
                                         , Zigzagfiltration_simplex_iterator() );
  }



public:

 /* 
  * Add an edge in the complex, its two vertices (if missing) 
  * and all the missing 
  * simplices of its star, defined to maintain the property 
  * of the complex to be a flag complex.
  *
  * In term of edges in the graph, inserting edge u,v only affects N^+(u), 
  * even if u and v where missing.
  *
  * For a new node with label v, we first do a local expansion for 
  * computing the 
  * children of this new node, and then a standard expansion for its children. 
  * Nodes with label v (and their subtrees) already in the tree 
  * do not get affected.
  *
  * Nodes with label u get affected only if a Node with label v is in their 
  * siblings set. 
  * We then try to insert "ponctually" v all over the subtree rooted 
  * at Node(u). Each 
  * insertion of a Node with v label induces a local expansion at this 
  * Node (as explained above)
  * and a sequence of "ponctual" insertion of Node(v) in the subtree 
  * rooted at sibling nodes 
  * of the new node, on its left.
  *
  * @param[in] u,v              Vertex_handle representing the new edge
  * @param[in] zz_filtration    Must be empty. Represents at the end the 
  *                             successive simplex insertions agreeing 
  *                             with the 
  *                             zigzag filtration order.
  */
  void zz_add_edge( Vertex_handle                   u 
                  , Vertex_handle                   v
                  , Filtration_value                fil
                  , int                             dim_max
                  , std::vector< Simplex_handle > & zz_filtration ) 
  {
    if(v < u) { std::swap(u,v); } //so as u < v

    static std::vector< std::pair<Siblings *, Vertex_handle> > 
                                                             zz_filtration_tmp;
    //check whether u, v and u,v are in the tree, insert them if necessary
    auto res_ins = root_.members().emplace(v,Node(&root_,fil));
    if(res_ins.second) { zz_filtration_tmp.emplace_back(&root_, v); }

    res_ins = root_.members().emplace(u,Node(&root_,fil)); 
    if(res_ins.second) { zz_filtration_tmp.emplace_back(&root_, u); }

    Simplex_handle sh_u = res_ins.first; //iterator still valid
    if( !has_children(sh_u) ) 
    { sh_u->second.assign_children( new Siblings(&root_, u) ); }
    //check whether the edge {u,v} is already in the complex
    if( sh_u->second.children()->members().find(v) != 
        sh_u->second.children()->members().end() ) 
    { std::cout << "Edge already in complex \n"; return; } 

    //expand the subtree rooted at new node with label v
    zz_punctual_expansion( v                        //insert v in sib
                         , sh_u->second.children()  //sib
                         , fil      //fil value of edge {u,v} in zz filtration
                         , dim_max - 1          //dim_max - dim edges (1)
                         // , u 
                         , zz_filtration_tmp );
    //for all siblings containing u (except root), check if it also contains v
    for( auto sib_u : vh_siblings_[u] )
    { 
      // Siblings * self_sib = self_siblings(node_u, u);
      auto sh_v = sib_u->members().find(v); 
      if(sh_v != sib_u->members().end()) //contains v
      {
        int curr_dim = dimension(sh_v);
        if(curr_dim < dim_max) 
        {
          auto sh_u = sib_u->members().find(u);
          if(!has_children(sh_u->second, u)) { sh_u->second.assign_children(new Siblings(sib_u, u)); } 
          zz_punctual_expansion( v
                               , sh_u->second.children()
                               , fil
                               , dim_max - curr_dim -1 //dim_max - current dim );
                               // , u_d.first 
                               , zz_filtration_tmp ); //u on top         
        }
      }
    }
//the complex is not modified anymore,
// hence the iterators are not invalidated from now on.
//convert zz_filtration_tmp into sequence of Simplex_handle
    for(auto sib_v : zz_filtration_tmp) {
      //todo remove
      if(sib_v.first->members().find(sib_v.second) == sib_v.first->members().end()) 
      { std::cout << "Error here \n"; }
      ////////
      zz_filtration.push_back(sib_v.first->members().find(sib_v.second));
    }
    //sort zz_filtration in filtration order, is_before... is a total ordering
    std::sort( zz_filtration.begin()
             , zz_filtration.end(),   is_before_in_filtration(this));
    zz_filtration_tmp.clear();
 }

private:
/* Insert a Node with label v in the set of siblings sib. sib does not contain 
 * v a priori. As a consequence, the new node with v label has a subtree computed 
 * via zz_local_expansion. All Node (with label x) in sib with a smaller label 
 * may in turn need a local_expansion by v iff Nˆ+(x) contains v.
 */
  void zz_punctual_expansion( Vertex_handle    v
                            , Siblings *       sib
                            , Filtration_value fil
                            , int              k //k == dim_max - dimension simplices in sib
                            , std::vector< std::pair<Siblings *, Vertex_handle> > & zz_filtration_tmp )
  {
    // std::cout << "zz_punctual_expansion \n";
// std::cout << "A1 \n";
    auto res_ins_v = sib->members().emplace(v,Node(sib,fil)); 
// std::cout << "A2 \n";
    if(!res_ins_v.second) 
    { std::cout << "Node v already here in zz_punctual_expansion... \n"; 
      return; } //TODO rm
// std::cout << "A3 " << v << "  " << vh_siblings_.size() << " \n";    
    vh_siblings_[v].insert(sib);
// std::cout << "A4 \n";    
    zz_filtration_tmp.emplace_back(sib, v);

// std::cout << "A \n";
    if(k == 0) { return; }
// std::cout << "B \n";
    //create the subtree of new Node(v)
    zz_local_expansion( res_ins_v.first
                      , sib
                      , fil
                      , k 
                      , zz_filtration_tmp ); 

// std::cout << "C \n";
    //punctual expansion in nodes on the left of v
    for( auto sh = sib->members().begin(); sh != res_ins_v.first; ++sh )
    { //if v belongs to N^+(x), punctual expansion
      Simplex_handle root_sh = find_vertex(sh->first);
      if( has_children(root_sh) && 
          root_sh->second.children()->members().find(v) != root_sh->second.children()->members().end() ) 
      {
        if(!has_children(sh)) { sh->second.assign_children(new Siblings(sib, sh->first)); } 
        zz_punctual_expansion( v
                             , sh->second.children()
                             , fil
                             , k-1
                             , zz_filtration_tmp );
      }
    }

// std::cout << "D \n";
    // std::cout << "done zz_punctual_expansion \n";
  }

/* After the insertion of edge {u,v}, expansion of a subtree rooted at v, where the 
 * Node with label v has just been inserted. sh has no children here. No need to 
 * add things in the cofaces of u,v. */
  void zz_local_expansion( Simplex_handle   sh_v     //Node with label v which has just been inserted
                         , Siblings       * curr_sib //Siblings containing the node
                         , Filtration_value fil_uv   //Filtration value of the edge uv in the zz filtration
                         , int              k        //Stopping condition for recursion based on dimension maximal
                         // , Vertex_handle    u        //on top of the current tree
                         , std::vector< std::pair<Siblings *, Vertex_handle > > & zz_filtration_tmp ) //<-depth first insertion
  {
    // std::cout << "zz_local_expansion \n";
    //if (dimension_ > k) { dimension_ = k; }    // TODO
    if (k == 0)           { return; }             

    Simplex_handle next_it = sh_v;    ++next_it; 
    static std::vector< std::pair<Vertex_handle, Node> > inter;  // static, not thread-safe.
    Simplex_handle root_sh_v = find_vertex(sh_v->first);
    if(!has_children(root_sh_v)) { return; } 

    zz_intersection( inter
                   , next_it
                   , curr_sib->members().end()
                   , root_sh_v->second.children()->members().begin()
                   , root_sh_v->second.children()->members().end()
                   , fil_uv );
    
    if(!inter.empty()) 
    {
      // this->num_simplices_ += inter.size();
      Siblings * new_sib = new Siblings(curr_sib, sh_v->first, inter);

      //link in cofaces and vh ds
      for( auto new_sh = new_sib->members().begin(); 
           new_sh != new_sib->members().end(); ++new_sh )
      { 
        // new_sh->second.assign_children(new_sib); //done in siblings.h
        vh_siblings_[new_sh->first].insert(new_sib); //new in vh ds
        zz_filtration_tmp.emplace_back(new_sib, new_sh->first);
      }
      sh_v->second.assign_children(new_sib);
      inter.clear();

      zz_siblings_expansion(new_sib, fil_uv, k-1, zz_filtration_tmp ); //recursive standard expansion for the rest of the subtree
    }
    else { sh_v->second.assign_children(curr_sib); inter.clear(); }
  }


  //TODO boost::container::ordered_unique_range_t in the creation of a Siblings
  // we also need to insert all nodes in the cofaces ds and the siblings_uv.
  // copy constructor? what about the hooks then?

  /* Global expansion of a subtree in the simplex tree.
   *
   * The filtration value is absolute and defined by "Filtration_value fil". 
   * The new Node are also connected appropriately in the cofaces and siblings_uv
   * data structures.
   */
  void zz_siblings_expansion( Siblings       * siblings  // must contain elements
                            , Filtration_value fil
                            , int              k  //==max_dim expansion - dimension curr siblings 
                            // , Vertex_handle    u  //top of the current tree
                            , std::vector< std::pair<Siblings *, Vertex_handle> > & zz_filtration_tmp )
  {
    // if (dimension_ > k) { dimension_ = k; } //TODO  
    if (k == 0)            { return; }            

    Dictionary_it next = ++(siblings->members().begin());

    static std::vector< std::pair<Vertex_handle, Node> > inter;  // static, not thread-safe.
    for( Dictionary_it s_h = siblings->members().begin();
         next != siblings->members().end(); ++s_h, ++next) 
    {
      Simplex_handle root_sh = find_vertex(s_h->first);
      if( has_children(root_sh) )
      {
        zz_intersection( inter                      // output intersection
                       , next                       // begin
                       , siblings->members().end()  // end
                       , root_sh->second.children()->members().begin()
                       , root_sh->second.children()->members().end()
                       , fil   );

        if ( !inter.empty() ) 
        {
          // this->num_simplices_ += inter.size();
          Siblings * new_sib = new Siblings( siblings    // oncles
                                           , s_h->first  // parent
                                           , inter);     // boost::container::ordered_unique_range_t

          for( auto new_sh = new_sib->members().begin(); 
               new_sh != new_sib->members().end(); ++new_sh )
          { 
            new_sh->second.assign_children(new_sib);
            vh_siblings_[new_sh->first].insert(new_sib); //new in vh ds
            zz_filtration_tmp.emplace_back(new_sib, new_sh->first);
          }

          inter.clear();
          s_h->second.assign_children(new_sib);
          zz_siblings_expansion(new_sib, fil, k - 1, zz_filtration_tmp);
        }
        else { s_h->second.assign_children(siblings); inter.clear();}  // ensure the children property
      }
    }
  }
  /* \brief Intersects Dictionary 1 [begin1;end1) with Dictionary 2 [begin2,end2)
   * and assigns Filtration_value fil to the Nodes. */
  void zz_intersection( std::vector<std::pair<Vertex_handle, Node> > & intersection
                      , Dictionary_it                                  begin1
                      , Dictionary_it                                  end1
                      , Dictionary_it                                  begin2
                      , Dictionary_it                                  end2
                      , Filtration_value                               fil ) 
  {
    if (begin1 == end1 || begin2 == end2) { return; } 
    while (true) 
    {
      if (begin1->first < begin2->first) 
      { ++begin1; if(begin1 == end1) { return; } }
      else 
      {
        if (begin1->first > begin2->first) 
        { ++begin2; if(begin2 == end2) { return; } }
        else // begin1->first == begin2->first
        {
          intersection.emplace_back( begin1->first, Node( NULL, fil ) ); //TODO NULL ? or self siblings?
          ++begin1; ++begin2;
          if (begin1 == end1 || begin2 == end2) { return; }
        }
      }
    }
  }

//TODO description -> adapted for lazy_remove_edge
  bool reverse_reverse_lexicographic_order(Simplex_handle sh1, Simplex_handle sh2) {
    Simplex_vertex_range rg1 = simplex_vertex_range(sh1);
    Simplex_vertex_range rg2 = simplex_vertex_range(sh2);
    Simplex_vertex_iterator it1 = rg1.begin();
    Simplex_vertex_iterator it2 = rg2.begin();
    while (it1 != rg1.end() && it2 != rg2.end()) {
      if (*it1 == *it2) { ++it1; ++it2; } 
      else { return *it1 < *it2; }//they are different
    }
    return !((it1 != rg1.end()) && (it2 == rg2.end()));//is negated afterward
  }

  /* StrictWeakOrdering, for the simplices, defined by the filtration.
   *
   * It corresponds to the partial order
   * induced by the filtration values, with ties resolved using reverse lexicographic order.
   * Reverse lexicographic order has the property to always consider the subsimplex of a simplex
   * to be smaller. The filtration function must be monotonic. */
  struct is_after_in_filtration 
  {
    explicit is_after_in_filtration(Simplex_tree * st)
        : st_(st) {}

    bool operator()(const Simplex_handle sh1, const Simplex_handle sh2) const {
      //no consideration of filtration value when removing
      return !(st_->reverse_reverse_lexicographic_order(sh1, sh2));  // is sh1 a proper subface of sh2
    }

    Simplex_tree * st_;
  };

/* Given a simplex handle in a simplex tree cpx_, traverse the tree upwards 
 * to find vertex u_. Does not test sh itself. Used for filter_iterator when 
 * traversing vh_siblings_.
 */
 template<class SimplexTree>
 class does_simplex_contain_vertex {
    typedef typename SimplexTree::Vertex_handle Vertex_handle;
    typedef typename SimplexTree::Simplex_handle Simplex_handle;

  public:
    does_simplex_contain_vertex() : cpx_(NULL), u_(-1), v_(-1) {}

    does_simplex_contain_vertex(SimplexTree* cpx, Vertex_handle u, Vertex_handle v) 
    : cpx_(cpx), u_(u), v_(v) {}

//node must contain v_ as label
    bool operator()(typename SimplexTree::Siblings *curr_sib) {
      // auto curr_sib = cpx_->self_siblings(node,v_);
      while(curr_sib->oncles() != NULL) {
        if(curr_sib->parent() == u_) {return true;}
        curr_sib = curr_sib->oncles();
      }
      return false;
    }
  private:
    SimplexTree *cpx_;
    Vertex_handle u_;
    Vertex_handle v_;
  };


//TODO heuristic for ordering the removal by most recent insertion?
public:
  /* Computes, in a filtration order, all simplices that ought to be removed 
  * if the edge {u,v} were to disappear (puts them in zz_filtration). 
  * This method does NOT modify the simplex tree. */
  void zz_lazy_remove_edge( Vertex_handle u
                          , Vertex_handle v
                          , std::vector< Simplex_handle > & zz_filtration ) 
  {
    // std::cout << "lazy rm edge \n";

    if(v < u) { std::swap(u,v); } //so as u < v
    typedef does_simplex_contain_vertex< Simplex_tree >     cof_predicate;
    typedef boost::filter_iterator< cof_predicate
                                  , typename Siblings_vector::iterator > Edge_coface_iterator;
    cof_predicate cof_pred(this,u,v);
    Edge_coface_iterator beg_it = Edge_coface_iterator( cof_pred
                                                      , vh_siblings_[v].begin()
                                                      , vh_siblings_[v].end() );
    Edge_coface_iterator end_it = Edge_coface_iterator( cof_pred
                                                      , vh_siblings_[v].end()
                                                      , vh_siblings_[v].end() );    
    for(auto sib_it = beg_it; sib_it != end_it; ++sib_it) //all siblings that contain a coface of {u,v} with v as max vertex.
    {
      auto sh_v = (*sib_it)->members().find(v); 
      if(has_children(sh_v))
      { zz_lazy_remove_edge_depth_first_traversal(sh_v->second.children(), zz_filtration); }
      zz_filtration.push_back(sh_v);
    }

    std::sort(zz_filtration.begin(), zz_filtration.end(), is_after_in_filtration(this));
  }

  void zz_lazy_empty_complex(std::vector< Simplex_handle > & zz_filtration)
  {
    for(auto sh : complex_simplex_range()) { zz_filtration.push_back(sh); }
    std::sort(zz_filtration.begin(), zz_filtration.end(), is_after_in_filtration(this));
  }


private:
  /* Reads recursively the tree in a depth first fashion. Non terminal recursion.*/
  void zz_lazy_remove_edge_depth_first_traversal( Siblings * sib
                                                , std::vector< Simplex_handle > & zz_filtration )
  {
    for(auto sh = sib->members().begin(); sh != sib->members().end(); ++sh)
    {
      if(has_children(sh)) 
      { zz_lazy_remove_edge_depth_first_traversal(sh->second.children(), zz_filtration); }
      zz_filtration.push_back(sh);
    }
  }


/***********************************************/
/****************** End ************************/
/***********************************************/




 private:
  Vertex_handle null_vertex_;
  /** \brief Upper bound on the filtration values of the simplices.*/
  Filtration_value threshold_;
  /** \brief Total number of simplices in the complex, without the empty simplex.*/
  /** \brief Set of simplex tree Nodes representing the vertices.*/
  Siblings root_;
  /** \brief Simplices ordered according to a filtration.*/
  std::vector<Simplex_handle> filtration_vect_;
  /** \brief Upper bound on the dimension of the simplicial complex.*/
  int dimension_;
};

// Print a Simplex_tree in os.
template<typename...T>
std::ostream& operator<<(std::ostream & os, Simplex_tree<T...> & st) {
  for (auto sh : st.filtration_simplex_range()) {
    os << st.dimension(sh) << " ";
    for (auto v : st.simplex_vertex_range(sh)) {
      os << v << " ";
    }
    os << st.filtration(sh) << "\n";  // TODO(VR): why adding the key ?? not read ?? << "     " << st.key(sh) << " \n";
  }
  return os;
}

template<typename...T>
std::istream& operator>>(std::istream & is, Simplex_tree<T...> & st) {
  typedef Simplex_tree<T...> ST;
  std::vector<typename ST::Vertex_handle> simplex;
  typename ST::Filtration_value fil;
  typename ST::Filtration_value max_fil = 0;
  int max_dim = -1;
  while (read_simplex(is, simplex, fil)) {
    // read all simplices in the file as a list of vertices
    // Warning : simplex_size needs to be casted in int - Can be 0
    int dim = static_cast<int> (simplex.size() - 1);
    if (max_dim < dim) {
      max_dim = dim;
    }
    if (max_fil < fil) {
      max_fil = fil;
    }
    // insert every simplex in the simplex tree
    st.insert_simplex(simplex, fil);
    simplex.clear();
  }
  st.set_dimension(max_dim);
  st.set_filtration(max_fil);

  return is;
}

/// Model of SimplexTreeOptions.
struct Simplex_tree_options_full_featured {
  typedef linear_indexing_tag Indexing_tag;
  typedef int Vertex_handle;
  typedef double Filtration_value;
  typedef int Simplex_key;
  static const bool store_key = true;
  static const bool store_filtration = true;
  static const bool contiguous_vertices = false;
  static const bool link_simplices_through_max_vertex = false;
};

/** Model of SimplexTreeOptions, faster than
  `Simplex_tree_options_full_featured` but note the unsafe
  `contiguous_vertices` option. */
struct Simplex_tree_options_fast_persistence {
  typedef linear_indexing_tag Indexing_tag;
  typedef int Vertex_handle;
  typedef float Filtration_value;
  typedef int Simplex_key;
  static const bool store_key = true;
  static const bool store_filtration = true;
  static const bool contiguous_vertices = true;
  static const bool link_simplices_through_max_vertex = false;
};

/** Model of SimplexTreeOptions, with a zigzag_indexing_tag.
    Note the unsafe `contiguous_vertices` option. */
struct Simplex_tree_options_zigzag_persistence {
  typedef zigzag_indexing_tag Indexing_tag;
  typedef int Vertex_handle;
  typedef float Filtration_value;
  typedef int Simplex_key;
  static const bool store_key = true;
  static const bool store_filtration = true;
  static const bool contiguous_vertices = true;
  static const bool link_simplices_through_max_vertex = true;
};

/** @} */  // end defgroup simplex_tree

}  // namespace Gudhi

#endif  // SIMPLEX_TREE_H_
