/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Clement Jamin
 *
 *    Copyright (C) 2017 INRIA
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

#ifndef PERIODIC_KD_TREE_SEARCH_H_
#define PERIODIC_KD_TREE_SEARCH_H_

#include <gudhi/Periodic_Euclidean_distance.h>

#include <CGAL/K_neighbor_search.h>
#include <CGAL/Incremental_neighbor_search.h>
#include <CGAL/Search_traits.h>
#include <CGAL/Search_traits_adapter.h>
#include <CGAL/Fuzzy_sphere.h>
#include <CGAL/property_map.h>

#include <boost/property_map/property_map.hpp>
#include <boost/iterator/counting_iterator.hpp>
#include <boost/mpl/if.hpp>

#include <cstddef>
#include <vector>

namespace Gudhi {
namespace spatial_searching {


  /**
  * \class Periodic_kd_tree_search Periodic_kd_tree_search.h gudhi/Periodic_kd_tree_search.h
  * \brief Spatial tree data structure to perform periodic (approximate) nearest neighbor search.
  *
  * \ingroup spatial_searching
  *
  * \details
  * The class Periodic_kd_tree_search is a tree-based data structure, based on
  * <a target="_blank" href="http://doc.cgal.org/latest/Spatial_searching/index.html">CGAL dD spatial searching data structures</a>.
  * It provides a simplified API to perform periodic (approximate) nearest neighbor searches. Contrary to CGAL default behavior, the tree
  * does not store the points themselves, but stores indices.
  *
  * There are two types of queries: the <i>k-nearest neighbor query</i>, where <i>k</i> is fixed and the <i>k</i> nearest 
  * points are computed right away,
  * and the <i>incremental nearest neighbor query</i>, where no number of neighbors is provided during the call, as the
  * neighbors will be computed incrementally when the iterator on the range is incremented.
  *
  * \tparam Search_traits must be a model of the <a target="_blank"
  *   href="http://doc.cgal.org/latest/Spatial_searching/classSearchTraits.html">SearchTraits</a>
  *   concept, such as the <a target="_blank"
  *   href="http://doc.cgal.org/latest/Kernel_d/classCGAL_1_1Epick__d.html">CGAL::Epick_d</a> class, which
  *   can be static if you know the ambiant dimension at compile-time, or dynamic if you don't.
  * \tparam Point_range is the type of the range that provides the points.
  *   It must be a range whose iterator type is a `RandomAccessIterator`.
  * \tparam Split_strategy allows to choose between different splitting strategies: 
  * `Periodic_splitter_enum::SLIDING_MIDPOINT` (default), `Periodic_splitter_enum::MEDIAN_OF_MAX_SPREAD`, 
  * `Periodic_splitter_enum::MIDPOINT_OF_MAX_SPREAD`.
  */

enum class Periodic_splitter_enum { SLIDING_MIDPOINT, MEDIAN_OF_MAX_SPREAD, MIDPOINT_OF_MAX_SPREAD };

template <typename Search_traits, typename Point_range, Periodic_splitter_enum Split_strategy = Periodic_splitter_enum::SLIDING_MIDPOINT>
class Periodic_kd_tree_search {
  typedef boost::iterator_property_map<
    typename Point_range::const_iterator,
    CGAL::Identity_property_map<std::ptrdiff_t> >           Point_property_map;

public:
  /// The Traits.
  typedef Search_traits                                     Traits;
  /// Number type used for distances.
  typedef typename Traits::FT                               FT;
  /// The point type.
  typedef typename Point_range::value_type                  Point;

  typedef CGAL::Search_traits_adapter<
    std::ptrdiff_t,
    Point_property_map,
    Search_traits>                                          STraits;
  typedef Periodic_Euclidean_distance<Search_traits>        Base_distance;
  typedef CGAL::Distance_adapter<
    std::ptrdiff_t,
    Point_property_map,
    Base_distance >                                         Periodic_distance;
  
  typedef typename boost::mpl::if_c <
    Split_strategy == Periodic_splitter_enum::MEDIAN_OF_MAX_SPREAD,
    CGAL::Median_of_max_spread<STraits>,
    typename boost::mpl::if_c <
      Split_strategy == Periodic_splitter_enum::MIDPOINT_OF_MAX_SPREAD,
      CGAL::Midpoint_of_max_spread<STraits>,
      CGAL::Sliding_midpoint<STraits> // default in CGAL -- best choice in average
    >::type
  >::type                                                   Splitter;

  typedef CGAL::Kd_tree<
    STraits, Splitter, CGAL::Tag_true>      Tree;

  typedef CGAL::K_neighbor_search<
    STraits, Periodic_distance, Splitter, Tree>           K_neighbor_search;
  /// \brief The range returned by a k-nearest neighbor search.
  /// Its value type is `std::pair<std::size_t, FT>` where `first` is the index
  /// of a point P and `second` is the squared distance between P and the query point.
  typedef K_neighbor_search                                 KNS_range;

  typedef CGAL::Incremental_neighbor_search<
    STraits, Periodic_distance, Splitter, Tree>           Incremental_neighbor_search;
  /// \brief The range returned by an incremental nearest neighbor search.
  /// Its value type is `std::pair<std::size_t, FT>` where `first` is the index
  /// of a point P and `second` is the squared distance between P and the query point.
  typedef Incremental_neighbor_search                       INS_range;

  typedef CGAL::Fuzzy_sphere<STraits>                       Fuzzy_sphere;

  /// \brief Constructor
  /// @param[in] points Const reference to the point range. This range
  /// is not copied, so it should not be destroyed or modified afterwards.
  /// @param[in] domain The periodic domain.
  Periodic_kd_tree_search(Point_range const& points, typename Base_distance::Iso_box_d const& domain)
  : m_points(points),
    m_tree(boost::counting_iterator<std::ptrdiff_t>(0),
           boost::counting_iterator<std::ptrdiff_t>(points.size()),
           Splitter(),
           STraits(std::begin(points))),
    m_domain(domain) {
    // Build the tree now (we don't want to wait for the first query)
    m_tree.build();
  }

  /// \brief Constructor
  /// @param[in] points Const reference to the point range. This range
  /// is not copied, so it should not be destroyed or modified afterwards.
  /// @param[in] domain The periodic domain.
  /// @param[in] only_these_points Specifies the indices of the points that
  /// should be actually inserted into the tree. The other points are ignored.
  template <typename Point_indices_range>
  Periodic_kd_tree_search(
    Point_range const& points,
    Point_indices_range const& only_these_points,
    typename Base_distance::Iso_box_d const& domain)
    : m_points(points),
      m_tree(
        only_these_points.begin(), only_these_points.end(),
        Splitter(),
        STraits(std::begin(points))),
      m_domain(domain) {
    // Build the tree now (we don't want to wait for the first query)
    m_tree.build();
  }

  /// \brief Constructor
  /// @param[in] points Const reference to the point range. This range
  /// is not copied, so it should not be destroyed or modified afterwards.
  /// @param[in] begin_idx, past_the_end_idx Define the subset of the points that
  /// should be actually inserted into the tree. The other points are ignored.
  Periodic_kd_tree_search(
    Point_range const& points,
    std::size_t begin_idx, std::size_t past_the_end_idx,
    typename Base_distance::Iso_box_d const& domain)
  : m_points(points),
    m_tree(
      boost::counting_iterator<std::ptrdiff_t>(begin_idx),
      boost::counting_iterator<std::ptrdiff_t>(past_the_end_idx),
      Splitter(),
      STraits(std::begin(points))),
    m_domain(domain) {
    // Build the tree now (we don't want to wait for the first query)
    m_tree.build();
  }

  // Be careful, this function invalidates the tree,
  // which will be recomputed at the next query
  void insert(std::ptrdiff_t point_idx) {
    m_tree.insert(point_idx);
  }

  /// \brief Search for the k-nearest neighbors from a query point.
  /// @param[in] p The query point.
  /// @param[in] k Number of nearest points to search.
  /// @param[in] sorted Indicates if the computed sequence of k-nearest neighbors needs to be sorted.
  /// @param[in] eps Approximation factor.
  /// @return A range (whose `value_type` is `std::size_t`) containing the k-nearest neighbors.
  KNS_range query_k_nearest_neighbors(
    Point const& p,
    unsigned int k,
    bool sorted = true,
    FT eps = FT(0)) const {
    // Initialize the search structure, and search all N points
    // Note that we need to pass the Distance explicitly since it needs to
    // know the property map
    K_neighbor_search search(
      m_tree,
      p,
      k,
      eps,
      true,
      Periodic_distance(std::begin(m_points), Base_distance(m_domain, m_tree.traits())),
      sorted);

    return search;
  }

  /// \brief Search incrementally for the nearest neighbors from a query point.
  /// @param[in] p The query point.
  /// @param[in] eps Approximation factor.
  /// @return A range (whose `value_type` is `std::size_t`) containing the 
  /// neighbors sorted by their distance to p.
  /// All the neighbors are not computed by this function, but they will be
  /// computed incrementally when the iterator on the range is incremented.
  INS_range query_incremental_nearest_neighbors(Point const& p, FT eps = FT(0)) const {
    // Initialize the search structure, and search all N points
    // Note that we need to pass the Distance explicitly since it needs to
    // know the property map
    Incremental_neighbor_search search(
      m_tree,
      p,
      eps,
      true,
      Periodic_distance(std::begin(m_points), Base_distance(m_domain, m_tree.traits())) );

    return search;
  }

  /// \brief Search for all the neighbors in a ball.
  /// @param[in] p The query point.
  /// @param[in] radius The search radius
  /// @param[out] it The points that lie inside the sphere of center `p` and radius `radius`.
  ///                Note: `it` is used this way: `*it++ = each_point`.
  /// @param[in] eps Approximation factor.
  // TODO: needs a Periodic_fuzzy_sphere
  /*template <typename OutputIterator>
  void near_search(
    Point const& p,
    FT radius,
    OutputIterator it,
    FT eps = FT(0)) const {
    
    m_tree.search(it, Fuzzy_sphere(p, radius, eps, m_tree.traits()));
  }*/

  /// \brief Search for any neighbor in a ball.
  /// @param[in] p The query point.
  /// @param[in] radius The search radius
  /// @param[in] eps Approximation factor.
  /// @return The index of a point approximately contained by the sphere of center `p`
  ///         and radius `radius`, or -1 if no point could be found.
  // TODO: needs a Periodic_fuzzy_sphere
  /*std::ptrdiff_t any_near_neighbor(
    Point const& p,
    FT radius,
    FT eps = FT(0)) const {

    auto ret = m_tree.search_any_point(Fuzzy_sphere(p, radius, eps, m_tree.traits()));
    return (ret ? *ret : -1);
  }*/

  int tree_depth() const
  {
    return m_tree.root()->depth();
  }

 private:
  Point_range const& m_points;
  Tree m_tree;
  typename Base_distance::Iso_box_d m_domain;
};

}  // namespace spatial_searching
}  // namespace Gudhi

#endif  // PERIODIC_KD_TREE_SEARCH_H_
