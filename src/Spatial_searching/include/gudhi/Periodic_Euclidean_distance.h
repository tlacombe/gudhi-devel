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


#ifndef PERIODIC_EUCLIDEAN_DISTANCE_H_
#define PERIODIC_EUCLIDEAN_DISTANCE_H_

#include <CGAL/Kd_tree_rectangle.h>
#include <CGAL/number_utils.h>
#include <CGAL/internal/Get_dimension_tag.h>
#include <vector>

namespace Gudhi {
namespace spatial_searching {

  template <class SearchTraits>
  class Periodic_Euclidean_distance {
  public:

    typedef typename SearchTraits::FT                                   FT;
    typedef typename SearchTraits::Point_d                              Point_d;
    typedef typename SearchTraits::Iso_box_d                            Iso_box_d;
    typedef typename SearchTraits::Construct_min_vertex_d               Construct_min_vertex_d;
    typedef typename SearchTraits::Construct_max_vertex_d               Construct_max_vertex_d;
    typedef typename SearchTraits::Construct_cartesian_const_iterator_d Construct_cartesian_const_iterator_d;
    typedef typename SearchTraits::Cartesian_const_iterator_d           Cartesian_const_iterator_d;
    typedef Point_d Query_item;

    typedef typename CGAL::internal::Get_dimension_tag<SearchTraits>::Dimension D;

    // default constructor
    Periodic_Euclidean_distance(
      const Iso_box_d &domain_,
      const SearchTraits& traits_ = SearchTraits())
    : m_traits(traits_)
    , m_construct_it(m_traits.construct_cartesian_const_iterator_d_object())
    , m_domain(domain_)
    {
      m_min_vertex = m_traits.construct_min_vertex_d_object()(m_domain);
      m_max_vertex = m_traits.construct_max_vertex_d_object()(m_domain);
      m_min_begin = m_construct_it(m_min_vertex);
      m_max_begin = m_construct_it(m_max_vertex);
    }

    inline FT transformed_distance(const Query_item& q, const Point_d& p) const
    {
      Cartesian_const_iterator_d p_begin = m_construct_it(p), p_end = m_construct_it(p, 0);
      return transformed_distance_from_coordinates(q, p_begin, p_end);
    }

    template <typename Coord_iterator>
    inline FT transformed_distance_from_coordinates(const Query_item& q,
      Coord_iterator it_coord_begin, Coord_iterator it_coord_end) const
    {
      return transformed_distance_from_coordinates(q, it_coord_begin, it_coord_end, D());
    }

    // Static dim = 2 loop unrolled
    template <typename Coord_iterator>
    inline FT transformed_distance_from_coordinates(const Query_item& q,
      Coord_iterator it_coord_begin, Coord_iterator /*unused*/,
      CGAL::Dimension_tag<2>) const
    {
      Cartesian_const_iterator_d qit = m_construct_it(q);
      FT distance = CGAL::square(periodic_one_coord_distance(0, *qit, *it_coord_begin));
      qit++; it_coord_begin++;
      distance += CGAL::square(periodic_one_coord_distance(1, *qit, *it_coord_begin));
      return distance;
    }

    // Static dim = 3 loop unrolled
    template <typename Coord_iterator>
    inline FT transformed_distance_from_coordinates(const Query_item& q,
      Coord_iterator it_coord_begin, Coord_iterator /*unused*/,
      CGAL::Dimension_tag<3>) const
    {
      Cartesian_const_iterator_d qit = m_construct_it(q);
      FT distance = CGAL::square(periodic_one_coord_distance(0, *qit, *it_coord_begin));
      qit++; it_coord_begin++;
      distance += CGAL::square(periodic_one_coord_distance(1, *qit, *it_coord_begin));
      qit++; it_coord_begin++;
      distance += CGAL::square(periodic_one_coord_distance(2, *qit, *it_coord_begin));
      return distance;
    }

    // Other cases: static dim > 3 or dynamic dim
    template <typename Coord_iterator, typename Dim>
    inline FT transformed_distance_from_coordinates(const Query_item& q,
      Coord_iterator it_coord_begin, Coord_iterator /*unused*/,
      Dim) const
    {
      FT distance = FT(0);
      Cartesian_const_iterator_d qit = m_construct_it(q), qe = m_construct_it(q, 1);
      for (int i = 0; qit != qe; ++qit, ++it_coord_begin, ++i)
      {
        FT diff = periodic_one_coord_distance(i, *qit, *it_coord_begin);
        distance += diff*diff;
      }
      return distance;
    }

    // During the computation, if the partially-computed distance `pcd` gets greater or equal
    // to `stop_if_geq_to_this`, the computation is stopped and `pcd` is returned
    template <typename Coord_iterator>
    inline FT interruptible_transformed_distance(const Query_item& q,
      Coord_iterator it_coord_begin, Coord_iterator /*unused*/,
      FT stop_if_geq_to_this) const
    {
      FT distance = FT(0);
      Cartesian_const_iterator_d qit = m_construct_it(q), qe = m_construct_it(q, 1);
      int dim = -1;
      if (qe - qit >= 6)
      {
        // Every 4 coordinates, the current partially-computed distance
        // is compared to stop_if_geq_to_this
        // Note: the concept SearchTraits specifies that Cartesian_const_iterator_d 
        //       must be a random-access iterator
        Cartesian_const_iterator_d qe_minus_5 = qe - 5;
        for (;;)
        {
          FT diff = periodic_one_coord_distance(++dim, *qit, *it_coord_begin);
          distance += diff*diff;
          ++qit; ++it_coord_begin;
          diff = periodic_one_coord_distance(++dim, *qit, *it_coord_begin);
          distance += diff*diff;
          ++qit; ++it_coord_begin;
          diff = periodic_one_coord_distance(++dim, *qit, *it_coord_begin);
          distance += diff*diff;
          ++qit; ++it_coord_begin;
          diff = periodic_one_coord_distance(++dim, *qit, *it_coord_begin);
          distance += diff*diff;
          ++qit, ++it_coord_begin;

          if (distance >= stop_if_geq_to_this)
            return distance;

          if (qit >= qe_minus_5)
            break;
        }
      }
      for (; qit != qe; ++qit, ++it_coord_begin)
      {
        FT diff = periodic_one_coord_distance(++dim, *qit, *it_coord_begin);
        distance += diff*diff;
      }
      return distance;
    }

    inline FT min_distance_to_rectangle(const Query_item& q,
      const CGAL::Kd_tree_rectangle<FT, D>& r) const {
      FT distance = FT(0);
      Cartesian_const_iterator_d qit = m_construct_it(q),
        qe = m_construct_it(q, 1);
      for (unsigned int i = 0; qit != qe; i++, qit++) {
        distance += CGAL::square(min_distance_to_periodic_interval(i, *qit, r.min_coord(i), r.max_coord(i)));
      }
      return distance;
    }

    inline FT min_distance_to_rectangle(const Query_item& q,
      const CGAL::Kd_tree_rectangle<FT, D>& r, std::vector<FT>& dists) const {
      FT distance = FT(0);
      Cartesian_const_iterator_d qit = m_construct_it(q),
        qe = m_construct_it(q, 1);
      for (unsigned int i = 0; qit != qe; i++, qit++) {
        FT dx = min_distance_to_periodic_interval(i, *qit, r.min_coord(i), r.max_coord(i));
        dists[i] = dx;
        distance += CGAL::square(dx);
      }
      return distance;
    }

    // Does not make sense for periodic spaces
    inline FT max_distance_to_rectangle(const Query_item& q,
      const CGAL::Kd_tree_rectangle<FT, D>& r) const {
      return 0.;
    }

    // Does not make sense for periodic spaces
    inline FT max_distance_to_rectangle(const Query_item& q,
      const CGAL::Kd_tree_rectangle<FT, D>& r, std::vector<FT>& dists) const {
      return 0.;
    }

    inline FT transformed_distance(FT d) const {
      return d*d;
    }

    inline FT inverse_of_transformed_distance(FT d) const {
      return CGAL::sqrt(d);
    }

private:
  // Warning: the returned "distance" might be negative
  inline FT periodic_one_coord_distance(int dim, FT x_i, FT x_j) const
  {
    FT dx = x_j - x_i;
    FT min_x = m_min_begin[dim];
    FT max_x = m_max_begin[dim];
    FT size_x = max_x - min_x;
    if (dx > 0.5 * size_x)
      dx -= size_x;
    else if (dx <= -0.5 * size_x)
      dx += size_x;

    return dx;
  }

  inline FT min_distance_to_periodic_interval(int dim, FT c, FT low, FT hi) const
  {
    FT min_x = m_min_begin[dim];
    FT max_x = m_max_begin[dim];
    FT size_x = max_x - min_x;

    if (c < low)
      // Check if the interval on the left is not closer
      return (std::min)(low - c, c - (hi - size_x));
    else if (c > high)
      // Check if the interval on the right is not closer
      return (std::min)(c - hi, low + size_x - c);
    else
      // c is inside the interval
      return 0.;
  }

  inline FT max_distance_to_periodic_interval(int dim, FT c, FT low, FT hi) const
  {
    FT min_x = m_min_begin[dim];
    FT max_x = m_max_begin[dim];
    FT size_x = max_x - min_x;
    FT dist_to_low = periodic_one_coord_distance(dim, c, low);
    FT dist_to_hi = periodic_one_coord_distance(dim, c, hi);

    if (dist_to_low < 0 && dist_to_hi > 0)
      return size_x/2;

    return (std::max)(std::abs(dist_to_low), std::abs(dist_to_hi));
  }

  SearchTraits                          m_traits;
  Construct_cartesian_const_iterator_d  m_construct_it;
  Iso_box_d                             m_domain;
  Point_d                               m_min_vertex;
  Point_d                               m_max_vertex;
  Cartesian_const_iterator_d            m_min_begin;
  Cartesian_const_iterator_d            m_max_begin;
}; // class Periodic_Euclidean_distance


}  // namespace spatial_searching
}  // namespace Gudhi

#endif // PERIODIC_EUCLIDEAN_DISTANCE_H_
