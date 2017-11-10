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

#ifndef FUNCTOR_GUDHI_PERIODIC_KD_TREE_SEARCH_
#define FUNCTOR_GUDHI_PERIODIC_KD_TREE_SEARCH_

#include <gudhi/Periodic_kd_tree_search.h>
#include <CGAL/Epick_d.h>

#include <vector>
#include <utility>

namespace gss = Gudhi::spatial_searching;

template <typename Kernel, gss::Periodic_splitter_enum Splitting_strategy = gss::Periodic_splitter_enum::SLIDING_MIDPOINT>
class GUDHI_Periodic_kd_tree_search
{
  typedef Kernel                                      K;
  typedef typename K::Point_d                         Point;
  typedef std::vector<Point>                          Points;

  typedef gss::Periodic_kd_tree_search<K, Points, Splitting_strategy> Points_ds;

public:
  GUDHI_Periodic_kd_tree_search(Points const& points, double /*epsilon*/, typename Kernel::Iso_box_d domain)
    : m_tree(points, domain)
    , m_k(K().point_dimension_d_object()(*points.begin()))
    , m_original_points(points)
  {}

  // Returns the sum of the indices
  std::size_t query_k_nearest_neighbors(
    Point const& p,
    unsigned int k,
    double eps = 0.,
    std::vector<std::pair<std::size_t, double>> *result = NULL) const
  {
    Points_ds::KNS_range neighbors = m_tree.query_k_nearest_neighbors(p, k, true, eps);

#ifdef PRINT_FOUND_NEIGHBORS
    std::cerr << "Query: " << p[0] << "; " << p[1] << "\n";
    for (auto nb : neighbors)
      std::cerr << "  " << nb.first << " : " << nb.second << "\n";
#endif

    std::size_t sum = 0;
    if (result) {
      for (auto nb : neighbors)
      {
        sum += nb.first;
        result->push_back(std::make_pair(nb.first, nb.second));
      }
    }
    else {
      for (auto nb : neighbors)
        sum += nb.first;
    }
    return sum;
  }

  // Returns the sum of the indices
  /*std::size_t query_all_near_neighbors(
    Point const& p,
    double radius,
    double eps = 0.,
    std::vector<std::pair<std::size_t, double>> *result = NULL) const
  {
    std::vector<std::size_t> neighbors;
    m_tree.near_search(p, radius, std::back_inserter(neighbors), eps);

#ifdef PRINT_FOUND_NEIGHBORS
    std::cerr << "Query:\n";
    for (auto nb : neighbors) {
      double sqdist = m_k.squared_distance_d_object()(m_original_points[nb], p);
      std::cerr << "  " << nb << " : " << sqdist << "\n";
    }
#endif

    std::size_t sum = 0;
    if (result) {
      for (auto nb : neighbors)
      {
        sum += nb;
        double sqdist = m_k.squared_distance_d_object()(m_original_points[nb], p);
        result->push_back(std::make_pair(nb, sqdist));
      }
    }
    else {
      for (auto nb : neighbors)
        sum += nb;
    }

    return sum;
  }

  std::ptrdiff_t query_any_near_neighbor(
    Point const& p,
    double radius,
    double eps = 0.,
    std::pair<std::ptrdiff_t, double> *result = NULL) const
  {
    std::ptrdiff_t ret = m_tree.any_near_neighbor(p, radius, eps);

#ifdef PRINT_FOUND_NEIGHBORS
    if (ret == -1)
      std::cerr << "Any near neighbor: none\n";
    else
      std::cerr << "Any near neighbor: " << ret << ": "
        << m_k.squared_distance_d_object()(m_original_points[ret], p) << "\n";
#endif


    if (result) {
      double sqdist = (ret == -1 ? -1. : m_k.squared_distance_d_object()(m_original_points[ret], p));
      *result = std::make_pair(ret, sqdist);
    }

    return ret;
  }*/

  int tree_depth() const
  {
    return m_tree.tree_depth();
  }

private:
  Points_ds                   m_tree;
  K                           m_k;
  Points const&               m_original_points;
};

#endif // FUNCTOR_GUDHI_PERIODIC_KD_TREE_SEARCH_
