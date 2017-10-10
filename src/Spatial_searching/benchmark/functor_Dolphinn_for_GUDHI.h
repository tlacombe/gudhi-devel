/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Clement Jamin
 *
 *    Copyright (C) 2016 INRIA
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

#ifndef FUNCTOR_DOLPHINN_FOR_GUDHI_
#define FUNCTOR_DOLPHINN_FOR_GUDHI_

#include <gudhi/Dolphinn.h>
#include <CGAL/Epick_d.h>

#include <cmath>
#include <algorithm> // for std::max
#include <vector>
#include <utility>

namespace dolphinn = Gudhi::dolphinn;

template <typename Kernel>
class Dolphinn_for_GUDHI
{
  typedef Kernel                                      K;
  typedef typename K::Point_d                         Point;
  typedef std::vector<Point>                          Points;

  typedef double                                      Coord_type;
  typedef std::vector<Coord_type>                     Dolphinn_point;

public:

  // hashing_method: if positive, the parameter of Stable Distribution, if 0, the LSH method used is the hyperplanes. 
  Dolphinn_for_GUDHI(Points const& points, Coord_type /*epsilon*/, double hashing_method = 0., std::size_t max_pts_to_search = 0)
    : m_dim(K().point_dimension_d_object()(*points.begin()))
    , m_k(m_dim)
    , m_points(create_points_vector(points))
    , m_dolphi(m_points, points.size(), m_dim, floor(log2(points.size())/2), hashing_method)
    , m_original_points(points)
    , m_max_pts_to_search(max_pts_to_search)
  {}

  ~Dolphinn_for_GUDHI()
  {
  }

  std::size_t query_k_nearest_neighbors(
    Point const& p,
    unsigned int k,
    Coord_type eps = 0.,
    std::vector<std::pair<std::size_t, double>> *result = NULL) // CJTODO: const
  {
    std::vector<std::vector<std::pair<int, double>>> neighbors;
    std::vector<std::pair<int, double>> dummy;
    neighbors.push_back(dummy);

    std::vector<Dolphinn_point> queries;
    queries.push_back(create_point(p));

    m_dolphi.k_nearest_neighbors_query(
      queries, queries.size(), k, 
      m_max_pts_to_search == 0 ? m_points.size()/100 + k : m_max_pts_to_search, 
      neighbors);

#ifdef PRINT_FOUND_NEIGHBORS
    std::cerr << "Query:\n";
    for (auto nb : *neighbors.begin())
      std::cerr << "  " << nb.first << " : " << nb.second << "\n";
#endif

    std::size_t sum = 0;
    if (result) {
      for (auto nb : *neighbors.begin())
      {
        sum += nb.first;
        result->push_back(nb);
      }
    }
    else {
      for (auto nb : *neighbors.begin())
        sum += nb.first;
    }
    return sum;
  }

  // result = <point_index, squared_distance>
  std::ptrdiff_t query_any_near_neighbor(
    Point const& p,
    double radius,
    Coord_type /*eps*/ = 0.,
    std::pair<std::ptrdiff_t, double> *result = NULL) /*const CJTODO*/
  {
    std::vector<int> neighbors(1);
    *neighbors.begin() = -1;

    std::vector<Dolphinn_point> queries;
    queries.push_back(create_point(p));

    m_dolphi.radius_query(
      queries, queries.size(), radius, 
      m_max_pts_to_search == 0 ? (std::max)(static_cast<std::size_t>(10), m_points.size() / 100 : m_max_pts_to_search),
      neighbors);

    int result_point_index = *neighbors.begin();

#ifdef PRINT_FOUND_NEIGHBORS
    if (result_point_index == -1)
      std::cerr << "Any near neighbor: none\n";
    else
      std::cerr << "Any near neighbor: " << result_point_index << ": "
        << m_k.squared_distance_d_object()(m_original_points[result_point_index], p) << "\n";
#endif

    if (result) {
      double sqdist = (result_point_index == -1 ? -1. : m_k.squared_distance_d_object()(m_original_points[result_point_index], p));
      *result = std::make_pair(result_point_index, sqdist);
    }

    return result_point_index;
  }

  std::size_t query_all_near_neighbors(
    Point const& p,
    double radius,
    double /*eps*/ = 0.,
    std::vector<std::pair<std::size_t, double>> *result = NULL) // CJTODO: const
  {
    std::vector<std::vector<int>> neighbors(1);

    std::vector<Dolphinn_point> queries;
    queries.push_back(create_point(p));

    m_dolphi.all_radius_query(
      queries, queries.size(), radius,
      m_max_pts_to_search == 0 ? (std::max)(static_cast<std::size_t>(10), m_points.size() / 100) : m_max_pts_to_search,
      neighbors);

#ifdef PRINT_FOUND_NEIGHBORS
    std::cerr << "Query:\n";
    for (auto nb : neighbors[0]) {
      double sqdist = m_k.squared_distance_d_object()(m_original_points[nb], p);
      std::cerr << "  " << nb << " : " << sqdist << "\n";
    }
#endif

    std::size_t sum = 0;
    if (result) {
      for (auto nb : neighbors[0])
      {
        sum += nb;
        double sqdist = m_k.squared_distance_d_object()(m_original_points[nb], p);
        result->push_back(std::make_pair(nb, sqdist));
      }
    }
    else {
      for (auto nb : neighbors[0])
        sum += nb;
    }

    return sum;
  }

  int tree_depth() const
  {
    return -1; // Not provided by the library.
  }

private:
  Dolphinn_point create_point(Point const& p) const
  {
    Dolphinn_point pt;
    pt.reserve(m_dim);

    int i = 0;
    for (auto coord = p.cartesian_begin(); coord != p.cartesian_end(); ++coord, ++i)
      pt.push_back(*coord);

    return pt;
  }

  template <typename Point_range>
  std::vector<Dolphinn_point> create_points_vector(Point_range const& input_points) const
  {
    std::vector<Dolphinn_point> points;
    points.reserve(input_points.size());

    for (auto const& p : input_points)
    {
      Dolphinn_point pt;
      pt.reserve(m_dim);
      int j = 0;
      for (auto coord = p.cartesian_begin(); coord != p.cartesian_end(); ++coord, ++j)
        pt.push_back(*coord);

      points.push_back(pt);
    }

    return points;
  }

  int m_dim;
  K m_k;
  std::vector<Dolphinn_point> m_points;
  dolphinn::Dolphinn<Coord_type> m_dolphi;
  Points const& m_original_points;
  std::size_t m_max_pts_to_search;
};

#endif // FUNCTOR_DOLPHINN_FOR_GUDHI_
