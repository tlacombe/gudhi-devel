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

#ifndef FUNCTOR_GUDHI_KD_TREE_SEARCH_
#define FUNCTOR_GUDHI_KD_TREE_SEARCH_

#include <gudhi/Kd_tree_search.h>
#include <CGAL/Epick_d.h>

#include <vector>
#include <utility>

namespace gss = Gudhi::spatial_searching;

class GUDHI_Kd_tree_search
{
  typedef CGAL::Epick_d<CGAL::Dynamic_dimension_tag>  K;
  typedef K::Point_d                                  Point;
  typedef std::vector<Point>                          Points;

  typedef gss::Kd_tree_search<K, Points>              Points_ds;

public:
  GUDHI_Kd_tree_search(Points const& points, double /*epsilon*/ = 0.)
    : m_tree(points)
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
    std::cerr << "Query:\n";
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

  int tree_depth() const
  {
    return m_tree.tree_depth();
  }

private:
  Points_ds m_tree;
};

#endif // FUNCTOR_GUDHI_KD_TREE_SEARCH_
