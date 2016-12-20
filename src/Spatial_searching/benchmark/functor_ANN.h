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

#ifndef FUNCTOR_GUDHI_ANN_
#define FUNCTOR_GUDHI_ANN_

#include <ANN/ANN.h>
#include <CGAL/Epick_d.h>

#include <vector>
#include <utility>

class Ann
{
  typedef CGAL::Epick_d<CGAL::Dynamic_dimension_tag>  K;
  typedef typename K::Point_d                         Point;
  typedef std::vector<Point>                          Points;

  typedef double                                      Coord_type;
  
public:
  
  // checks: this parameter will be ignored (recomputed) if ground_truth != NULL
  Ann(Points const& points, Coord_type epsilon,
    std::vector<std::vector<std::pair<std::size_t, double>>> const *ground_truth = NULL,
    std::vector<Point> const* gt_queries = NULL)
    : m_points(create_points_vector(points))
    , m_dim(K().point_dimension_d_object()(*points.begin()))
    , m_tree(m_points, points.size(), m_dim)
  {}

  std::size_t query_k_nearest_neighbors(
    Point const& p,
    unsigned int k,
    Coord_type eps = 0.,
    std::vector<std::pair<std::size_t, double>> *result = NULL)
  {
    std::vector<int> neighbors_indices(k);
    std::vector<double> neighbors_sq_distances(k);

    m_tree.annkSearch(                // search
      create_point(p),                // query point
      k,                              // number of near neighbors
      neighbors_indices.data(),       // nearest neighbors (returned)
      neighbors_sq_distances.data(),  // distance (returned)
      eps);                             // error bound


#ifdef PRINT_FOUND_NEIGHBORS
    std::cerr << "Query:\n";
    auto dist_it = neighbors_sq_distances.begin();
    for (auto i : neighbors_indices)
    {
      std::cerr << "  " << i << " : " << *dist_it << "\n";
      ++dist_it;
    }
#endif

    std::size_t sum = 0;
    if (result) {
      auto dist_it = neighbors_sq_distances.begin();
      for (auto i : neighbors_indices) {
        sum += i;
        result->push_back(std::make_pair(i, *dist_it));
        ++dist_it;
      }
    }
    else {
      for (auto i : neighbors_indices)
        sum += i;
    }
    return sum;
  }

private:
  ANNpoint create_point(Point const& p) const
  {
    ANNpoint pt = annAllocPt(m_dim);

    int i = 0;
    for (auto coord = p.cartesian_begin() ; coord != p.cartesian_end(); ++coord, ++i)
      pt[i] = *coord;

    return pt;
  }

  template <typename Point_range>
  ANNpointArray create_points_vector(Point_range const& input_points) const
  {
    ANNpointArray points = annAllocPts(input_points.size(), m_dim);

    int i = 0;
    for (auto const& p : input_points)
    {
      int j = 0;
      for (auto coord = p.cartesian_begin(); coord != p.cartesian_end(); ++coord, ++j)
        points[i][j] = *coord;

      ++i;
    }

    return points;
  }

  int m_dim;
  ANNpointArray m_points;
  ANNkd_tree m_tree;
};

#endif // FUNCTOR_GUDHI_ANN_
