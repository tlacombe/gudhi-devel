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

#ifndef FUNCTOR_GUDHI_FLANN_
#define FUNCTOR_GUDHI_FLANN_

#include <flann/flann.hpp>
#include <CGAL/Epick_d.h>

#include <vector>

class Flann
{
  typedef CGAL::Epick_d<CGAL::Dynamic_dimension_tag>  K;
  typedef typename K::Point_d                         Point;
  typedef std::vector<Point>                          Points;

  typedef gss::Kd_tree_search<K, Points>              Points_ds;

public:
  Flann(Points const& points, double /*epsilon*/, flann::IndexParams const& index_params)
    : m_points(create_points_vector(points))
    , m_index(m_points, index_params)
  {
    m_index.buildIndex();
  }

  ~Flann()
  {
    delete[] m_points.ptr();
  }

  // WARNING: eps is ONLY used be KDTreeSingleIndex and KDTreeCuda3dIndex
  void query_k_nearest_neighbors(
    Point const& p,
    unsigned int k,
    double eps = 0.) const
  {
    std::vector<std::vector<int>> indices;
    std::vector<std::vector<double>> dists;
    indices.reserve(1);
    dists.reserve(1);

    m_index.knnSearch(create_point(p), indices, dists, k, flann::SearchParams(32, eps));

#ifdef PRINT_FOUND_NEIGHBORS
    std::cerr << "Query:\n";
    auto dist_it = dists[0].begin();
    for (auto i : indices[0])
    {
      std::cerr << "  " << i << " : " << *dist_it << "\n";
      ++dist_it;
    }
#endif
  }

private:
  flann::Matrix<double> create_point(Point const& p) const
  {
    int dim = K().point_dimension_d_object()(p);
    flann::Matrix<double> pt(new double[dim], 1, dim);

    int c = 0;
    for (auto coord = p.cartesian_begin() ; coord != p.cartesian_end(); ++coord, ++c)
      pt[0][c] = *coord;

    return pt;
  }

  flann::Matrix<double> create_points_vector(std::vector<Point> const& input_points) const
  {
    int dim = K().point_dimension_d_object()(*input_points.begin());
    flann::Matrix<double> pts(new double[input_points.size()*dim], input_points.size(), dim);

    int r = 0;
    for (auto const& p : input_points)
    {
      int c = 0;
      for (auto coord = p.cartesian_begin() ; coord != p.cartesian_end(); ++coord, ++c)
        pts[r][c] = *coord;

      ++r;
    }

    return pts;
  }


  flann::Matrix<double>             m_points;
  flann::Index<flann::L2<double>>   m_index;
};

#endif // FUNCTOR_GUDHI_FLANN_
