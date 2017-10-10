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

#ifndef FUNCTOR_COVER_TREE_MANZIL_
#define FUNCTOR_COVER_TREE_MANZIL_

#include <cover_tree_manzil/cover_tree.h>
#include <CGAL/Epick_d.h>

#include <vector>
#include <utility>

class Cover_tree_Manzil
{
  typedef CGAL::Epick_d<CGAL::Dynamic_dimension_tag>  K;
  typedef typename K::Point_d                         Point;
  typedef std::vector<Point>                          Points;

  typedef double                                      Coord_type;
  
public:

  Cover_tree_Manzil(Points const& points, Coord_type)
    : m_tree(create_points_vector(points))
  {
  }

  // WARNING: eps is IGNORED
  std::size_t query_k_nearest_neighbors(
    Point const& p,
    unsigned int k,
    Coord_type eps = 0., // ignored
    std::vector<std::pair<std::size_t, double>> *result = NULL) const
  {
    std::vector<Eigen::VectorXd> points_found = m_tree.nearNeighbors(create_point(p), k);

#ifdef PRINT_FOUND_NEIGHBORS
    std::cerr << "Query: " << p << "\n";
    for (auto p : points_found)
      std::cerr << "  " << p << "\n";
#endif
    return 0;
  }

  int tree_depth() const
  {
    return -1;
  }

private:
  Eigen::VectorXd create_point(Point const& p) const
  {
    int dim = K().point_dimension_d_object()(p);
    Eigen::VectorXd pt(dim);

    int c = 0;
    for (auto coord = p.cartesian_begin() ; coord != p.cartesian_end(); ++coord, ++c)
      pt(c) = *coord;

    return pt;
  }

  template <typename Point_range>
  std::vector<Eigen::VectorXd> create_points_vector(Point_range const& input_points) const
  {
    std::vector<Eigen::VectorXd> points;
    points.reserve(input_points.size());

    for (auto const& p : input_points)
      points.emplace_back(create_point(p));

    return points;
  }


  CoverTree m_tree;
};

#endif // FUNCTOR_COVER_TREE_MANZIL_
