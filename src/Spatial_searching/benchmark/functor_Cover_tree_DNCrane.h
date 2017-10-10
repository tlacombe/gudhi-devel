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

#ifndef FUNCTOR_COVER_TREE_DNCRANE_
#define FUNCTOR_COVER_TREE_DNCRANE_

#include <Cover_Tree.h>
#include <CGAL/Epick_d.h>

#include <vector>
#include <utility>

class Cover_tree_DNCrane
{
  typedef CGAL::Epick_d<CGAL::Dynamic_dimension_tag>  K;
  typedef typename K::Point_d                         Point;
  typedef std::vector<Point>                          Points;

  typedef double                                      Coord_type;
  
public:

  struct MyPoint : public Point
  {
    MyPoint(Point const& p)
      : Point(p)
    {}

    double distance(const MyPoint& p) const
    {
      return std::sqrt(K().squared_distance_d_object()(*this, p));
    }

    bool operator==(const MyPoint& p) const
    {
      return static_cast<Point>(*this) == static_cast<Point>(p);
    }
  };

  Cover_tree_DNCrane(Points const& points, Coord_type)
    : m_tree(10000., create_points_vector(points)) // CJOTOD max distance = ????
  {
  }

  // WARNING: eps is IGNORED
  std::size_t query_k_nearest_neighbors(
    Point const& p,
    unsigned int k,
    Coord_type eps = 0., // ignored
    std::vector<std::pair<std::size_t, double>> *result = NULL) const
  {
    std::vector<MyPoint> points_found = m_tree.kNearestNeighbors(create_point(p), k);

#ifdef PRINT_FOUND_NEIGHBORS
    std::cerr << "Query: " << p << "\n";
    for (auto pt : points_found)
      std::cerr << "  " << pt << " - sq. dist = " << create_point(p).distance(pt) * create_point(p).distance(pt) << "\n";
#endif
    return 0;
  }

  int tree_depth() const
  {
    return -1;
  }

private:
  MyPoint create_point(Point const& p) const
  {
    return MyPoint(p);
  }

  template <typename Point_range>
  std::vector<MyPoint> create_points_vector(Point_range const& input_points) const
  {
    std::vector<MyPoint> points;
    points.reserve(input_points.size());

    for (auto const& p : input_points)
      points.emplace_back(create_point(p));

    return points;
  }


  CoverTree<MyPoint> m_tree;
};

#endif // FUNCTOR_COVER_TREE_DNCRANE_
