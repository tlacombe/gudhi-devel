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

#ifndef FUNCTOR_SBL_PROXIMITY_FOREST_
#define FUNCTOR_SBL_PROXIMITY_FOREST_

#include <SBL/GT/ANN_metric_space_proximity_tree.hpp>
#include <SBL/GT/ANN_meta.hpp>

#include <CGAL/Epick_d.h>
#include <CGAL/iterator.h>

#include <boost/iterator/counting_iterator.hpp>
#include <boost/iterator/zip_iterator.hpp>

#include <vector>
#include <utility>

namespace sbl = SBL::GT;

class SBL_Proximity_Forest
{
  typedef CGAL::Epick_d<CGAL::Dynamic_dimension_tag>               K;
  typedef K::Point_d                                               Point;
  typedef K::FT                                                    FT;
  typedef std::vector<Point>                                       Points;
  typedef std::pair<int, Point const*>                             SBL_point;

  class Squared_distance
  {
    K::Squared_distance_d m_sq_dist;

  public:
    typedef K::FT FT;
    typedef SBL_point Point;

    FT operator()(Point p, Point q)const
    {
      return m_sq_dist(p.second[p.first], q.second[q.first]);
    }
  };

  typedef sbl::T_ANN_metric_space_proximity_tree<Squared_distance> Proximity_tree;
  typedef sbl::T_ANN_meta<Proximity_tree>                          Proximity_forest;

public:
  SBL_Proximity_Forest(Points const& points, double /*epsilon*/ = 0., unsigned int num_trees = 10)
    : m_points(points)
    , m_forest(m_sq_dist_functor, num_trees)
  {
    std::vector<SBL_point> sbl_points;
    sbl_points.reserve(m_points.size());
    for (int i = 0; i < m_points.size(); ++i)
      sbl_points.push_back(std::make_pair(i, m_points.data()));
    m_forest.insert(sbl_points.begin(), sbl_points.end());
    /*m_forest.insert(
      boost::make_zip_iterator(
        boost::make_tuple(boost::counting_iterator<std::ptrdiff_t>(0), CGAL::Const_oneset_iterator<Point const*>(m_points.data()))),
      boost::make_zip_iterator(
        boost::make_tuple(boost::counting_iterator<std::ptrdiff_t>(points.size()), CGAL::Const_oneset_iterator<Point const*>(m_points.data()))));*/
    m_forest.set_query_type(Proximity_forest::ANN_BY_K);
  }

  // Returns the sum of the indices
  std::size_t query_k_nearest_neighbors(
    Point const& p,
    unsigned int k,
    double eps = 0.,
    std::vector<std::pair<std::size_t, double>> *result = NULL)
  {
    m_forest.set_number_of_neighbors(k);
    m_forest.set_range(3);

    std::vector<SBL_point> neighbors;
    neighbors.reserve(k);
    SBL_point sbl_p = std::make_pair(0, &p);
    m_forest(sbl_p, true, std::back_inserter(neighbors));

    std::size_t sum = 0;
    if (result) {
      for (auto n : neighbors) {
        sum += n.first;
        result->push_back(std::make_pair(n.first, m_k.squared_distance_d_object()(m_points[n.first], p)));
      }
    }
    else {
      for (auto n : neighbors)
        sum += n.first;
    }

#ifdef PRINT_FOUND_NEIGHBORS
    std::cerr << "Query:\n";
    for (auto n : neighbors) {
      std::cerr << "  " << n.first << " : " << m_k.squared_distance_d_object()(m_points[n.first], p) << "\n";
    }
#endif

    return sum;
  }

  int tree_depth() const
  {
    return -1;
  }

private:
  Points const& m_points;
  Squared_distance m_sq_dist_functor;
  Proximity_forest m_forest;
  K m_k;
};

#endif // FUNCTOR_SBL_PROXIMITY_FOREST_
