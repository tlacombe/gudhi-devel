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

#ifndef FUNCTOR_NANOFLANN_
#define FUNCTOR_NANOFLANN_

#include <nanoflann.hpp>
#include <CGAL/Epick_d.h>

#include <vector>
#include <utility>

 // "dataset to kd-tree" adaptor class
template <typename K, typename Point_container_>
class Point_cloud_adaptator__nanoflann
{
public:
  typedef typename Point_container_::value_type         Point;
  typedef K                                             Kernel;

  /// The constructor that sets the data set source
  Point_cloud_adaptator__nanoflann(Point_container_ const& points, Kernel const& k)
    : m_points(points), m_k(k)
  {}

  /// CRTP helper method
  inline Point_container_ const& points() const
  {
    return m_points;
  }
  inline Point_container_& points()
  {
    return m_points;
  }

  // Must return the number of data points
  inline size_t kdtree_get_point_count() const
  {
    return m_points.size();
  }

  // Returns the distance between the vector "p1[0:size-1]"
  // and the data point with index "idx_p2" stored in the class:
  inline double kdtree_distance(
    const double *p1, const size_t idx_p2, size_t size) const
  {
    Point sp(p1, p1 + size);
    return m_k.squared_distance_d_object()(sp, points()[idx_p2]);
  }

  // Returns the dim'th component of the idx'th point in the class:
  // Since this is inlined and the "dim" argument is typically an
  // immediate value, the "if/else's" are actually solved at compile time.
  inline double kdtree_get_pt(const size_t idx, int dim) const
  {
    return m_k.compute_coordinate_d_object()(points()[idx], dim);
  }

  // Optional bounding-box computation: return false to default to a standard
  // bbox computation loop.
  // Return true if the BBOX was already computed by the class and returned
  // in "bb" so it can be avoided to redo it again.
  // Look at bb.size() to find out the expected dimensionality
  // (e.g. 2 or 3 for point clouds)
  template <class Bbox>
  bool kdtree_get_bbox(Bbox &bb) const
  {
    return false;
  }

  Kernel const& kernel() const
  {
    return m_k;
  }

protected:
  Point_container_ const& m_points; //!< A const ref to the data set origin
  Kernel const& m_k;      //!< A const ref to the kernel

};

//////////////////////////////////////////////////////////////////////

class Nanoflann
{
  typedef CGAL::Epick_d<CGAL::Dynamic_dimension_tag>  K;
  typedef typename K::Point_d                         Point;
  typedef std::vector<Point>                          Points;

public:
  /// Constructor
  /// "points" must not be empty
  Nanoflann(Points const& points, double /*epsilon*/)
    : m_adaptor(points, m_k),
      m_kd_tree(
        m_k.point_dimension_d_object()(*points.begin()),
        m_adaptor,
        nanoflann::KDTreeSingleIndexAdaptorParams(10 /* max leaf */)
      )
  {
    m_kd_tree.buildIndex();
  }

  std::size_t query_k_nearest_neighbors(
    Point const& p,
    unsigned int k,
    double eps = 0.,
    std::vector<std::pair<std::size_t, double>> *result = NULL) const
  {
    size_t *neighbor_indices = new size_t[k];
    double *squared_distance = new double[k];

    std::vector<double> sp_vec(
      m_adaptor.kernel().construct_cartesian_const_iterator_d_object()(p),
      m_adaptor.kernel().construct_cartesian_const_iterator_d_object()(p, 0));
    nanoflann::KNNResultSet<double> result_set(k);
    result_set.init(neighbor_indices, squared_distance);
    m_kd_tree.findNeighbors(result_set,
      &sp_vec[0],
      nanoflann::SearchParams(32, eps));

#ifdef PRINT_FOUND_NEIGHBORS
    std::cerr << "Query:\n";
    for (int i = 0 ; i < k ; ++i)
      std::cerr << "  " << neighbor_indices[i] << " : " << squared_distance[i] << "\n";
#endif

    std::size_t sum = 0;
    if (result) {
      for (int i = 0; i < k; ++i) {
        sum += neighbor_indices[i];
        result->push_back(std::make_pair(neighbor_indices[i], squared_distance[i]));
      }
    }
    else {
      for (int i = 0; i < k; ++i)
        sum += neighbor_indices[i];
    }
    return sum;
  }

  /*void query_ball(const Point &sp,
    const double radius,
    std::vector<std::pair<std::size_t, double> > &neighbors,
    bool sort_output = true)
  {
    std::vector<double> sp_vec(
      m_adaptor.kernel().construct_cartesian_const_iterator_d_object()(sp),
      m_adaptor.kernel().construct_cartesian_const_iterator_d_object()(sp, 0));
    m_kd_tree.radiusSearch(&sp_vec[0],
      radius,
      neighbors,
      nanoflann::SearchParams(32, 0.f, sort_output));

    std::cout << "radiusSearch(num="<< neighbors.size() <<"): \n";
    for (const auto idx_and_dist : neighbors)
    {
    std::cout << "  * neighbor_indices = " << idx_and_dist.first
    << " (out_dist_sqr = " << idx_and_dist.second << ")"
    << std::endl;
    }
  }*/

private:
  typedef Point_cloud_adaptator__nanoflann<K, Points> Adaptor;
  typedef nanoflann::KDTreeSingleIndexAdaptor<
    nanoflann::L2_Simple_Adaptor<double, Adaptor>,
    Adaptor,
    -1 // dim
  > Kd_tree;

  K m_k;
  Adaptor m_adaptor;
  Kd_tree m_kd_tree;
};

#endif // FUNCTOR_NANOFLANN_
