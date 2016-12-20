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

#ifndef FUNCTOR_GUDHI_FALCONN_
#define FUNCTOR_GUDHI_FALCONN_

#include <falconn/lsh_nn_table.h>
#include <CGAL/Epick_d.h>

#include <vector>
#include <utility>

class Falconn
{
  typedef CGAL::Epick_d<CGAL::Dynamic_dimension_tag>  K;
  typedef typename K::Point_d                         Point;
  typedef std::vector<Point>                          Points;

  typedef double                                      Coord_type;
  typedef falconn::DenseVector<Coord_type>            Falconn_point;
  
public:
  
  // checks: this parameter will be ignored (recomputed) if ground_truth != NULL
  Falconn(Points const& points, Coord_type /*epsilon*/, int num_probes,
    std::vector<std::vector<std::pair<std::size_t, double>>> const *ground_truth = NULL,
    std::vector<Point> const* gt_queries = NULL)
    : m_dim(K().point_dimension_d_object()(*points.begin()))
    , m_k(m_dim)
    , m_points(create_points_vector(points))
    , m_table(falconn::construct_table<Falconn_point>(m_points, get_params(m_points.size(), m_dim)))
    , m_original_points(points)
  {
    m_table->set_num_probes(num_probes); // CJTODO
  }

  std::size_t query_k_nearest_neighbors(
    Point const& p,
    unsigned int k,
    Coord_type eps = 0.,
    std::vector<std::pair<std::size_t, double>> *result = NULL)
  {
    std::vector<int> neighbors_indices;
    neighbors_indices.reserve(k);
    std::vector<double> neighbors_sq_distances;
    neighbors_sq_distances.reserve(k);

    m_table->find_k_nearest_neighbors(create_point(p), k, &neighbors_indices);

    std::size_t sum = 0;
    if (result) {
      for (auto i : neighbors_indices) {
        sum += i;
        result->push_back(std::make_pair(i, m_k.squared_distance_d_object()(m_original_points[i], p)));
      }
    }
    else {
      for (auto i : neighbors_indices)
        sum += i;
    }

#ifdef PRINT_FOUND_NEIGHBORS
    std::cerr << "Query:\n";
    for (auto i : neighbors_indices) {
      std::cerr << "  " << i << " : " << m_k.squared_distance_d_object()(m_original_points[i], p) << "\n";
    }
#endif

    return sum;
  }

private:
  falconn::LSHConstructionParameters get_params(std::size_t num_points, int dim)
  {
    falconn::LSHConstructionParameters params =
      falconn::get_default_parameters<Falconn_point>(num_points, dim, falconn::DistanceFunction::EuclideanSquared, true);
    params.num_setup_threads = 1;
    return params;
  }

  Falconn_point create_point(Point const& p) const
  {
    Falconn_point pt(m_dim);

    int i = 0;
    for (auto coord = p.cartesian_begin() ; coord != p.cartesian_end(); ++coord, ++i)
      pt[i] = *coord;

    return pt;
  }

  template <typename Point_range>
  std::vector<Falconn_point> create_points_vector(Point_range const& input_points) const
  {
    std::vector<Falconn_point> points;
    points.reserve(input_points.size());

    for (auto const& p : input_points)
      points.push_back(create_point(p));

    return points;
  }

  int m_dim;
  K m_k;
  std::vector<Falconn_point> m_points;
  std::unique_ptr<falconn::LSHNearestNeighborTable<Falconn_point, int>> m_table;
  Points const& m_original_points;
};

#endif // FUNCTOR_GUDHI_FALCONN_
