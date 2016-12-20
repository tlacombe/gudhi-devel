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
#include <utility>

class Flann
{
  typedef CGAL::Epick_d<CGAL::Dynamic_dimension_tag>  K;
  typedef typename K::Point_d                         Point;
  typedef std::vector<Point>                          Points;

  typedef double                                      Coord_type;
  
public:
  
  // checks: this parameter will be ignored (recomputed) if ground_truth != NULL
  Flann(Points const& points, Coord_type epsilon, flann::IndexParams const& index_params,
    int checks = 32,
    std::vector<std::vector<std::pair<std::size_t, double>>> const *ground_truth = NULL,
    std::vector<Point> const* gt_queries = NULL)
    : m_points(create_points_vector(points))
    , m_index(m_points, index_params)
    , m_checks(checks)
  {
    m_index.buildIndex();
    if (ground_truth)
      auto_tune_precision(epsilon, *gt_queries, *ground_truth);
  }

  ~Flann()
  {
    delete[] m_points.ptr();
  }

  void set_checks(int checks)
  {
    m_checks = checks;
  }

  template <typename Point_range, typename Range_of_query_results>
  void auto_tune_precision(
    Coord_type desired_eps,
    Point_range const& queries,
    Range_of_query_results const& ground_truth)
  {
    std::cerr << "Auto-tuning FLANN precision...\n";
    m_checks = 16; // Will start at 32 since we begin by doubling it (see below)
    double worst_eps = 1.;
    int k = ground_truth.begin()->size();

    auto queries2 = create_points_vector(queries);

    while (worst_eps > desired_eps) {
      flann::SearchParams sp(m_checks, 0.);

      m_checks *= 2;

      std::vector<std::vector<int>> indices;
      std::vector<std::vector<Coord_type>> dists;
      indices.reserve(queries.size());
      dists.reserve(queries.size());

      m_index.knnSearch(queries2, indices, dists, k, sp);

      // Transform the result to the correct format
      std::vector<std::vector<std::pair<std::size_t, double>>> results(
        indices.size(), std::vector<std::pair<std::size_t, double>>());
      auto sqdists_of_one_query_it = dists.begin();
      auto result_it = results.begin();
      for (auto indices_of_one_query : indices) {
        auto sqdist_it = sqdists_of_one_query_it->begin();
        for (auto index : indices_of_one_query) {
          result_it->push_back(std::make_pair(index, *sqdist_it));
          ++sqdist_it;
        }

        ++sqdists_of_one_query_it;
        ++result_it;
      }

      // Check for worst epsilon
      worst_eps = compute_actual_precision(results, ground_truth).first;
    }
    std::cerr << "m_checks = " << m_checks << "\n";
    std::cerr << "Auto-tuning FLANN precision: done!\n";
  }

  // WARNING: eps is ONLY used by KDTreeSingleIndex
  std::size_t query_k_nearest_neighbors(
    Point const& p,
    unsigned int k,
    Coord_type eps = 0.,
    std::vector<std::pair<std::size_t, double>> *result = NULL) const
  {
    std::vector<std::vector<int>> indices;
    std::vector<std::vector<Coord_type>> dists;
    indices.reserve(1);
    dists.reserve(1);

    flann::SearchParams sp(m_checks, eps);
    sp.cores = 1;
    m_index.knnSearch(create_point(p), indices, dists, k, sp);

#ifdef PRINT_FOUND_NEIGHBORS
    std::cerr << "Query:\n";
    auto dist_it = dists[0].begin();
    for (auto i : indices[0])
    {
      std::cerr << "  " << i << " : " << *dist_it << "\n";
      ++dist_it;
    }
#endif

    std::size_t sum = 0;
    if (result) {
      auto dist_it = dists[0].begin();
      for (auto i : indices[0]) {
        sum += i;
        result->push_back(std::make_pair(i, *dist_it));
        ++dist_it;
      }
    }
    else {
      for (auto i : indices[0])
        sum += i;
    }
    return sum;
  }

private:
  flann::Matrix<Coord_type> create_point(Point const& p) const
  {
    int dim = K().point_dimension_d_object()(p);
    flann::Matrix<Coord_type> pt(new Coord_type[dim], 1, dim);

    int c = 0;
    for (auto coord = p.cartesian_begin() ; coord != p.cartesian_end(); ++coord, ++c)
      pt[0][c] = *coord;

    return pt;
  }

  template <typename Point_range>
  flann::Matrix<Coord_type> create_points_vector(Point_range const& input_points) const
  {
    int dim = K().point_dimension_d_object()(*input_points.begin());
    flann::Matrix<Coord_type> pts(new Coord_type[input_points.size()*dim], input_points.size(), dim);

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


  flann::Matrix<Coord_type>             m_points;
  flann::Index<flann::L2<Coord_type>>   m_index;
  int m_checks;
};

#endif // FUNCTOR_GUDHI_FLANN_
