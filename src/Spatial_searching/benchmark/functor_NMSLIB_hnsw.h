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

#ifdef GUDHI_NMSLIB_IS_AVAILABLE

#ifndef FUNCTOR_NMSLIB_HNSW_
#define FUNCTOR_NMSLIB_HNSW_

// NMSLIB
#include <params.h>
#include <object.h>
#include <knnquery.h>
#include <space/space_vector.h>
#include <factory/space/space_lp.h>
#include <method/hnsw.h>

#include <memory>
#include <vector>

namespace nms = similarity;

class NMSLIB_hnsw
{
  typedef CGAL::Epick_d<CGAL::Dynamic_dimension_tag>      K;
  typedef typename K::Point_d                             Point;
  typedef std::vector<Point>                              Points;

  typedef float                                          Dist_type; // CJTODO
  typedef std::unique_ptr<nms::SpaceLp<Dist_type>>        Space;
  typedef nms::Hnsw<Dist_type>                            Hnsw;

public:
  NMSLIB_hnsw(Points const& points, double /*epsilon*/)
  : m_space(std::unique_ptr<nms::SpaceLp<Dist_type>>(static_cast<nms::SpaceLp<Dist_type>*>(nms::CreateL2<Dist_type>(nms::AnyParams())))),
    m_points(create_points_vector(points)),
    m_hnsw(/*PrintProgress =*/ false, *m_space, m_points)
  {
    nms::AnyParams index_params({
      "M=10",
      "efConstruction=20",
      "indexThreadQty=1",
      "searchMethod=0" });
    m_hnsw.CreateIndex(index_params);

    nms::AnyParams query_time_params({
      //"algoType=2",
      "efSearch=100" }); /// ???
    m_hnsw.SetQueryTimeParams(query_time_params);
  }

  void query_k_nearest_neighbors(
    Point const& p,
    unsigned int k,
    double eps = 0.) const
  {
    int dummy = 0;
    nms::Object *pt = create_point(p);
    nms::KNNQuery<Dist_type> q(*m_space, pt, k, eps);
    m_hnsw.Search(&q, dummy);
#ifdef PRINT_FOUND_NEIGHBORS
    q.Print();
#endif
    delete pt;
  }

private:

  nms::Object *create_point(Point const& p, int index = 0) const
  {
    int dummy = 0;
    // CJTODO: this temporary object creation is costly
    return m_space->CreateObjFromVect(
      index, dummy, std::vector<Dist_type>(p.cartesian_begin(), p.cartesian_end()));
  }

  nms::ObjectVector create_points_vector(std::vector<Point> const& input_points) const
  {
    nms::ObjectVector points;
    points.reserve(input_points.size());

    int c = 0;
    for (auto const& p : input_points)
    {
      points.push_back(create_point(p, c));
      ++c;
    }

    return points;
  }

  Space m_space;
  Hnsw m_hnsw;
  nms::ObjectVector m_points;
};

#endif // FUNCTOR_NMSLIB_HNSW_
#endif // GUDHI_NMSLIB_IS_AVAILABLE