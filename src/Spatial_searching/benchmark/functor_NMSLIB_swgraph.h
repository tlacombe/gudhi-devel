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

#ifndef FUNCTOR_NMSLIB_SWGRAPH_
#define FUNCTOR_NMSLIB_SWGRAPH_

// NMSLIB
#include <params.h>
#include <object.h>
#include <knnquery.h>
#include <knnqueue.h>
#include <space/space_vector.h>
#include <factory/space/space_lp.h>
#include <method/small_world_rand.h>

#include <CGAL/Epick_d.h>

#include <memory>
#include <vector>
#include <utility>

namespace nms = similarity;

class NMSLIB_swgraph
{
  typedef CGAL::Epick_d<CGAL::Dynamic_dimension_tag>      K;
  typedef typename K::Point_d                             Point;
  typedef std::vector<Point>                              Points;

  typedef double                                          Dist_type;
  typedef std::unique_ptr<nms::SpaceLp<Dist_type>>        Space;
  typedef nms::SmallWorldRand<Dist_type>                  SWgraph;

public:
  NMSLIB_swgraph(Points const& points, double /*epsilon*/)
  : m_space(std::unique_ptr<nms::SpaceLp<Dist_type>>(static_cast<nms::SpaceLp<Dist_type>*>(nms::CreateL2<Dist_type>(nms::AnyParams())))),
    m_points(create_points_vector(points)),
    m_swgraph(/*PrintProgress =*/ false, *m_space, m_points)
  {
    nms::AnyParams index_params({
      "NN=10",
      "efConstruction=20",
      "indexThreadQty=1" });
    m_swgraph.CreateIndex(index_params);

    nms::AnyParams query_time_params({
      "algoType=v1merge", // or "old"
      "efSearch=100" });
    m_swgraph.SetQueryTimeParams(query_time_params);
  }

  std::size_t query_k_nearest_neighbors(
    Point const& p,
    unsigned int k,
    double eps = 0.,
    std::vector<std::pair<std::size_t, double>> *result = NULL) const
  {
    int dummy = 0;
    nms::Object *pt = create_point(p);
    nms::KNNQuery<Dist_type> q(*m_space, pt, k, eps);
    m_swgraph.Search(&q, dummy);
    delete pt;

#ifdef PRINT_FOUND_NEIGHBORS
    std::cerr << "Query:\n";
    nms::KNNQueue<Dist_type> *res2 = q.Result()->Clone();
    while (!res2->Empty()) {
      std::cerr << "  " << res2->TopObject()->id() 
        << " : " << res2->TopDistance()*res2->TopDistance() << "\n";
      res2->Pop();
    }
#endif

    std::size_t sum = 0;
    nms::KNNQueue<Dist_type> *res = q.Result()->Clone();
    if (result) {
      result->resize(k);
      int i = k - 1;
      while (!res->Empty()) {
        sum += res->TopObject()->id();
        (*result)[i] = std::make_pair(res->TopObject()->id(), res->TopDistance()*res->TopDistance());
        res->Pop();
        --i;
      }
    }
    else {
      while (!res->Empty()) {
        sum += res->TopObject()->id();
        res->Pop();
      }
    }
    return sum;
  }

private:

  nms::Object *create_point(Point const& p, int index = 0) const
  {
    int dummy = 0;
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
  nms::ObjectVector m_points;
  SWgraph m_swgraph;
};

#endif // FUNCTOR_NMSLIB_SWGRAPH_
#endif // GUDHI_NMSLIB_IS_AVAILABLE