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
#include <knnqueue.h>
#include <space/space_vector.h>
#include <factory/space/space_lp.h>
#include <method/hnsw.h>

#include <CGAL/Epick_d.h>

#include <memory>
#include <vector>
#include <utility>

namespace nms = similarity;

class NMSLIB_hnsw
{
  typedef CGAL::Epick_d<CGAL::Dynamic_dimension_tag>      K;
  typedef typename K::Point_d                             Point;
  typedef std::vector<Point>                              Points;

  typedef double                                          Dist_type;
  typedef std::unique_ptr<nms::SpaceLp<Dist_type>>        Space;
  typedef nms::Hnsw<Dist_type>                            Hnsw;

public:
  // M, post, etc.: see below
  NMSLIB_hnsw(Points const& points, double /*epsilon*/, int M = 16, int post = 0, int efConstruction = 200, int efSearch = 100)
  : m_space(std::unique_ptr<nms::SpaceLp<Dist_type>>(static_cast<nms::SpaceLp<Dist_type>*>(nms::CreateL2<Dist_type>(nms::AnyParams())))),
    m_points(create_points_vector(points)),
    m_hnsw(/*PrintProgress =*/ false, *m_space, m_points)
  {
    std::cerr << "M=" << M << ", post=" << post << ", efConstruction=" << efConstruction << ", efSearch=" << efSearch << "\n";

    // ann-benchmark uses: (M, post, efSearch)
    // (32, 2, [20, 30, 40, 50, 60, 70, 80, 90, 100, 120, 140, 160, 200, 300, 400]),
    // (20, 2, [2, 5, 10, 15, 20, 30, 40, 50, 70, 80, 120, 200, 400]),
    //  (12, 0, [1, 2, 5, 10, 15, 20, 30, 40, 50, 70, 80, 120]),
    //  (4, 0, [1, 2, 5, 10, 20, 30, 50, 70, 90, 120]),
    //  (8, 0, [1, 2, 5, 10, 20, 30, 50, 70, 90, 120, 160])
    nms::AnyParams index_params({
      // 5-100 (higer = better).
      // Performance impact: 
      // * built time = small impact
      // * query time = medium impact
      // The size of the initial set of potential neighbors for the indexing
      // phase. The set may be further pruned so that the overall number
      // of neighbors does not exceed maxM0(for the ground layer)
      // or maxM (for all layers but the ground one).
      std::string("M=") + std::to_string(M),
      // 0 or 2.
      // The post=2 adds an additional post-processing step to symmetrize 
      // the index. It is not documented yet in nmslib and do not presented in
      // the paper on HNSW. In short, two indexes is built with different order 
      // of the data(direct and reverse), with following union of the produced 
      // connections. It does not affect directly the query time parameters, such as ef.
      // On overall, it leads to a roughly twice as long construction time, while 
      // in the end adding extra search performance at high recalls(up to several
      // tens of percent in tests on 1M sift at ~0.999 recall, the gain smaller for 
      // less recall).So if you want to get maximum performance at search you 
      // should use post = 2.
      std::string("post=") + std::to_string(post),
      // 100-2000 (higher = better).
      // Performance impact: 
      // * built time = big impact (~linear)
      // * query time = no impact
      // The depth of the search that is used to find neighbors during indexing
      // (this parameter is used only for the search in the ground layer).
      // *** ann-benchmark uses: 400 ***
      std::string("efConstruction=") + std::to_string(efConstruction),
      // Number of threads.
      "indexThreadQty=1",
      // Undocumented. Default: 0.
      "searchMethod=0" });
    m_hnsw.CreateIndex(index_params);

    nms::AnyParams query_time_params({
      // "old" or "v1merge". Undocumented.
      //"algoType=v1merge",

      // 100-2000 (higher = better).
      // Performance impact: 
      // * built time = no impact
      // * query time = big impact (almost linear, depending on other params)
      // Same as "ef". The search depth: specifically, a sub-search is stopped, when it
      // cannot find a point closer than efSearch points (seen so far) closest to the query.
      std::string("efSearch=") + std::to_string(efSearch)});
    m_hnsw.SetQueryTimeParams(query_time_params);
  }

  ~NMSLIB_hnsw()
  {
    for (auto p : m_points)
      delete p;
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
    m_hnsw.Search(&q, dummy);
    delete pt;

#ifdef PRINT_FOUND_NEIGHBORS
    std::cerr << "Query:\n";
    nms::KNNQueue<Dist_type> *res2 = q.Result()->Clone();
    while (!res2->Empty()) {
      std::cerr << "  " << res2->TopObject()->id() << " : " << res2->TopDistance()*res2->TopDistance() << "\n";
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

  int tree_depth() const
  {
    return -1;
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
  nms::ObjectVector m_points;
  Hnsw m_hnsw;
};

#endif // FUNCTOR_NMSLIB_HNSW_
#endif // GUDHI_NMSLIB_IS_AVAILABLE