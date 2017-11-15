/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Clement Jamin
 *
 *    Copyright (C) 2017 INRIA
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

#ifndef FUNCTOR_CGAL_KD_TREE_STORING_POINTS_
#define FUNCTOR_CGAL_KD_TREE_STORING_POINTS_

#include <CGAL/Search_traits.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>

#include <vector>
#include <utility>


// Recommended splitters:
// * CGAL::Sliding_midpoint: default in CGAL -- best choice for low dimension
// * CGAL::Median_of_max_spread: best choice for medium dimension
// * See CGAL doc for more...

template <typename Kernel, typename Splitter_ = CGAL::Median_of_max_spread<K> >
class CGAL_Kd_tree_storing_points
{
  typedef Kernel                                          K;
  typedef typename K::FT                                  FT;
  typedef typename K::Point_d                             Point;
  typedef std::vector<Point>                              Points;

  typedef K                                               STraits;

  typedef Splitter_                                       Splitter;
  typedef CGAL::Kd_tree<
    STraits, Splitter, CGAL::Tag_true>                    Tree;

  typedef CGAL::Orthogonal_k_neighbor_search<STraits,
    CGAL::Euclidean_distance<STraits>, Splitter, Tree>    Neighbor_search;

public:
  CGAL_Kd_tree_storing_points(Points const& points, double /*epsilon*/ = 0.)
    : m_tree(points.begin(), points.end())
  {}

  // Returns 0 because the indices are not returned in this version
  std::size_t query_k_nearest_neighbors(
    Point const& p,
    unsigned int k,
    double eps = 0.,
    std::vector<std::pair<std::size_t, double>> *result = NULL) const
  {
    Neighbor_search search(m_tree, p, k, eps);

#ifdef PRINT_FOUND_NEIGHBORS
    std::cerr << "Query:\n";
    for (auto nb : search)
      std::cerr << "  " << nb.first << " : " << nb.second << "\n";
#endif

    std::size_t sum = 0;
    if (result) {
      for (auto nb : search)
      {
        //sum += nb.first;
        result->push_back(std::make_pair(0, nb.second));
      }
    }
    else {
      //for (auto nb : search)
        //sum += nb.first;
    }
    return sum;
  }

  int tree_depth() const
  {
    return m_tree.root()->depth();
  }

private:
  Tree m_tree;
};

#endif // FUNCTOR_CGAL_KD_TREE_STORING_POINTS_
