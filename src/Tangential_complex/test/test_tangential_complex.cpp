/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Clement Jamin
 *
 *    Copyright (C) 2016  INRIA Sophia-Antipolis (France)
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

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Tangential_complex - test tangential complex
#include <boost/test/unit_test.hpp>

#include <gudhi/Tangential_complex.h>
#include <gudhi/sparsify_point_set.h>

#include <CGAL/Epick_d.h>
#include <CGAL/Random.h>

#include <array>
#include <vector>

namespace subsampl = Gudhi::subsampling;
namespace tc = Gudhi::tangential_complex;

BOOST_AUTO_TEST_CASE(test_Spatial_tree_data_structure)
{
  typedef CGAL::Epick_d<CGAL::Dynamic_dimension_tag>              Kernel;
  typedef Kernel::FT                                              FT;
  typedef Kernel::Point_d                                         Point;
  typedef Kernel::Vector_d                                        Vector;
  typedef tc::Tangential_complex<
    Kernel, CGAL::Dynamic_dimension_tag,
    CGAL::Parallel_tag>                                           TC;

  const int INTRINSIC_DIM = 2;
  const int AMBIENT_DIM = 3;
  const int NUM_POINTS = 50;

  Kernel k;

  // Generate points on a 2-sphere
  CGAL::Random_points_on_sphere_d<Point> generator(AMBIENT_DIM, 3.);
  std::vector<Point> points;
  points.reserve(NUM_POINTS);
  for (int i = 0; i < NUM_POINTS; ++i)
    points.push_back(*generator++);

  // Compute the TC
  TC tc(points, INTRINSIC_DIM, 0.01, k);
  tc.compute_tangential_complex();

  // Try to fix inconsistencies
  unsigned int num_perturb_steps;
  std::size_t initial_num_inconsistent_local_tr;
  std::size_t best_num_inconsistent_local_tr;
  std::size_t final_num_inconsistent_local_tr;
  auto perturb_ret = tc.fix_inconsistencies_using_perturbation(
    num_perturb_steps, initial_num_inconsistent_local_tr,
    best_num_inconsistent_local_tr, final_num_inconsistent_local_tr,
    60); // give it 60 seconds to succeed

  BOOST_CHECK(perturb_ret == TC_FIXED);

  // Export the TC into a Simplex_tree
  Gudhi::Simplex_tree<> stree;
  tc.export_complex(stree);

}
