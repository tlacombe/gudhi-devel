/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Siargey Kachanovich
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

// #ifdef _DEBUG
// # define TBB_USE_THREADING_TOOL
// #endif

// #define BOOST_TEST_DYN_LINK
// #define BOOST_TEST_MODULE "test_choose_farthest_point"
//#include <boost/test/unit_test.hpp>
//#include <boost/mpl/list.hpp>

#include <gudhi/choose_by_farthest_point.h>
#include <gudhi/feder_greene_clustering.h>
#include <vector>
#include <iterator>
#include <gudhi/Clock.h>

#include <CGAL/Epick_d.h>

typedef CGAL::Epick_d<CGAL::Dynamic_dimension_tag>                K;
typedef typename K::FT                                            FT;
typedef typename K::Point_d                                       Point_d;


//BOOST_AUTO_TEST_CASE(test_choose_farthest_point)
int main() {
  std::vector< Point_d > points, results, results2;
  //FT r1, r2;
  K k;
  Clock t;
  int width = 10;
  // Add grid points (810000 points)
  for (FT i = 0; i < width; i += 1.0)
    for (FT j = 0; j < width; j += 1.0)
      for (FT k = 0; k < width; k += 1.0)
        for (FT l = 0; l < width; l += 1.0)
          points.push_back(Point_d(std::vector<FT>({i, j, k, l})));

  std::cout << "Point set size: " << points.size() << std::endl;
  unsigned final_size = 5000, numeral = 1;
  std::cout << "Test   New     Old\n";
  while (final_size < 5001) {
    std::cout << final_size << ": ";
    results.clear();
    t.begin();
    Gudhi::subsampling::feder_greene_clustering(k, points, final_size, 0, std::back_inserter(results));
    t.end();
    std::cout << t.num_seconds()  << " s, ";
  
    // std::cout << "New algorithm result:\n";
    // for (auto p: results)
    //   std::cout << p << std::endl;
    //std:: cout << "r1^2 = " << r1 << std::endl;
    
    results2.clear();
    t.begin();
    Gudhi::subsampling::choose_by_farthest_point(k, points, final_size, 0, std::back_inserter(results2));
    t.end();
    std::cout << t.num_seconds()  << " s" << std::endl;
    
    // std::cout << "Old algorithm result:\n";
    // for (auto p: results2)
    //   std::cout << p << std::endl;
    //std::cout << "r2^2 = " << r2 << std::endl;
  
    // assert(results.size() == final_size);
    // assert(results2.size() == final_size);
    // assert(r1 == r2);

    switch (numeral) {
    case 1: numeral = 2; final_size *= 2; break;
    case 2: numeral = 5; final_size = final_size/2*5; break;
    case 5: numeral = 1; final_size *= 2; break;
    default: assert(false);
    }
    
  }
}
