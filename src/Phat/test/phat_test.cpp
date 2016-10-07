/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Pawel Dlotko
 *
 *    Copyright (C) 2016  INRIA Saclay (France)
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


#include <gudhi/reader_utils.h>

//cubical complex include
#include <gudhi/Bitmap_cubical_complex_base.h>
#include <gudhi/Bitmap_cubical_complex_periodic_boundary_conditions_base.h>
#include <gudhi/Bitmap_cubical_complex.h>
#include <gudhi/Compute_persistence_with_phat.h>

//Rips compelx and simpelx tree:
#include <gudhi/graph_simplicial_complex.h>
#include <gudhi/distance_functions.h>
#include <gudhi/Simplex_tree.h>
#include <boost/program_options.hpp>

//unit tests
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE "phat_bindings"
#include <boost/test/unit_test.hpp>


//standard stuff
#include <iostream>
#include <sstream>
#include <vector>
#include <string>
#include <limits>

using namespace Gudhi;
using namespace Gudhi::Cubical_complex;
using namespace Gudhi::phat_interface;
using namespace std;

BOOST_AUTO_TEST_CASE(phat_simplicial_twist) {
  bool dbg = false;
  double EPSILON = 0.000005;

  std::vector< std::pair<double, double> > dim0;
  dim0.push_back(std::make_pair(0, 0.0168145));
  dim0.push_back(std::make_pair(0, 0.0831613));
  dim0.push_back(std::make_pair(0, 0.0831613));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));

  std::vector< std::pair<double, double> > dim1;
  dim1.push_back(std::make_pair(0.0999583, 1.73485));
  std::vector< std::vector<double> > betti_numbers;
  std::vector<double> bn_dim_0;
  bn_dim_0.push_back(0);
  std::vector<double> bn_dim_1;
  std::vector<double> bn_dim_2;
  betti_numbers.push_back(bn_dim_0);
  betti_numbers.push_back(bn_dim_1);
  betti_numbers.push_back(bn_dim_2);


  std::string filepoints = "plane_circle";
  std::string filediag;
  double threshold = 2;
  int dim_max = 2;

  // Extract the points from the file filepoints
  typedef std::vector<double> Point_t;
  std::vector< Point_t > points;
  read_points(filepoints, points);

  // Compute the proximity graph of the points
  Graph_t prox_graph = compute_proximity_graph(points, threshold, euclidean_distance<Point_t>);

  // Construct the Rips complex in a Simplex Tree
  typedef Simplex_tree<Simplex_tree_options_fast_persistence> ST;
  ST st;
  // insert the proximity graph in the simplex tree
  st.insert_graph(prox_graph);
  // expand the graph until dimension dim_max
  st.expansion(dim_max);

  // Sort the simplices in the order of the filtration
  st.initialize_filtration();

  Compute_persistence_with_phat< ST > phat(&st);

  phat::persistence_pairs pairs = phat.compute_persistence_pairs_twist_reduction();
  std::pair< std::vector< std::vector<float> >, std::vector< std::vector< std::pair<float, float> > > > persistence = phat.get_the_intervals(pairs);




  //compare Betti numbers:
  if (dbg) {
    cerr << "betti_numbers.size() : " << betti_numbers.size() << endl;
    cerr << "persistence.first.size() : " << persistence.first.size() << endl;
  }
  BOOST_CHECK(betti_numbers.size() == persistence.first.size());
  for (size_t dim = 0; dim != persistence.first.size() - 1; ++dim) {
    if (dbg) {
      cerr << "betti_numbers[dim].size() : " << betti_numbers[dim].size() << endl;
      cerr << "persistence.first[dim].size() : " << persistence.first[dim].size() << endl;
    }
    BOOST_CHECK(betti_numbers[dim].size() == persistence.first[dim].size());
    for (size_t i = 0; i != persistence.first[dim].size(); ++i) {
      if (dbg) {
        cerr << "persistence.first[dim][i] : " << persistence.first[dim][i] << endl;
        cerr << "betti_numbers[dim][i] : " << betti_numbers[dim][i] << endl;
      }
      BOOST_CHECK(persistence.first[dim][i] == betti_numbers[dim][i]);
    }
  }



  //compare persistence :
  //first check if we get persistence in the same dimensions:
  BOOST_CHECK(persistence.second.size() == 2);
  //check dimension 0:
  BOOST_CHECK(persistence.second[0].size() == dim0.size());
  for (size_t i = 0; i != dim0.size(); ++i) {
    if (dbg) {
      cerr << "persistence.second[0][i] : " << persistence.second[0][i].first << " " << persistence.second[0][i].second << endl;
      cerr << "dim0[i] : " << dim0[i].first << " " << dim0[i].second << endl;
    }
    BOOST_CHECK((fabs(persistence.second[0][i].first - (double) dim0[i].first) < EPSILON));
    BOOST_CHECK((fabs(persistence.second[0][i].second - (double) dim0[i].second) < EPSILON));
  }
  //check dimension 1:
  BOOST_CHECK(persistence.second[1].size() == dim1.size());
  for (size_t i = 0; i != dim1.size(); ++i) {
    if (dbg) {
      cerr << "persistence.second[1][i] : " << persistence.second[1][i].first << " " << persistence.second[1][i].second << endl;
      cerr << "dim1[i] : " << dim1[i].first << " " << dim1[i].second << endl;

      cerr << "fabs(persistence.second[1][i].second - (double)dim1[i].second ) : " << fabs(persistence.second[1][i].second - (double) dim1[i].second) << endl;
    }
    BOOST_CHECK((fabs(persistence.second[1][i].first - (double) dim1[i].first) < EPSILON));
    BOOST_CHECK((fabs(persistence.second[1][i].second - (double) dim1[i].second) < EPSILON));
  }



}

BOOST_AUTO_TEST_CASE(phat_simplicial_chunk) {
  bool dbg = false;
  double EPSILON = 0.000005;

  std::vector< std::pair<double, double> > dim0;
  dim0.push_back(std::make_pair(0, 0.0168145));
  dim0.push_back(std::make_pair(0, 0.0831613));
  dim0.push_back(std::make_pair(0, 0.0831613));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));

  std::vector< std::pair<double, double> > dim1;
  dim1.push_back(std::make_pair(0.0999583, 1.73485));
  std::vector< std::vector<double> > betti_numbers;
  std::vector<double> bn_dim_0;
  bn_dim_0.push_back(0);
  std::vector<double> bn_dim_1;
  std::vector<double> bn_dim_2;
  betti_numbers.push_back(bn_dim_0);
  betti_numbers.push_back(bn_dim_1);
  betti_numbers.push_back(bn_dim_2);


  std::string filepoints = "plane_circle";
  std::string filediag;
  double threshold = 2;
  int dim_max = 2;

  // Extract the points from the file filepoints
  typedef std::vector<double> Point_t;
  std::vector< Point_t > points;
  read_points(filepoints, points);

  // Compute the proximity graph of the points
  Graph_t prox_graph = compute_proximity_graph(points, threshold, euclidean_distance<Point_t>);

  // Construct the Rips complex in a Simplex Tree
  typedef Simplex_tree<Simplex_tree_options_fast_persistence> ST;
  ST st;
  // insert the proximity graph in the simplex tree
  st.insert_graph(prox_graph);
  // expand the graph until dimension dim_max
  st.expansion(dim_max);

  // Sort the simplices in the order of the filtration
  st.initialize_filtration();

  Compute_persistence_with_phat< ST > phat(&st);

  phat::persistence_pairs pairs = phat.compute_persistence_pairs_dualized_chunk_reduction();
  std::pair< std::vector< std::vector<float> >, std::vector< std::vector< std::pair<float, float> > > > persistence = phat.get_the_intervals(pairs);




  //compare Betti numbers:
  if (dbg) {
    cerr << "betti_numbers.size() : " << betti_numbers.size() << endl;
    cerr << "persistence.first.size() : " << persistence.first.size() << endl;
  }
  BOOST_CHECK(betti_numbers.size() == persistence.first.size());
  for (size_t dim = 0; dim != persistence.first.size() - 1; ++dim) {
    if (dbg) {
      cerr << "betti_numbers[dim].size() : " << betti_numbers[dim].size() << endl;
      cerr << "persistence.first[dim].size() : " << persistence.first[dim].size() << endl;
    }
    BOOST_CHECK(betti_numbers[dim].size() == persistence.first[dim].size());
    for (size_t i = 0; i != persistence.first[dim].size(); ++i) {
      if (dbg) {
        cerr << "persistence.first[dim][i] : " << persistence.first[dim][i] << endl;
        cerr << "betti_numbers[dim][i] : " << betti_numbers[dim][i] << endl;
      }
      BOOST_CHECK(persistence.first[dim][i] == betti_numbers[dim][i]);
    }
  }



  //compare persistence :
  //first check if we get persistence in the same dimensions:
  BOOST_CHECK(persistence.second.size() == 2);
  //check dimension 0:
  BOOST_CHECK(persistence.second[0].size() == dim0.size());
  for (size_t i = 0; i != dim0.size(); ++i) {
    if (dbg) {
      cerr << "persistence.second[0][i] : " << persistence.second[0][i].first << " " << persistence.second[0][i].second << endl;
      cerr << "dim0[i] : " << dim0[i].first << " " << dim0[i].second << endl;
    }
    BOOST_CHECK((fabs(persistence.second[0][i].first - (double) dim0[i].first) < EPSILON));
    BOOST_CHECK((fabs(persistence.second[0][i].second - (double) dim0[i].second) < EPSILON));
  }
  //check dimension 1:
  BOOST_CHECK(persistence.second[1].size() == dim1.size());
  for (size_t i = 0; i != dim1.size(); ++i) {
    if (dbg) {
      cerr << "persistence.second[1][i] : " << persistence.second[1][i].first << " " << persistence.second[1][i].second << endl;
      cerr << "dim1[i] : " << dim1[i].first << " " << dim1[i].second << endl;

      cerr << "fabs(persistence.second[1][i].second - (double)dim1[i].second ) : " << fabs(persistence.second[1][i].second - (double) dim1[i].second) << endl;
    }
    BOOST_CHECK((fabs(persistence.second[1][i].first - (double) dim1[i].first) < EPSILON));
    BOOST_CHECK((fabs(persistence.second[1][i].second - (double) dim1[i].second) < EPSILON));
  }



}

BOOST_AUTO_TEST_CASE(phat_simplicial_standard_reductions) {
  bool dbg = false;
  double EPSILON = 0.000005;

  std::vector< std::pair<double, double> > dim0;
  dim0.push_back(std::make_pair(0, 0.0168145));
  dim0.push_back(std::make_pair(0, 0.0831613));
  dim0.push_back(std::make_pair(0, 0.0831613));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));

  std::vector< std::pair<double, double> > dim1;
  dim1.push_back(std::make_pair(0.0999583, 1.73485));
  std::vector< std::vector<double> > betti_numbers;
  std::vector<double> bn_dim_0;
  bn_dim_0.push_back(0);
  std::vector<double> bn_dim_1;
  std::vector<double> bn_dim_2;
  betti_numbers.push_back(bn_dim_0);
  betti_numbers.push_back(bn_dim_1);
  betti_numbers.push_back(bn_dim_2);


  std::string filepoints = "plane_circle";
  std::string filediag;
  double threshold = 2;
  int dim_max = 2;

  // Extract the points from the file filepoints
  typedef std::vector<double> Point_t;
  std::vector< Point_t > points;
  read_points(filepoints, points);

  // Compute the proximity graph of the points
  Graph_t prox_graph = compute_proximity_graph(points, threshold, euclidean_distance<Point_t>);

  // Construct the Rips complex in a Simplex Tree
  typedef Simplex_tree<Simplex_tree_options_fast_persistence> ST;
  ST st;
  // insert the proximity graph in the simplex tree
  st.insert_graph(prox_graph);
  // expand the graph until dimension dim_max
  st.expansion(dim_max);

  // Sort the simplices in the order of the filtration
  st.initialize_filtration();

  Compute_persistence_with_phat< ST > phat(&st);

  phat::persistence_pairs pairs = phat.compute_persistence_pairs_standard_reduction();
  std::pair< std::vector< std::vector<float> >, std::vector< std::vector< std::pair<float, float> > > > persistence = phat.get_the_intervals(pairs);




  //compare Betti numbers:
  if (dbg) {
    cerr << "betti_numbers.size() : " << betti_numbers.size() << endl;
    cerr << "persistence.first.size() : " << persistence.first.size() << endl;
  }
  BOOST_CHECK(betti_numbers.size() == persistence.first.size());
  for (size_t dim = 0; dim != persistence.first.size() - 1; ++dim) {
    if (dbg) {
      cerr << "betti_numbers[dim].size() : " << betti_numbers[dim].size() << endl;
      cerr << "persistence.first[dim].size() : " << persistence.first[dim].size() << endl;
    }
    BOOST_CHECK(betti_numbers[dim].size() == persistence.first[dim].size());
    for (size_t i = 0; i != persistence.first[dim].size(); ++i) {
      if (dbg) {
        cerr << "persistence.first[dim][i] : " << persistence.first[dim][i] << endl;
        cerr << "betti_numbers[dim][i] : " << betti_numbers[dim][i] << endl;
      }
      BOOST_CHECK(persistence.first[dim][i] == betti_numbers[dim][i]);
    }
  }



  //compare persistence :
  //first check if we get persistence in the same dimensions:
  BOOST_CHECK(persistence.second.size() == 2);
  //check dimension 0:
  BOOST_CHECK(persistence.second[0].size() == dim0.size());
  for (size_t i = 0; i != dim0.size(); ++i) {
    if (dbg) {
      cerr << "persistence.second[0][i] : " << persistence.second[0][i].first << " " << persistence.second[0][i].second << endl;
      cerr << "dim0[i] : " << dim0[i].first << " " << dim0[i].second << endl;
    }
    BOOST_CHECK((fabs(persistence.second[0][i].first - (double) dim0[i].first) < EPSILON));
    BOOST_CHECK((fabs(persistence.second[0][i].second - (double) dim0[i].second) < EPSILON));
  }
  //check dimension 1:
  BOOST_CHECK(persistence.second[1].size() == dim1.size());
  for (size_t i = 0; i != dim1.size(); ++i) {
    if (dbg) {
      cerr << "persistence.second[1][i] : " << persistence.second[1][i].first << " " << persistence.second[1][i].second << endl;
      cerr << "dim1[i] : " << dim1[i].first << " " << dim1[i].second << endl;

      cerr << "fabs(persistence.second[1][i].second - (double)dim1[i].second ) : " << fabs(persistence.second[1][i].second - (double) dim1[i].second) << endl;
    }
    BOOST_CHECK((fabs(persistence.second[1][i].first - (double) dim1[i].first) < EPSILON));
    BOOST_CHECK((fabs(persistence.second[1][i].second - (double) dim1[i].second) < EPSILON));
  }



}

BOOST_AUTO_TEST_CASE(phat_simplicial_spectral_sequence) {
  bool dbg = false;
  double EPSILON = 0.000005;

  std::vector< std::pair<double, double> > dim0;
  dim0.push_back(std::make_pair(0, 0.0168145));
  dim0.push_back(std::make_pair(0, 0.0831613));
  dim0.push_back(std::make_pair(0, 0.0831613));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));
  dim0.push_back(std::make_pair(0, 0.0999583));

  std::vector< std::pair<double, double> > dim1;
  dim1.push_back(std::make_pair(0.0999583, 1.73485));
  std::vector< std::vector<double> > betti_numbers;
  std::vector<double> bn_dim_0;
  bn_dim_0.push_back(0);
  std::vector<double> bn_dim_1;
  std::vector<double> bn_dim_2;
  betti_numbers.push_back(bn_dim_0);
  betti_numbers.push_back(bn_dim_1);
  betti_numbers.push_back(bn_dim_2);


  std::string filepoints = "plane_circle";
  std::string filediag;
  double threshold = 2;
  int dim_max = 2;

  // Extract the points from the file filepoints
  typedef std::vector<double> Point_t;
  std::vector< Point_t > points;
  read_points(filepoints, points);

  // Compute the proximity graph of the points
  Graph_t prox_graph = compute_proximity_graph(points, threshold, euclidean_distance<Point_t>);

  // Construct the Rips complex in a Simplex Tree
  typedef Simplex_tree<Simplex_tree_options_fast_persistence> ST;
  ST st;
  // insert the proximity graph in the simplex tree
  st.insert_graph(prox_graph);
  // expand the graph until dimension dim_max
  st.expansion(dim_max);

  // Sort the simplices in the order of the filtration
  st.initialize_filtration();

  Compute_persistence_with_phat< ST > phat(&st);

  phat::persistence_pairs pairs = phat.compute_persistence_pairs_spectral_sequence_reduction();
  std::pair< std::vector< std::vector<float> >, std::vector< std::vector< std::pair<float, float> > > > persistence = phat.get_the_intervals(pairs);




  //compare Betti numbers:
  if (dbg) {
    cerr << "betti_numbers.size() : " << betti_numbers.size() << endl;
    cerr << "persistence.first.size() : " << persistence.first.size() << endl;
  }
  BOOST_CHECK(betti_numbers.size() == persistence.first.size());
  for (size_t dim = 0; dim != persistence.first.size() - 1; ++dim) {
    if (dbg) {
      cerr << "betti_numbers[dim].size() : " << betti_numbers[dim].size() << endl;
      cerr << "persistence.first[dim].size() : " << persistence.first[dim].size() << endl;
    }
    BOOST_CHECK(betti_numbers[dim].size() == persistence.first[dim].size());
    for (size_t i = 0; i != persistence.first[dim].size(); ++i) {
      if (dbg) {
        cerr << "persistence.first[dim][i] : " << persistence.first[dim][i] << endl;
        cerr << "betti_numbers[dim][i] : " << betti_numbers[dim][i] << endl;
      }
      BOOST_CHECK(persistence.first[dim][i] == betti_numbers[dim][i]);
    }
  }



  //compare persistence :
  //first check if we get persistence in the same dimensions:
  BOOST_CHECK(persistence.second.size() == 2);
  //check dimension 0:
  BOOST_CHECK(persistence.second[0].size() == dim0.size());
  for (size_t i = 0; i != dim0.size(); ++i) {
    if (dbg) {
      cerr << "persistence.second[0][i] : " << persistence.second[0][i].first << " " << persistence.second[0][i].second << endl;
      cerr << "dim0[i] : " << dim0[i].first << " " << dim0[i].second << endl;
    }
    BOOST_CHECK((fabs(persistence.second[0][i].first - (double) dim0[i].first) < EPSILON));
    BOOST_CHECK((fabs(persistence.second[0][i].second - (double) dim0[i].second) < EPSILON));
  }
  //check dimension 1:
  BOOST_CHECK(persistence.second[1].size() == dim1.size());
  for (size_t i = 0; i != dim1.size(); ++i) {
    if (dbg) {
      cerr << "persistence.second[1][i] : " << persistence.second[1][i].first << " " << persistence.second[1][i].second << endl;
      cerr << "dim1[i] : " << dim1[i].first << " " << dim1[i].second << endl;

      cerr << "fabs(persistence.second[1][i].second - (double)dim1[i].second ) : " << fabs(persistence.second[1][i].second - (double) dim1[i].second) << endl;
    }
    BOOST_CHECK((fabs(persistence.second[1][i].first - (double) dim1[i].first) < EPSILON));
    BOOST_CHECK((fabs(persistence.second[1][i].second - (double) dim1[i].second) < EPSILON));
  }



}

BOOST_AUTO_TEST_CASE(phat_cubical_twist) {
  double EPSILON = 0.0000005;

  std::vector< std::pair<double, double> > dim0;
  dim0.push_back(std::make_pair(0.000356789, 0.0209834));
  dim0.push_back(std::make_pair(0.00137761, 0.178343));
  dim0.push_back(std::make_pair(0.00296197, 0.182727));
  dim0.push_back(std::make_pair(0.00406146, 0.138371));
  dim0.push_back(std::make_pair(0.00572209, 0.0095895));
  dim0.push_back(std::make_pair(0.00724978, 0.088271));
  dim0.push_back(std::make_pair(0.00736106, 0.138713));
  dim0.push_back(std::make_pair(0.00779113, 0.176245));
  dim0.push_back(std::make_pair(0.00841857, 0.088271));
  dim0.push_back(std::make_pair(0.0101254, 0.0272611));
  dim0.push_back(std::make_pair(0.0133762, 0.11487));
  dim0.push_back(std::make_pair(0.0134719, 0.116554));
  dim0.push_back(std::make_pair(0.0145639, 0.155647));
  dim0.push_back(std::make_pair(0.0157889, 0.0683329));
  dim0.push_back(std::make_pair(0.0161422, 0.0476644));
  dim0.push_back(std::make_pair(0.0173804, 0.0743558));
  dim0.push_back(std::make_pair(0.0184607, 0.25154));
  dim0.push_back(std::make_pair(0.0186499, 0.339905));
  dim0.push_back(std::make_pair(0.0247987, 0.0559076));
  dim0.push_back(std::make_pair(0.0250133, 0.10707));
  dim0.push_back(std::make_pair(0.0251446, 0.222457));
  dim0.push_back(std::make_pair(0.026607, 0.12124));
  dim0.push_back(std::make_pair(0.0279297, 0.226968));
  dim0.push_back(std::make_pair(0.0294929, 0.0743558));
  dim0.push_back(std::make_pair(0.0321475, 0.0959911));
  dim0.push_back(std::make_pair(0.0338132, 0.15912));
  dim0.push_back(std::make_pair(0.0345952, 0.0956828));
  dim0.push_back(std::make_pair(0.0350849, 0.0476644));
  dim0.push_back(std::make_pair(0.0354295, 0.138721));
  dim0.push_back(std::make_pair(0.0358397, 0.138721));
  dim0.push_back(std::make_pair(0.0369627, 0.114244));
  dim0.push_back(std::make_pair(0.0382787, 0.153073));
  dim0.push_back(std::make_pair(0.0389954, 0.2347));
  dim0.push_back(std::make_pair(0.0424284, 0.0873569));
  dim0.push_back(std::make_pair(0.0460244, 0.0628887));
  dim0.push_back(std::make_pair(0.0460566, 0.100457));
  dim0.push_back(std::make_pair(0.0548351, 0.0956828));
  dim0.push_back(std::make_pair(0.0568281, 0.108871));
  dim0.push_back(std::make_pair(0.0598285, 0.0894398));
  dim0.push_back(std::make_pair(0.0643188, 0.180772));
  dim0.push_back(std::make_pair(0.067139, 0.217477));
  dim0.push_back(std::make_pair(0.0685163, 0.149041));
  dim0.push_back(std::make_pair(0.069365, 0.0840998));
  dim0.push_back(std::make_pair(0.0725224, 0.0873569));
  dim0.push_back(std::make_pair(0.0726488, 0.141462));
  dim0.push_back(std::make_pair(0.0750035, 0.156528));
  dim0.push_back(std::make_pair(0.0784113, 0.136386));
  dim0.push_back(std::make_pair(0.082815, 0.211298));
  dim0.push_back(std::make_pair(0.0857338, 0.186306));
  dim0.push_back(std::make_pair(0.0944081, 0.133989));
  dim0.push_back(std::make_pair(0.0958755, 0.238023));
  dim0.push_back(std::make_pair(0.114808, 0.124519));
  dim0.push_back(std::make_pair(0.142382, 0.20245));
  dim0.push_back(std::make_pair(0.147678, 0.161484));
  dim0.push_back(std::make_pair(0.158924, 0.217477));

  std::vector< std::pair<double, double> > dim1;
  dim1.push_back(std::make_pair(0.0896383, 0.307383));
  dim1.push_back(std::make_pair(0.0896383, 0.574247));
  dim1.push_back(std::make_pair(0.0922039, 0.337222));
  dim1.push_back(std::make_pair(0.0952979, 0.635619));
  dim1.push_back(std::make_pair(0.10707, 0.65752));
  dim1.push_back(std::make_pair(0.112554, 0.437864));
  dim1.push_back(std::make_pair(0.114244, 0.322239));
  dim1.push_back(std::make_pair(0.115086, 0.31284));
  dim1.push_back(std::make_pair(0.131682, 0.609287));
  dim1.push_back(std::make_pair(0.138669, 0.15837));
  dim1.push_back(std::make_pair(0.138669, 0.212147));
  dim1.push_back(std::make_pair(0.138669, 0.299792));
  dim1.push_back(std::make_pair(0.146283, 0.375957));
  dim1.push_back(std::make_pair(0.150632, 0.354247));
  dim1.push_back(std::make_pair(0.150632, 0.597383));
  dim1.push_back(std::make_pair(0.152943, 0.283281));
  dim1.push_back(std::make_pair(0.15912, 0.642942));
  dim1.push_back(std::make_pair(0.160004, 0.683229));
  dim1.push_back(std::make_pair(0.165983, 0.287496));
  dim1.push_back(std::make_pair(0.165983, 0.307383));
  dim1.push_back(std::make_pair(0.171803, 0.190303));
  dim1.push_back(std::make_pair(0.177224, 0.559145));
  dim1.push_back(std::make_pair(0.177959, 0.219433));
  dim1.push_back(std::make_pair(0.177959, 0.481775));
  dim1.push_back(std::make_pair(0.179939, 0.426291));
  dim1.push_back(std::make_pair(0.180772, 0.600957));
  dim1.push_back(std::make_pair(0.182727, 0.593896));
  dim1.push_back(std::make_pair(0.184009, 0.504145));
  dim1.push_back(std::make_pair(0.193217, 0.714115));
  dim1.push_back(std::make_pair(0.20074, 0.26084));
  dim1.push_back(std::make_pair(0.203173, 0.493206));
  dim1.push_back(std::make_pair(0.203173, 0.538447));
  dim1.push_back(std::make_pair(0.203361, 0.618832));
  dim1.push_back(std::make_pair(0.203361, 0.656591));
  dim1.push_back(std::make_pair(0.205348, 0.427394));
  dim1.push_back(std::make_pair(0.205348, 0.556077));
  dim1.push_back(std::make_pair(0.213038, 0.642942));
  dim1.push_back(std::make_pair(0.21696, 0.397078));
  dim1.push_back(std::make_pair(0.21696, 0.493206));
  dim1.push_back(std::make_pair(0.218123, 0.812171));
  dim1.push_back(std::make_pair(0.222869, 0.249479));
  dim1.push_back(std::make_pair(0.2241, 0.692066));
  dim1.push_back(std::make_pair(0.224358, 0.667901));
  dim1.push_back(std::make_pair(0.224437, 0.254382));
  dim1.push_back(std::make_pair(0.224437, 0.393736));
  dim1.push_back(std::make_pair(0.224437, 0.562944));
  dim1.push_back(std::make_pair(0.226809, 0.342489));
  dim1.push_back(std::make_pair(0.227838, 0.58166));
  dim1.push_back(std::make_pair(0.227838, 0.671021));
  dim1.push_back(std::make_pair(0.229861, 0.259465));
  dim1.push_back(std::make_pair(0.233203, 0.504145));
  dim1.push_back(std::make_pair(0.2347, 0.253456));
  dim1.push_back(std::make_pair(0.2347, 0.401604));
  dim1.push_back(std::make_pair(0.23542, 0.509011));
  dim1.push_back(std::make_pair(0.23542, 0.656926));
  dim1.push_back(std::make_pair(0.238023, 0.306904));
  dim1.push_back(std::make_pair(0.238187, 0.560137));
  dim1.push_back(std::make_pair(0.242099, 0.375768));
  dim1.push_back(std::make_pair(0.251213, 0.419192));
  dim1.push_back(std::make_pair(0.251213, 0.470896));
  dim1.push_back(std::make_pair(0.25154, 0.441861));
  dim1.push_back(std::make_pair(0.252778, 0.683845));
  dim1.push_back(std::make_pair(0.253932, 0.303766));
  dim1.push_back(std::make_pair(0.25521, 0.318305));
  dim1.push_back(std::make_pair(0.25521, 0.366637));
  dim1.push_back(std::make_pair(0.256709, 0.552841));
  dim1.push_back(std::make_pair(0.256709, 0.566244));
  dim1.push_back(std::make_pair(0.256923, 0.662304));
  dim1.push_back(std::make_pair(0.261521, 0.634892));
  dim1.push_back(std::make_pair(0.270009, 0.642002));
  dim1.push_back(std::make_pair(0.271325, 0.329129));
  dim1.push_back(std::make_pair(0.27284, 0.801019));
  dim1.push_back(std::make_pair(0.273637, 0.283281));
  dim1.push_back(std::make_pair(0.273637, 0.505287));
  dim1.push_back(std::make_pair(0.276398, 0.667667));
  dim1.push_back(std::make_pair(0.279702, 0.365257));
  dim1.push_back(std::make_pair(0.280491, 0.579041));
  dim1.push_back(std::make_pair(0.284502, 0.405565));
  dim1.push_back(std::make_pair(0.284932, 0.377696));
  dim1.push_back(std::make_pair(0.291091, 0.624392));
  dim1.push_back(std::make_pair(0.291315, 0.443989));
  dim1.push_back(std::make_pair(0.293683, 0.522518));
  dim1.push_back(std::make_pair(0.293683, 0.741924));
  dim1.push_back(std::make_pair(0.307746, 0.588104));
  dim1.push_back(std::make_pair(0.318734, 0.498793));
  dim1.push_back(std::make_pair(0.32605, 0.52335));
  dim1.push_back(std::make_pair(0.335347, 0.58166));
  dim1.push_back(std::make_pair(0.336204, 0.71353));
  dim1.push_back(std::make_pair(0.342399, 0.595119));
  dim1.push_back(std::make_pair(0.348409, 0.664241));
  dim1.push_back(std::make_pair(0.357035, 0.70087));
  dim1.push_back(std::make_pair(0.357874, 0.476921));
  dim1.push_back(std::make_pair(0.36049, 0.640687));
  dim1.push_back(std::make_pair(0.363861, 0.494538));
  dim1.push_back(std::make_pair(0.363861, 0.494538));
  dim1.push_back(std::make_pair(0.364503, 0.640687));
  dim1.push_back(std::make_pair(0.364503, 0.640687));
  dim1.push_back(std::make_pair(0.364571, 0.661572));
  dim1.push_back(std::make_pair(0.367343, 0.597987));
  dim1.push_back(std::make_pair(0.387655, 0.569132));
  dim1.push_back(std::make_pair(0.397369, 0.644901));
  dim1.push_back(std::make_pair(0.397369, 0.659947));
  dim1.push_back(std::make_pair(0.397386, 0.431372));
  dim1.push_back(std::make_pair(0.39905, 0.440457));
  dim1.push_back(std::make_pair(0.399317, 0.566299));
  dim1.push_back(std::make_pair(0.4005, 0.406971));
  dim1.push_back(std::make_pair(0.400725, 0.484885));
  dim1.push_back(std::make_pair(0.405866, 0.426268));
  dim1.push_back(std::make_pair(0.407038, 0.478693));
  dim1.push_back(std::make_pair(0.409665, 0.429602));
  dim1.push_back(std::make_pair(0.418432, 0.476169));
  dim1.push_back(std::make_pair(0.420472, 0.478693));
  dim1.push_back(std::make_pair(0.424782, 0.626694));
  dim1.push_back(std::make_pair(0.424782, 0.629479));
  dim1.push_back(std::make_pair(0.438875, 0.504145));
  dim1.push_back(std::make_pair(0.439381, 0.539878));
  dim1.push_back(std::make_pair(0.439816, 0.527948));
  dim1.push_back(std::make_pair(0.442226, 0.529129));
  dim1.push_back(std::make_pair(0.446357, 0.454291));
  dim1.push_back(std::make_pair(0.452574, 0.527445));
  dim1.push_back(std::make_pair(0.457876, 0.504913));
  dim1.push_back(std::make_pair(0.476169, 0.529129));
  dim1.push_back(std::make_pair(0.479762, 0.51536));
  dim1.push_back(std::make_pair(0.479762, 0.528202));
  dim1.push_back(std::make_pair(0.479762, 0.612093));
  dim1.push_back(std::make_pair(0.499352, 0.552841));
  dim1.push_back(std::make_pair(0.515218, 0.579041));
  dim1.push_back(std::make_pair(0.515218, 0.696697));
  dim1.push_back(std::make_pair(0.51536, 0.661572));
  dim1.push_back(std::make_pair(0.516364, 0.579478));
  dim1.push_back(std::make_pair(0.524837, 0.747096));
  dim1.push_back(std::make_pair(0.530507, 0.611131));
  dim1.push_back(std::make_pair(0.540148, 0.560137));
  dim1.push_back(std::make_pair(0.544518, 0.547481));
  dim1.push_back(std::make_pair(0.547481, 0.595119));
  dim1.push_back(std::make_pair(0.552251, 0.635619));
  dim1.push_back(std::make_pair(0.567289, 0.595748));
  dim1.push_back(std::make_pair(0.571876, 0.593348));
  dim1.push_back(std::make_pair(0.597383, 0.698937));
  dim1.push_back(std::make_pair(0.629459, 0.65752));
  dim1.push_back(std::make_pair(0.634643, 0.667667));
  dim1.push_back(std::make_pair(0.634892, 0.667901));


  std::vector< std::pair<double, double> > dim2;
  dim2.push_back(std::make_pair(0.159766, 0.514062));
  dim2.push_back(std::make_pair(0.397435, 0.857648));
  dim2.push_back(std::make_pair(0.424396, 0.494389));
  dim2.push_back(std::make_pair(0.426268, 0.831535));
  dim2.push_back(std::make_pair(0.443989, 0.737309));
  dim2.push_back(std::make_pair(0.521027, 0.55268));
  dim2.push_back(std::make_pair(0.5294, 0.827161));
  dim2.push_back(std::make_pair(0.531269, 0.975176));
  dim2.push_back(std::make_pair(0.538447, 0.893562));
  dim2.push_back(std::make_pair(0.559145, 0.78499));
  dim2.push_back(std::make_pair(0.560137, 0.764799));
  dim2.push_back(std::make_pair(0.562944, 0.672131));
  dim2.push_back(std::make_pair(0.568065, 0.833694));
  dim2.push_back(std::make_pair(0.582169, 0.614194));
  dim2.push_back(std::make_pair(0.595119, 0.99503));
  dim2.push_back(std::make_pair(0.6178, 0.645186));
  dim2.push_back(std::make_pair(0.633861, 0.894727));
  dim2.push_back(std::make_pair(0.65058, 0.720305));
  dim2.push_back(std::make_pair(0.654454, 0.882628));
  dim2.push_back(std::make_pair(0.656591, 0.911192));
  dim2.push_back(std::make_pair(0.663122, 0.835692));
  dim2.push_back(std::make_pair(0.677699, 0.828483));
  dim2.push_back(std::make_pair(0.684151, 0.714499));
  dim2.push_back(std::make_pair(0.691498, 0.79696));
  dim2.push_back(std::make_pair(0.692066, 0.836093));
  dim2.push_back(std::make_pair(0.696697, 0.924678));
  dim2.push_back(std::make_pair(0.698937, 0.888513));
  dim2.push_back(std::make_pair(0.708468, 0.825471));
  dim2.push_back(std::make_pair(0.714853, 0.998385));
  dim2.push_back(std::make_pair(0.723869, 0.972154));
  dim2.push_back(std::make_pair(0.725867, 0.900344));
  dim2.push_back(std::make_pair(0.726513, 0.972148));
  dim2.push_back(std::make_pair(0.728115, 0.853808));
  dim2.push_back(std::make_pair(0.730914, 0.95898));
  dim2.push_back(std::make_pair(0.734646, 0.940408));
  dim2.push_back(std::make_pair(0.735457, 0.784873));
  dim2.push_back(std::make_pair(0.737998, 0.989994));
  dim2.push_back(std::make_pair(0.749057, 0.796175));
  dim2.push_back(std::make_pair(0.755644, 0.760529));
  dim2.push_back(std::make_pair(0.763753, 0.770437));
  dim2.push_back(std::make_pair(0.775107, 0.801097));
  dim2.push_back(std::make_pair(0.776465, 0.808596));
  dim2.push_back(std::make_pair(0.78092, 0.836048));
  dim2.push_back(std::make_pair(0.784286, 0.957797));
  dim2.push_back(std::make_pair(0.792679, 0.834564));
  dim2.push_back(std::make_pair(0.792679, 0.940152));
  dim2.push_back(std::make_pair(0.794893, 0.919911));
  dim2.push_back(std::make_pair(0.808823, 0.986678));
  dim2.push_back(std::make_pair(0.813632, 0.914683));
  dim2.push_back(std::make_pair(0.823974, 0.96578));
  dim2.push_back(std::make_pair(0.826181, 0.952838));
  dim2.push_back(std::make_pair(0.839664, 0.932946));
  dim2.push_back(std::make_pair(0.843357, 0.921802));
  dim2.push_back(std::make_pair(0.84446, 0.950299));
  dim2.push_back(std::make_pair(0.84446, 0.966199));
  dim2.push_back(std::make_pair(0.848269, 0.937031));
  dim2.push_back(std::make_pair(0.856414, 0.911481));
  dim2.push_back(std::make_pair(0.882767, 0.919131));
  dim2.push_back(std::make_pair(0.882767, 0.939881));
  dim2.push_back(std::make_pair(0.889066, 0.972132));
  dim2.push_back(std::make_pair(0.891205, 0.917646));
  dim2.push_back(std::make_pair(0.891205, 0.920893));
  dim2.push_back(std::make_pair(0.895655, 0.901842));
  dim2.push_back(std::make_pair(0.911343, 0.930285));
  dim2.push_back(std::make_pair(0.916196, 0.933545));
  dim2.push_back(std::make_pair(0.93033, 0.950769));
  dim2.push_back(std::make_pair(0.936075, 0.936608));


  std::vector< std::vector<double> > betti_numbers;
  std::vector<double> bn_dim_0;
  bn_dim_0.push_back(0.000103006);
  std::vector<double> bn_dim_1;
  std::vector<double> bn_dim_2;
  std::vector<double> bn_dim_3;
  betti_numbers.push_back(bn_dim_0);
  betti_numbers.push_back(bn_dim_1);
  betti_numbers.push_back(bn_dim_2);
  betti_numbers.push_back(bn_dim_3);



  Bitmap_cubical_complex< Bitmap_cubical_complex_base<double> > b("random_cubical_complex");
  Compute_persistence_with_phat< Bitmap_cubical_complex< Bitmap_cubical_complex_base<double> > > phat(&b);
  phat::persistence_pairs pairs = phat.compute_persistence_pairs_twist_reduction();
  std::pair< std::vector< std::vector<double> >, std::vector< std::vector< std::pair<double, double> > > > persistence = phat.get_the_intervals(pairs);



  //compare Betti numbers:  
  BOOST_CHECK(betti_numbers.size() == persistence.first.size());
  for (size_t dim = 0; dim != persistence.first.size() - 1; ++dim) {
    BOOST_CHECK(betti_numbers[dim].size() == persistence.first[dim].size());

    for (size_t i = 0; i != persistence.first[dim].size(); ++i) {
      BOOST_CHECK(fabs(persistence.first[dim][i] - betti_numbers[dim][i]) < EPSILON);
    }

  }



  //compare persistence :
  //first check if we get persistence in the same dimensions:
  BOOST_CHECK(persistence.second.size() == 3);

  //check dimension 0:
  BOOST_CHECK(persistence.second[0].size() == dim0.size());
  for (size_t i = 0; i != persistence.second[0].size(); ++i) {
    BOOST_CHECK((fabs(persistence.second[0][i].first - (double) dim0[i].first) < EPSILON));
    BOOST_CHECK((fabs(persistence.second[0][i].second - (double) dim0[i].second) < EPSILON));
  }

  //check dimension 1:
  BOOST_CHECK(persistence.second[1].size() == dim1.size());
  for (size_t i = 0; i != persistence.second[1].size(); ++i) {
    BOOST_CHECK((fabs(persistence.second[1][i].first - (double) dim1[i].first) < EPSILON));
    BOOST_CHECK((fabs(persistence.second[1][i].second - (double) dim1[i].second) < EPSILON));
  }

  //check dimension 2:
  BOOST_CHECK(persistence.second[2].size() == dim2.size());
  for (size_t i = 0; i != persistence.second[2].size(); ++i) {
    BOOST_CHECK((fabs(persistence.second[2][i].first - (double) dim2[i].first) < EPSILON));
    BOOST_CHECK((fabs(persistence.second[2][i].second - (double) dim2[i].second) < EPSILON));
  }
}

BOOST_AUTO_TEST_CASE(phat_cubical_chunk) {
  double EPSILON = 0.0000005;

  std::vector< std::pair<double, double> > dim0;
  dim0.push_back(std::make_pair(0.000356789, 0.0209834));
  dim0.push_back(std::make_pair(0.00137761, 0.178343));
  dim0.push_back(std::make_pair(0.00296197, 0.182727));
  dim0.push_back(std::make_pair(0.00406146, 0.138371));
  dim0.push_back(std::make_pair(0.00572209, 0.0095895));
  dim0.push_back(std::make_pair(0.00724978, 0.088271));
  dim0.push_back(std::make_pair(0.00736106, 0.138713));
  dim0.push_back(std::make_pair(0.00779113, 0.176245));
  dim0.push_back(std::make_pair(0.00841857, 0.088271));
  dim0.push_back(std::make_pair(0.0101254, 0.0272611));
  dim0.push_back(std::make_pair(0.0133762, 0.11487));
  dim0.push_back(std::make_pair(0.0134719, 0.116554));
  dim0.push_back(std::make_pair(0.0145639, 0.155647));
  dim0.push_back(std::make_pair(0.0157889, 0.0683329));
  dim0.push_back(std::make_pair(0.0161422, 0.0476644));
  dim0.push_back(std::make_pair(0.0173804, 0.0743558));
  dim0.push_back(std::make_pair(0.0184607, 0.25154));
  dim0.push_back(std::make_pair(0.0186499, 0.339905));
  dim0.push_back(std::make_pair(0.0247987, 0.0559076));
  dim0.push_back(std::make_pair(0.0250133, 0.10707));
  dim0.push_back(std::make_pair(0.0251446, 0.222457));
  dim0.push_back(std::make_pair(0.026607, 0.12124));
  dim0.push_back(std::make_pair(0.0279297, 0.226968));
  dim0.push_back(std::make_pair(0.0294929, 0.0743558));
  dim0.push_back(std::make_pair(0.0321475, 0.0959911));
  dim0.push_back(std::make_pair(0.0338132, 0.15912));
  dim0.push_back(std::make_pair(0.0345952, 0.0956828));
  dim0.push_back(std::make_pair(0.0350849, 0.0476644));
  dim0.push_back(std::make_pair(0.0354295, 0.138721));
  dim0.push_back(std::make_pair(0.0358397, 0.138721));
  dim0.push_back(std::make_pair(0.0369627, 0.114244));
  dim0.push_back(std::make_pair(0.0382787, 0.153073));
  dim0.push_back(std::make_pair(0.0389954, 0.2347));
  dim0.push_back(std::make_pair(0.0424284, 0.0873569));
  dim0.push_back(std::make_pair(0.0460244, 0.0628887));
  dim0.push_back(std::make_pair(0.0460566, 0.100457));
  dim0.push_back(std::make_pair(0.0548351, 0.0956828));
  dim0.push_back(std::make_pair(0.0568281, 0.108871));
  dim0.push_back(std::make_pair(0.0598285, 0.0894398));
  dim0.push_back(std::make_pair(0.0643188, 0.180772));
  dim0.push_back(std::make_pair(0.067139, 0.217477));
  dim0.push_back(std::make_pair(0.0685163, 0.149041));
  dim0.push_back(std::make_pair(0.069365, 0.0840998));
  dim0.push_back(std::make_pair(0.0725224, 0.0873569));
  dim0.push_back(std::make_pair(0.0726488, 0.141462));
  dim0.push_back(std::make_pair(0.0750035, 0.156528));
  dim0.push_back(std::make_pair(0.0784113, 0.136386));
  dim0.push_back(std::make_pair(0.082815, 0.211298));
  dim0.push_back(std::make_pair(0.0857338, 0.186306));
  dim0.push_back(std::make_pair(0.0944081, 0.133989));
  dim0.push_back(std::make_pair(0.0958755, 0.238023));
  dim0.push_back(std::make_pair(0.114808, 0.124519));
  dim0.push_back(std::make_pair(0.142382, 0.20245));
  dim0.push_back(std::make_pair(0.147678, 0.161484));
  dim0.push_back(std::make_pair(0.158924, 0.217477));

  std::vector< std::pair<double, double> > dim1;
  dim1.push_back(std::make_pair(0.0896383, 0.307383));
  dim1.push_back(std::make_pair(0.0896383, 0.574247));
  dim1.push_back(std::make_pair(0.0922039, 0.337222));
  dim1.push_back(std::make_pair(0.0952979, 0.635619));
  dim1.push_back(std::make_pair(0.10707, 0.65752));
  dim1.push_back(std::make_pair(0.112554, 0.437864));
  dim1.push_back(std::make_pair(0.114244, 0.322239));
  dim1.push_back(std::make_pair(0.115086, 0.31284));
  dim1.push_back(std::make_pair(0.131682, 0.609287));
  dim1.push_back(std::make_pair(0.138669, 0.15837));
  dim1.push_back(std::make_pair(0.138669, 0.212147));
  dim1.push_back(std::make_pair(0.138669, 0.299792));
  dim1.push_back(std::make_pair(0.146283, 0.375957));
  dim1.push_back(std::make_pair(0.150632, 0.354247));
  dim1.push_back(std::make_pair(0.150632, 0.597383));
  dim1.push_back(std::make_pair(0.152943, 0.283281));
  dim1.push_back(std::make_pair(0.15912, 0.642942));
  dim1.push_back(std::make_pair(0.160004, 0.683229));
  dim1.push_back(std::make_pair(0.165983, 0.287496));
  dim1.push_back(std::make_pair(0.165983, 0.307383));
  dim1.push_back(std::make_pair(0.171803, 0.190303));
  dim1.push_back(std::make_pair(0.177224, 0.559145));
  dim1.push_back(std::make_pair(0.177959, 0.219433));
  dim1.push_back(std::make_pair(0.177959, 0.481775));
  dim1.push_back(std::make_pair(0.179939, 0.426291));
  dim1.push_back(std::make_pair(0.180772, 0.600957));
  dim1.push_back(std::make_pair(0.182727, 0.593896));
  dim1.push_back(std::make_pair(0.184009, 0.504145));
  dim1.push_back(std::make_pair(0.193217, 0.714115));
  dim1.push_back(std::make_pair(0.20074, 0.26084));
  dim1.push_back(std::make_pair(0.203173, 0.493206));
  dim1.push_back(std::make_pair(0.203173, 0.538447));
  dim1.push_back(std::make_pair(0.203361, 0.618832));
  dim1.push_back(std::make_pair(0.203361, 0.656591));
  dim1.push_back(std::make_pair(0.205348, 0.427394));
  dim1.push_back(std::make_pair(0.205348, 0.556077));
  dim1.push_back(std::make_pair(0.213038, 0.642942));
  dim1.push_back(std::make_pair(0.21696, 0.397078));
  dim1.push_back(std::make_pair(0.21696, 0.493206));
  dim1.push_back(std::make_pair(0.218123, 0.812171));
  dim1.push_back(std::make_pair(0.222869, 0.249479));
  dim1.push_back(std::make_pair(0.2241, 0.692066));
  dim1.push_back(std::make_pair(0.224358, 0.667901));
  dim1.push_back(std::make_pair(0.224437, 0.254382));
  dim1.push_back(std::make_pair(0.224437, 0.393736));
  dim1.push_back(std::make_pair(0.224437, 0.562944));
  dim1.push_back(std::make_pair(0.226809, 0.342489));
  dim1.push_back(std::make_pair(0.227838, 0.58166));
  dim1.push_back(std::make_pair(0.227838, 0.671021));
  dim1.push_back(std::make_pair(0.229861, 0.259465));
  dim1.push_back(std::make_pair(0.233203, 0.504145));
  dim1.push_back(std::make_pair(0.2347, 0.253456));
  dim1.push_back(std::make_pair(0.2347, 0.401604));
  dim1.push_back(std::make_pair(0.23542, 0.509011));
  dim1.push_back(std::make_pair(0.23542, 0.656926));
  dim1.push_back(std::make_pair(0.238023, 0.306904));
  dim1.push_back(std::make_pair(0.238187, 0.560137));
  dim1.push_back(std::make_pair(0.242099, 0.375768));
  dim1.push_back(std::make_pair(0.251213, 0.419192));
  dim1.push_back(std::make_pair(0.251213, 0.470896));
  dim1.push_back(std::make_pair(0.25154, 0.441861));
  dim1.push_back(std::make_pair(0.252778, 0.683845));
  dim1.push_back(std::make_pair(0.253932, 0.303766));
  dim1.push_back(std::make_pair(0.25521, 0.318305));
  dim1.push_back(std::make_pair(0.25521, 0.366637));
  dim1.push_back(std::make_pair(0.256709, 0.552841));
  dim1.push_back(std::make_pair(0.256709, 0.566244));
  dim1.push_back(std::make_pair(0.256923, 0.662304));
  dim1.push_back(std::make_pair(0.261521, 0.634892));
  dim1.push_back(std::make_pair(0.270009, 0.642002));
  dim1.push_back(std::make_pair(0.271325, 0.329129));
  dim1.push_back(std::make_pair(0.27284, 0.801019));
  dim1.push_back(std::make_pair(0.273637, 0.283281));
  dim1.push_back(std::make_pair(0.273637, 0.505287));
  dim1.push_back(std::make_pair(0.276398, 0.667667));
  dim1.push_back(std::make_pair(0.279702, 0.365257));
  dim1.push_back(std::make_pair(0.280491, 0.579041));
  dim1.push_back(std::make_pair(0.284502, 0.405565));
  dim1.push_back(std::make_pair(0.284932, 0.377696));
  dim1.push_back(std::make_pair(0.291091, 0.624392));
  dim1.push_back(std::make_pair(0.291315, 0.443989));
  dim1.push_back(std::make_pair(0.293683, 0.522518));
  dim1.push_back(std::make_pair(0.293683, 0.741924));
  dim1.push_back(std::make_pair(0.307746, 0.588104));
  dim1.push_back(std::make_pair(0.318734, 0.498793));
  dim1.push_back(std::make_pair(0.32605, 0.52335));
  dim1.push_back(std::make_pair(0.335347, 0.58166));
  dim1.push_back(std::make_pair(0.336204, 0.71353));
  dim1.push_back(std::make_pair(0.342399, 0.595119));
  dim1.push_back(std::make_pair(0.348409, 0.664241));
  dim1.push_back(std::make_pair(0.357035, 0.70087));
  dim1.push_back(std::make_pair(0.357874, 0.476921));
  dim1.push_back(std::make_pair(0.36049, 0.640687));
  dim1.push_back(std::make_pair(0.363861, 0.494538));
  dim1.push_back(std::make_pair(0.363861, 0.494538));
  dim1.push_back(std::make_pair(0.364503, 0.640687));
  dim1.push_back(std::make_pair(0.364503, 0.640687));
  dim1.push_back(std::make_pair(0.364571, 0.661572));
  dim1.push_back(std::make_pair(0.367343, 0.597987));
  dim1.push_back(std::make_pair(0.387655, 0.569132));
  dim1.push_back(std::make_pair(0.397369, 0.644901));
  dim1.push_back(std::make_pair(0.397369, 0.659947));
  dim1.push_back(std::make_pair(0.397386, 0.431372));
  dim1.push_back(std::make_pair(0.39905, 0.440457));
  dim1.push_back(std::make_pair(0.399317, 0.566299));
  dim1.push_back(std::make_pair(0.4005, 0.406971));
  dim1.push_back(std::make_pair(0.400725, 0.484885));
  dim1.push_back(std::make_pair(0.405866, 0.426268));
  dim1.push_back(std::make_pair(0.407038, 0.478693));
  dim1.push_back(std::make_pair(0.409665, 0.429602));
  dim1.push_back(std::make_pair(0.418432, 0.476169));
  dim1.push_back(std::make_pair(0.420472, 0.478693));
  dim1.push_back(std::make_pair(0.424782, 0.626694));
  dim1.push_back(std::make_pair(0.424782, 0.629479));
  dim1.push_back(std::make_pair(0.438875, 0.504145));
  dim1.push_back(std::make_pair(0.439381, 0.539878));
  dim1.push_back(std::make_pair(0.439816, 0.527948));
  dim1.push_back(std::make_pair(0.442226, 0.529129));
  dim1.push_back(std::make_pair(0.446357, 0.454291));
  dim1.push_back(std::make_pair(0.452574, 0.527445));
  dim1.push_back(std::make_pair(0.457876, 0.504913));
  dim1.push_back(std::make_pair(0.476169, 0.529129));
  dim1.push_back(std::make_pair(0.479762, 0.51536));
  dim1.push_back(std::make_pair(0.479762, 0.528202));
  dim1.push_back(std::make_pair(0.479762, 0.612093));
  dim1.push_back(std::make_pair(0.499352, 0.552841));
  dim1.push_back(std::make_pair(0.515218, 0.579041));
  dim1.push_back(std::make_pair(0.515218, 0.696697));
  dim1.push_back(std::make_pair(0.51536, 0.661572));
  dim1.push_back(std::make_pair(0.516364, 0.579478));
  dim1.push_back(std::make_pair(0.524837, 0.747096));
  dim1.push_back(std::make_pair(0.530507, 0.611131));
  dim1.push_back(std::make_pair(0.540148, 0.560137));
  dim1.push_back(std::make_pair(0.544518, 0.547481));
  dim1.push_back(std::make_pair(0.547481, 0.595119));
  dim1.push_back(std::make_pair(0.552251, 0.635619));
  dim1.push_back(std::make_pair(0.567289, 0.595748));
  dim1.push_back(std::make_pair(0.571876, 0.593348));
  dim1.push_back(std::make_pair(0.597383, 0.698937));
  dim1.push_back(std::make_pair(0.629459, 0.65752));
  dim1.push_back(std::make_pair(0.634643, 0.667667));
  dim1.push_back(std::make_pair(0.634892, 0.667901));


  std::vector< std::pair<double, double> > dim2;
  dim2.push_back(std::make_pair(0.159766, 0.514062));
  dim2.push_back(std::make_pair(0.397435, 0.857648));
  dim2.push_back(std::make_pair(0.424396, 0.494389));
  dim2.push_back(std::make_pair(0.426268, 0.831535));
  dim2.push_back(std::make_pair(0.443989, 0.737309));
  dim2.push_back(std::make_pair(0.521027, 0.55268));
  dim2.push_back(std::make_pair(0.5294, 0.827161));
  dim2.push_back(std::make_pair(0.531269, 0.975176));
  dim2.push_back(std::make_pair(0.538447, 0.893562));
  dim2.push_back(std::make_pair(0.559145, 0.78499));
  dim2.push_back(std::make_pair(0.560137, 0.764799));
  dim2.push_back(std::make_pair(0.562944, 0.672131));
  dim2.push_back(std::make_pair(0.568065, 0.833694));
  dim2.push_back(std::make_pair(0.582169, 0.614194));
  dim2.push_back(std::make_pair(0.595119, 0.99503));
  dim2.push_back(std::make_pair(0.6178, 0.645186));
  dim2.push_back(std::make_pair(0.633861, 0.894727));
  dim2.push_back(std::make_pair(0.65058, 0.720305));
  dim2.push_back(std::make_pair(0.654454, 0.882628));
  dim2.push_back(std::make_pair(0.656591, 0.911192));
  dim2.push_back(std::make_pair(0.663122, 0.835692));
  dim2.push_back(std::make_pair(0.677699, 0.828483));
  dim2.push_back(std::make_pair(0.684151, 0.714499));
  dim2.push_back(std::make_pair(0.691498, 0.79696));
  dim2.push_back(std::make_pair(0.692066, 0.836093));
  dim2.push_back(std::make_pair(0.696697, 0.924678));
  dim2.push_back(std::make_pair(0.698937, 0.888513));
  dim2.push_back(std::make_pair(0.708468, 0.825471));
  dim2.push_back(std::make_pair(0.714853, 0.998385));
  dim2.push_back(std::make_pair(0.723869, 0.972154));
  dim2.push_back(std::make_pair(0.725867, 0.900344));
  dim2.push_back(std::make_pair(0.726513, 0.972148));
  dim2.push_back(std::make_pair(0.728115, 0.853808));
  dim2.push_back(std::make_pair(0.730914, 0.95898));
  dim2.push_back(std::make_pair(0.734646, 0.940408));
  dim2.push_back(std::make_pair(0.735457, 0.784873));
  dim2.push_back(std::make_pair(0.737998, 0.989994));
  dim2.push_back(std::make_pair(0.749057, 0.796175));
  dim2.push_back(std::make_pair(0.755644, 0.760529));
  dim2.push_back(std::make_pair(0.763753, 0.770437));
  dim2.push_back(std::make_pair(0.775107, 0.801097));
  dim2.push_back(std::make_pair(0.776465, 0.808596));
  dim2.push_back(std::make_pair(0.78092, 0.836048));
  dim2.push_back(std::make_pair(0.784286, 0.957797));
  dim2.push_back(std::make_pair(0.792679, 0.834564));
  dim2.push_back(std::make_pair(0.792679, 0.940152));
  dim2.push_back(std::make_pair(0.794893, 0.919911));
  dim2.push_back(std::make_pair(0.808823, 0.986678));
  dim2.push_back(std::make_pair(0.813632, 0.914683));
  dim2.push_back(std::make_pair(0.823974, 0.96578));
  dim2.push_back(std::make_pair(0.826181, 0.952838));
  dim2.push_back(std::make_pair(0.839664, 0.932946));
  dim2.push_back(std::make_pair(0.843357, 0.921802));
  dim2.push_back(std::make_pair(0.84446, 0.950299));
  dim2.push_back(std::make_pair(0.84446, 0.966199));
  dim2.push_back(std::make_pair(0.848269, 0.937031));
  dim2.push_back(std::make_pair(0.856414, 0.911481));
  dim2.push_back(std::make_pair(0.882767, 0.919131));
  dim2.push_back(std::make_pair(0.882767, 0.939881));
  dim2.push_back(std::make_pair(0.889066, 0.972132));
  dim2.push_back(std::make_pair(0.891205, 0.917646));
  dim2.push_back(std::make_pair(0.891205, 0.920893));
  dim2.push_back(std::make_pair(0.895655, 0.901842));
  dim2.push_back(std::make_pair(0.911343, 0.930285));
  dim2.push_back(std::make_pair(0.916196, 0.933545));
  dim2.push_back(std::make_pair(0.93033, 0.950769));
  dim2.push_back(std::make_pair(0.936075, 0.936608));


  std::vector< std::vector<double> > betti_numbers;
  std::vector<double> bn_dim_0;
  bn_dim_0.push_back(0.000103006);
  std::vector<double> bn_dim_1;
  std::vector<double> bn_dim_2;
  std::vector<double> bn_dim_3;
  betti_numbers.push_back(bn_dim_0);
  betti_numbers.push_back(bn_dim_1);
  betti_numbers.push_back(bn_dim_2);
  betti_numbers.push_back(bn_dim_3);



  Bitmap_cubical_complex< Bitmap_cubical_complex_base<double> > b("random_cubical_complex");
  Compute_persistence_with_phat< Bitmap_cubical_complex< Bitmap_cubical_complex_base<double> > > phat(&b);
  phat::persistence_pairs pairs = phat.compute_persistence_pairs_dualized_chunk_reduction();
  std::pair< std::vector< std::vector<double> >, std::vector< std::vector< std::pair<double, double> > > > persistence = phat.get_the_intervals(pairs);



  //compare Betti numbers:  
  BOOST_CHECK(betti_numbers.size() == persistence.first.size());
  for (size_t dim = 0; dim != persistence.first.size() - 1; ++dim) {
    BOOST_CHECK(betti_numbers[dim].size() == persistence.first[dim].size());

    for (size_t i = 0; i != persistence.first[dim].size(); ++i) {
      BOOST_CHECK(fabs(persistence.first[dim][i] - betti_numbers[dim][i]) < EPSILON);
    }

  }



  //compare persistence :
  //first check if we get persistence in the same dimensions:
  BOOST_CHECK(persistence.second.size() == 3);

  //check dimension 0:
  BOOST_CHECK(persistence.second[0].size() == dim0.size());
  for (size_t i = 0; i != persistence.second[0].size(); ++i) {
    BOOST_CHECK((fabs(persistence.second[0][i].first - (double) dim0[i].first) < EPSILON));
    BOOST_CHECK((fabs(persistence.second[0][i].second - (double) dim0[i].second) < EPSILON));
  }

  //check dimension 1:
  BOOST_CHECK(persistence.second[1].size() == dim1.size());
  for (size_t i = 0; i != persistence.second[1].size(); ++i) {
    BOOST_CHECK((fabs(persistence.second[1][i].first - (double) dim1[i].first) < EPSILON));
    BOOST_CHECK((fabs(persistence.second[1][i].second - (double) dim1[i].second) < EPSILON));
  }

  //check dimension 2:
  BOOST_CHECK(persistence.second[2].size() == dim2.size());
  for (size_t i = 0; i != persistence.second[2].size(); ++i) {
    BOOST_CHECK((fabs(persistence.second[2][i].first - (double) dim2[i].first) < EPSILON));
    BOOST_CHECK((fabs(persistence.second[2][i].second - (double) dim2[i].second) < EPSILON));
  }
}

BOOST_AUTO_TEST_CASE(phat_cubical_standard_reductions) {
  double EPSILON = 0.0000005;

  std::vector< std::pair<double, double> > dim0;
  dim0.push_back(std::make_pair(0.000356789, 0.0209834));
  dim0.push_back(std::make_pair(0.00137761, 0.178343));
  dim0.push_back(std::make_pair(0.00296197, 0.182727));
  dim0.push_back(std::make_pair(0.00406146, 0.138371));
  dim0.push_back(std::make_pair(0.00572209, 0.0095895));
  dim0.push_back(std::make_pair(0.00724978, 0.088271));
  dim0.push_back(std::make_pair(0.00736106, 0.138713));
  dim0.push_back(std::make_pair(0.00779113, 0.176245));
  dim0.push_back(std::make_pair(0.00841857, 0.088271));
  dim0.push_back(std::make_pair(0.0101254, 0.0272611));
  dim0.push_back(std::make_pair(0.0133762, 0.11487));
  dim0.push_back(std::make_pair(0.0134719, 0.116554));
  dim0.push_back(std::make_pair(0.0145639, 0.155647));
  dim0.push_back(std::make_pair(0.0157889, 0.0683329));
  dim0.push_back(std::make_pair(0.0161422, 0.0476644));
  dim0.push_back(std::make_pair(0.0173804, 0.0743558));
  dim0.push_back(std::make_pair(0.0184607, 0.25154));
  dim0.push_back(std::make_pair(0.0186499, 0.339905));
  dim0.push_back(std::make_pair(0.0247987, 0.0559076));
  dim0.push_back(std::make_pair(0.0250133, 0.10707));
  dim0.push_back(std::make_pair(0.0251446, 0.222457));
  dim0.push_back(std::make_pair(0.026607, 0.12124));
  dim0.push_back(std::make_pair(0.0279297, 0.226968));
  dim0.push_back(std::make_pair(0.0294929, 0.0743558));
  dim0.push_back(std::make_pair(0.0321475, 0.0959911));
  dim0.push_back(std::make_pair(0.0338132, 0.15912));
  dim0.push_back(std::make_pair(0.0345952, 0.0956828));
  dim0.push_back(std::make_pair(0.0350849, 0.0476644));
  dim0.push_back(std::make_pair(0.0354295, 0.138721));
  dim0.push_back(std::make_pair(0.0358397, 0.138721));
  dim0.push_back(std::make_pair(0.0369627, 0.114244));
  dim0.push_back(std::make_pair(0.0382787, 0.153073));
  dim0.push_back(std::make_pair(0.0389954, 0.2347));
  dim0.push_back(std::make_pair(0.0424284, 0.0873569));
  dim0.push_back(std::make_pair(0.0460244, 0.0628887));
  dim0.push_back(std::make_pair(0.0460566, 0.100457));
  dim0.push_back(std::make_pair(0.0548351, 0.0956828));
  dim0.push_back(std::make_pair(0.0568281, 0.108871));
  dim0.push_back(std::make_pair(0.0598285, 0.0894398));
  dim0.push_back(std::make_pair(0.0643188, 0.180772));
  dim0.push_back(std::make_pair(0.067139, 0.217477));
  dim0.push_back(std::make_pair(0.0685163, 0.149041));
  dim0.push_back(std::make_pair(0.069365, 0.0840998));
  dim0.push_back(std::make_pair(0.0725224, 0.0873569));
  dim0.push_back(std::make_pair(0.0726488, 0.141462));
  dim0.push_back(std::make_pair(0.0750035, 0.156528));
  dim0.push_back(std::make_pair(0.0784113, 0.136386));
  dim0.push_back(std::make_pair(0.082815, 0.211298));
  dim0.push_back(std::make_pair(0.0857338, 0.186306));
  dim0.push_back(std::make_pair(0.0944081, 0.133989));
  dim0.push_back(std::make_pair(0.0958755, 0.238023));
  dim0.push_back(std::make_pair(0.114808, 0.124519));
  dim0.push_back(std::make_pair(0.142382, 0.20245));
  dim0.push_back(std::make_pair(0.147678, 0.161484));
  dim0.push_back(std::make_pair(0.158924, 0.217477));

  std::vector< std::pair<double, double> > dim1;
  dim1.push_back(std::make_pair(0.0896383, 0.307383));
  dim1.push_back(std::make_pair(0.0896383, 0.574247));
  dim1.push_back(std::make_pair(0.0922039, 0.337222));
  dim1.push_back(std::make_pair(0.0952979, 0.635619));
  dim1.push_back(std::make_pair(0.10707, 0.65752));
  dim1.push_back(std::make_pair(0.112554, 0.437864));
  dim1.push_back(std::make_pair(0.114244, 0.322239));
  dim1.push_back(std::make_pair(0.115086, 0.31284));
  dim1.push_back(std::make_pair(0.131682, 0.609287));
  dim1.push_back(std::make_pair(0.138669, 0.15837));
  dim1.push_back(std::make_pair(0.138669, 0.212147));
  dim1.push_back(std::make_pair(0.138669, 0.299792));
  dim1.push_back(std::make_pair(0.146283, 0.375957));
  dim1.push_back(std::make_pair(0.150632, 0.354247));
  dim1.push_back(std::make_pair(0.150632, 0.597383));
  dim1.push_back(std::make_pair(0.152943, 0.283281));
  dim1.push_back(std::make_pair(0.15912, 0.642942));
  dim1.push_back(std::make_pair(0.160004, 0.683229));
  dim1.push_back(std::make_pair(0.165983, 0.287496));
  dim1.push_back(std::make_pair(0.165983, 0.307383));
  dim1.push_back(std::make_pair(0.171803, 0.190303));
  dim1.push_back(std::make_pair(0.177224, 0.559145));
  dim1.push_back(std::make_pair(0.177959, 0.219433));
  dim1.push_back(std::make_pair(0.177959, 0.481775));
  dim1.push_back(std::make_pair(0.179939, 0.426291));
  dim1.push_back(std::make_pair(0.180772, 0.600957));
  dim1.push_back(std::make_pair(0.182727, 0.593896));
  dim1.push_back(std::make_pair(0.184009, 0.504145));
  dim1.push_back(std::make_pair(0.193217, 0.714115));
  dim1.push_back(std::make_pair(0.20074, 0.26084));
  dim1.push_back(std::make_pair(0.203173, 0.493206));
  dim1.push_back(std::make_pair(0.203173, 0.538447));
  dim1.push_back(std::make_pair(0.203361, 0.618832));
  dim1.push_back(std::make_pair(0.203361, 0.656591));
  dim1.push_back(std::make_pair(0.205348, 0.427394));
  dim1.push_back(std::make_pair(0.205348, 0.556077));
  dim1.push_back(std::make_pair(0.213038, 0.642942));
  dim1.push_back(std::make_pair(0.21696, 0.397078));
  dim1.push_back(std::make_pair(0.21696, 0.493206));
  dim1.push_back(std::make_pair(0.218123, 0.812171));
  dim1.push_back(std::make_pair(0.222869, 0.249479));
  dim1.push_back(std::make_pair(0.2241, 0.692066));
  dim1.push_back(std::make_pair(0.224358, 0.667901));
  dim1.push_back(std::make_pair(0.224437, 0.254382));
  dim1.push_back(std::make_pair(0.224437, 0.393736));
  dim1.push_back(std::make_pair(0.224437, 0.562944));
  dim1.push_back(std::make_pair(0.226809, 0.342489));
  dim1.push_back(std::make_pair(0.227838, 0.58166));
  dim1.push_back(std::make_pair(0.227838, 0.671021));
  dim1.push_back(std::make_pair(0.229861, 0.259465));
  dim1.push_back(std::make_pair(0.233203, 0.504145));
  dim1.push_back(std::make_pair(0.2347, 0.253456));
  dim1.push_back(std::make_pair(0.2347, 0.401604));
  dim1.push_back(std::make_pair(0.23542, 0.509011));
  dim1.push_back(std::make_pair(0.23542, 0.656926));
  dim1.push_back(std::make_pair(0.238023, 0.306904));
  dim1.push_back(std::make_pair(0.238187, 0.560137));
  dim1.push_back(std::make_pair(0.242099, 0.375768));
  dim1.push_back(std::make_pair(0.251213, 0.419192));
  dim1.push_back(std::make_pair(0.251213, 0.470896));
  dim1.push_back(std::make_pair(0.25154, 0.441861));
  dim1.push_back(std::make_pair(0.252778, 0.683845));
  dim1.push_back(std::make_pair(0.253932, 0.303766));
  dim1.push_back(std::make_pair(0.25521, 0.318305));
  dim1.push_back(std::make_pair(0.25521, 0.366637));
  dim1.push_back(std::make_pair(0.256709, 0.552841));
  dim1.push_back(std::make_pair(0.256709, 0.566244));
  dim1.push_back(std::make_pair(0.256923, 0.662304));
  dim1.push_back(std::make_pair(0.261521, 0.634892));
  dim1.push_back(std::make_pair(0.270009, 0.642002));
  dim1.push_back(std::make_pair(0.271325, 0.329129));
  dim1.push_back(std::make_pair(0.27284, 0.801019));
  dim1.push_back(std::make_pair(0.273637, 0.283281));
  dim1.push_back(std::make_pair(0.273637, 0.505287));
  dim1.push_back(std::make_pair(0.276398, 0.667667));
  dim1.push_back(std::make_pair(0.279702, 0.365257));
  dim1.push_back(std::make_pair(0.280491, 0.579041));
  dim1.push_back(std::make_pair(0.284502, 0.405565));
  dim1.push_back(std::make_pair(0.284932, 0.377696));
  dim1.push_back(std::make_pair(0.291091, 0.624392));
  dim1.push_back(std::make_pair(0.291315, 0.443989));
  dim1.push_back(std::make_pair(0.293683, 0.522518));
  dim1.push_back(std::make_pair(0.293683, 0.741924));
  dim1.push_back(std::make_pair(0.307746, 0.588104));
  dim1.push_back(std::make_pair(0.318734, 0.498793));
  dim1.push_back(std::make_pair(0.32605, 0.52335));
  dim1.push_back(std::make_pair(0.335347, 0.58166));
  dim1.push_back(std::make_pair(0.336204, 0.71353));
  dim1.push_back(std::make_pair(0.342399, 0.595119));
  dim1.push_back(std::make_pair(0.348409, 0.664241));
  dim1.push_back(std::make_pair(0.357035, 0.70087));
  dim1.push_back(std::make_pair(0.357874, 0.476921));
  dim1.push_back(std::make_pair(0.36049, 0.640687));
  dim1.push_back(std::make_pair(0.363861, 0.494538));
  dim1.push_back(std::make_pair(0.363861, 0.494538));
  dim1.push_back(std::make_pair(0.364503, 0.640687));
  dim1.push_back(std::make_pair(0.364503, 0.640687));
  dim1.push_back(std::make_pair(0.364571, 0.661572));
  dim1.push_back(std::make_pair(0.367343, 0.597987));
  dim1.push_back(std::make_pair(0.387655, 0.569132));
  dim1.push_back(std::make_pair(0.397369, 0.644901));
  dim1.push_back(std::make_pair(0.397369, 0.659947));
  dim1.push_back(std::make_pair(0.397386, 0.431372));
  dim1.push_back(std::make_pair(0.39905, 0.440457));
  dim1.push_back(std::make_pair(0.399317, 0.566299));
  dim1.push_back(std::make_pair(0.4005, 0.406971));
  dim1.push_back(std::make_pair(0.400725, 0.484885));
  dim1.push_back(std::make_pair(0.405866, 0.426268));
  dim1.push_back(std::make_pair(0.407038, 0.478693));
  dim1.push_back(std::make_pair(0.409665, 0.429602));
  dim1.push_back(std::make_pair(0.418432, 0.476169));
  dim1.push_back(std::make_pair(0.420472, 0.478693));
  dim1.push_back(std::make_pair(0.424782, 0.626694));
  dim1.push_back(std::make_pair(0.424782, 0.629479));
  dim1.push_back(std::make_pair(0.438875, 0.504145));
  dim1.push_back(std::make_pair(0.439381, 0.539878));
  dim1.push_back(std::make_pair(0.439816, 0.527948));
  dim1.push_back(std::make_pair(0.442226, 0.529129));
  dim1.push_back(std::make_pair(0.446357, 0.454291));
  dim1.push_back(std::make_pair(0.452574, 0.527445));
  dim1.push_back(std::make_pair(0.457876, 0.504913));
  dim1.push_back(std::make_pair(0.476169, 0.529129));
  dim1.push_back(std::make_pair(0.479762, 0.51536));
  dim1.push_back(std::make_pair(0.479762, 0.528202));
  dim1.push_back(std::make_pair(0.479762, 0.612093));
  dim1.push_back(std::make_pair(0.499352, 0.552841));
  dim1.push_back(std::make_pair(0.515218, 0.579041));
  dim1.push_back(std::make_pair(0.515218, 0.696697));
  dim1.push_back(std::make_pair(0.51536, 0.661572));
  dim1.push_back(std::make_pair(0.516364, 0.579478));
  dim1.push_back(std::make_pair(0.524837, 0.747096));
  dim1.push_back(std::make_pair(0.530507, 0.611131));
  dim1.push_back(std::make_pair(0.540148, 0.560137));
  dim1.push_back(std::make_pair(0.544518, 0.547481));
  dim1.push_back(std::make_pair(0.547481, 0.595119));
  dim1.push_back(std::make_pair(0.552251, 0.635619));
  dim1.push_back(std::make_pair(0.567289, 0.595748));
  dim1.push_back(std::make_pair(0.571876, 0.593348));
  dim1.push_back(std::make_pair(0.597383, 0.698937));
  dim1.push_back(std::make_pair(0.629459, 0.65752));
  dim1.push_back(std::make_pair(0.634643, 0.667667));
  dim1.push_back(std::make_pair(0.634892, 0.667901));


  std::vector< std::pair<double, double> > dim2;
  dim2.push_back(std::make_pair(0.159766, 0.514062));
  dim2.push_back(std::make_pair(0.397435, 0.857648));
  dim2.push_back(std::make_pair(0.424396, 0.494389));
  dim2.push_back(std::make_pair(0.426268, 0.831535));
  dim2.push_back(std::make_pair(0.443989, 0.737309));
  dim2.push_back(std::make_pair(0.521027, 0.55268));
  dim2.push_back(std::make_pair(0.5294, 0.827161));
  dim2.push_back(std::make_pair(0.531269, 0.975176));
  dim2.push_back(std::make_pair(0.538447, 0.893562));
  dim2.push_back(std::make_pair(0.559145, 0.78499));
  dim2.push_back(std::make_pair(0.560137, 0.764799));
  dim2.push_back(std::make_pair(0.562944, 0.672131));
  dim2.push_back(std::make_pair(0.568065, 0.833694));
  dim2.push_back(std::make_pair(0.582169, 0.614194));
  dim2.push_back(std::make_pair(0.595119, 0.99503));
  dim2.push_back(std::make_pair(0.6178, 0.645186));
  dim2.push_back(std::make_pair(0.633861, 0.894727));
  dim2.push_back(std::make_pair(0.65058, 0.720305));
  dim2.push_back(std::make_pair(0.654454, 0.882628));
  dim2.push_back(std::make_pair(0.656591, 0.911192));
  dim2.push_back(std::make_pair(0.663122, 0.835692));
  dim2.push_back(std::make_pair(0.677699, 0.828483));
  dim2.push_back(std::make_pair(0.684151, 0.714499));
  dim2.push_back(std::make_pair(0.691498, 0.79696));
  dim2.push_back(std::make_pair(0.692066, 0.836093));
  dim2.push_back(std::make_pair(0.696697, 0.924678));
  dim2.push_back(std::make_pair(0.698937, 0.888513));
  dim2.push_back(std::make_pair(0.708468, 0.825471));
  dim2.push_back(std::make_pair(0.714853, 0.998385));
  dim2.push_back(std::make_pair(0.723869, 0.972154));
  dim2.push_back(std::make_pair(0.725867, 0.900344));
  dim2.push_back(std::make_pair(0.726513, 0.972148));
  dim2.push_back(std::make_pair(0.728115, 0.853808));
  dim2.push_back(std::make_pair(0.730914, 0.95898));
  dim2.push_back(std::make_pair(0.734646, 0.940408));
  dim2.push_back(std::make_pair(0.735457, 0.784873));
  dim2.push_back(std::make_pair(0.737998, 0.989994));
  dim2.push_back(std::make_pair(0.749057, 0.796175));
  dim2.push_back(std::make_pair(0.755644, 0.760529));
  dim2.push_back(std::make_pair(0.763753, 0.770437));
  dim2.push_back(std::make_pair(0.775107, 0.801097));
  dim2.push_back(std::make_pair(0.776465, 0.808596));
  dim2.push_back(std::make_pair(0.78092, 0.836048));
  dim2.push_back(std::make_pair(0.784286, 0.957797));
  dim2.push_back(std::make_pair(0.792679, 0.834564));
  dim2.push_back(std::make_pair(0.792679, 0.940152));
  dim2.push_back(std::make_pair(0.794893, 0.919911));
  dim2.push_back(std::make_pair(0.808823, 0.986678));
  dim2.push_back(std::make_pair(0.813632, 0.914683));
  dim2.push_back(std::make_pair(0.823974, 0.96578));
  dim2.push_back(std::make_pair(0.826181, 0.952838));
  dim2.push_back(std::make_pair(0.839664, 0.932946));
  dim2.push_back(std::make_pair(0.843357, 0.921802));
  dim2.push_back(std::make_pair(0.84446, 0.950299));
  dim2.push_back(std::make_pair(0.84446, 0.966199));
  dim2.push_back(std::make_pair(0.848269, 0.937031));
  dim2.push_back(std::make_pair(0.856414, 0.911481));
  dim2.push_back(std::make_pair(0.882767, 0.919131));
  dim2.push_back(std::make_pair(0.882767, 0.939881));
  dim2.push_back(std::make_pair(0.889066, 0.972132));
  dim2.push_back(std::make_pair(0.891205, 0.917646));
  dim2.push_back(std::make_pair(0.891205, 0.920893));
  dim2.push_back(std::make_pair(0.895655, 0.901842));
  dim2.push_back(std::make_pair(0.911343, 0.930285));
  dim2.push_back(std::make_pair(0.916196, 0.933545));
  dim2.push_back(std::make_pair(0.93033, 0.950769));
  dim2.push_back(std::make_pair(0.936075, 0.936608));


  std::vector< std::vector<double> > betti_numbers;
  std::vector<double> bn_dim_0;
  bn_dim_0.push_back(0.000103006);
  std::vector<double> bn_dim_1;
  std::vector<double> bn_dim_2;
  std::vector<double> bn_dim_3;
  betti_numbers.push_back(bn_dim_0);
  betti_numbers.push_back(bn_dim_1);
  betti_numbers.push_back(bn_dim_2);
  betti_numbers.push_back(bn_dim_3);



  Bitmap_cubical_complex< Bitmap_cubical_complex_base<double> > b("random_cubical_complex");
  Compute_persistence_with_phat< Bitmap_cubical_complex< Bitmap_cubical_complex_base<double> > > phat(&b);
  phat::persistence_pairs pairs = phat.compute_persistence_pairs_standard_reduction();
  std::pair< std::vector< std::vector<double> >, std::vector< std::vector< std::pair<double, double> > > > persistence = phat.get_the_intervals(pairs);



  //compare Betti numbers:  
  BOOST_CHECK(betti_numbers.size() == persistence.first.size());
  for (size_t dim = 0; dim != persistence.first.size() - 1; ++dim) {
    BOOST_CHECK(betti_numbers[dim].size() == persistence.first[dim].size());

    for (size_t i = 0; i != persistence.first[dim].size(); ++i) {
      BOOST_CHECK(fabs(persistence.first[dim][i] - betti_numbers[dim][i]) < EPSILON);
    }

  }



  //compare persistence :
  //first check if we get persistence in the same dimensions:
  BOOST_CHECK(persistence.second.size() == 3);

  //check dimension 0:
  BOOST_CHECK(persistence.second[0].size() == dim0.size());
  for (size_t i = 0; i != persistence.second[0].size(); ++i) {
    BOOST_CHECK((fabs(persistence.second[0][i].first - (double) dim0[i].first) < EPSILON));
    BOOST_CHECK((fabs(persistence.second[0][i].second - (double) dim0[i].second) < EPSILON));
  }

  //check dimension 1:
  BOOST_CHECK(persistence.second[1].size() == dim1.size());
  for (size_t i = 0; i != persistence.second[1].size(); ++i) {
    BOOST_CHECK((fabs(persistence.second[1][i].first - (double) dim1[i].first) < EPSILON));
    BOOST_CHECK((fabs(persistence.second[1][i].second - (double) dim1[i].second) < EPSILON));
  }

  //check dimension 2:
  BOOST_CHECK(persistence.second[2].size() == dim2.size());
  for (size_t i = 0; i != persistence.second[2].size(); ++i) {
    BOOST_CHECK((fabs(persistence.second[2][i].first - (double) dim2[i].first) < EPSILON));
    BOOST_CHECK((fabs(persistence.second[2][i].second - (double) dim2[i].second) < EPSILON));
  }
}

BOOST_AUTO_TEST_CASE(phat_cubical_spectral_sequence) {
  double EPSILON = 0.0000005;

  std::vector< std::pair<double, double> > dim0;
  dim0.push_back(std::make_pair(0.000356789, 0.0209834));
  dim0.push_back(std::make_pair(0.00137761, 0.178343));
  dim0.push_back(std::make_pair(0.00296197, 0.182727));
  dim0.push_back(std::make_pair(0.00406146, 0.138371));
  dim0.push_back(std::make_pair(0.00572209, 0.0095895));
  dim0.push_back(std::make_pair(0.00724978, 0.088271));
  dim0.push_back(std::make_pair(0.00736106, 0.138713));
  dim0.push_back(std::make_pair(0.00779113, 0.176245));
  dim0.push_back(std::make_pair(0.00841857, 0.088271));
  dim0.push_back(std::make_pair(0.0101254, 0.0272611));
  dim0.push_back(std::make_pair(0.0133762, 0.11487));
  dim0.push_back(std::make_pair(0.0134719, 0.116554));
  dim0.push_back(std::make_pair(0.0145639, 0.155647));
  dim0.push_back(std::make_pair(0.0157889, 0.0683329));
  dim0.push_back(std::make_pair(0.0161422, 0.0476644));
  dim0.push_back(std::make_pair(0.0173804, 0.0743558));
  dim0.push_back(std::make_pair(0.0184607, 0.25154));
  dim0.push_back(std::make_pair(0.0186499, 0.339905));
  dim0.push_back(std::make_pair(0.0247987, 0.0559076));
  dim0.push_back(std::make_pair(0.0250133, 0.10707));
  dim0.push_back(std::make_pair(0.0251446, 0.222457));
  dim0.push_back(std::make_pair(0.026607, 0.12124));
  dim0.push_back(std::make_pair(0.0279297, 0.226968));
  dim0.push_back(std::make_pair(0.0294929, 0.0743558));
  dim0.push_back(std::make_pair(0.0321475, 0.0959911));
  dim0.push_back(std::make_pair(0.0338132, 0.15912));
  dim0.push_back(std::make_pair(0.0345952, 0.0956828));
  dim0.push_back(std::make_pair(0.0350849, 0.0476644));
  dim0.push_back(std::make_pair(0.0354295, 0.138721));
  dim0.push_back(std::make_pair(0.0358397, 0.138721));
  dim0.push_back(std::make_pair(0.0369627, 0.114244));
  dim0.push_back(std::make_pair(0.0382787, 0.153073));
  dim0.push_back(std::make_pair(0.0389954, 0.2347));
  dim0.push_back(std::make_pair(0.0424284, 0.0873569));
  dim0.push_back(std::make_pair(0.0460244, 0.0628887));
  dim0.push_back(std::make_pair(0.0460566, 0.100457));
  dim0.push_back(std::make_pair(0.0548351, 0.0956828));
  dim0.push_back(std::make_pair(0.0568281, 0.108871));
  dim0.push_back(std::make_pair(0.0598285, 0.0894398));
  dim0.push_back(std::make_pair(0.0643188, 0.180772));
  dim0.push_back(std::make_pair(0.067139, 0.217477));
  dim0.push_back(std::make_pair(0.0685163, 0.149041));
  dim0.push_back(std::make_pair(0.069365, 0.0840998));
  dim0.push_back(std::make_pair(0.0725224, 0.0873569));
  dim0.push_back(std::make_pair(0.0726488, 0.141462));
  dim0.push_back(std::make_pair(0.0750035, 0.156528));
  dim0.push_back(std::make_pair(0.0784113, 0.136386));
  dim0.push_back(std::make_pair(0.082815, 0.211298));
  dim0.push_back(std::make_pair(0.0857338, 0.186306));
  dim0.push_back(std::make_pair(0.0944081, 0.133989));
  dim0.push_back(std::make_pair(0.0958755, 0.238023));
  dim0.push_back(std::make_pair(0.114808, 0.124519));
  dim0.push_back(std::make_pair(0.142382, 0.20245));
  dim0.push_back(std::make_pair(0.147678, 0.161484));
  dim0.push_back(std::make_pair(0.158924, 0.217477));

  std::vector< std::pair<double, double> > dim1;
  dim1.push_back(std::make_pair(0.0896383, 0.307383));
  dim1.push_back(std::make_pair(0.0896383, 0.574247));
  dim1.push_back(std::make_pair(0.0922039, 0.337222));
  dim1.push_back(std::make_pair(0.0952979, 0.635619));
  dim1.push_back(std::make_pair(0.10707, 0.65752));
  dim1.push_back(std::make_pair(0.112554, 0.437864));
  dim1.push_back(std::make_pair(0.114244, 0.322239));
  dim1.push_back(std::make_pair(0.115086, 0.31284));
  dim1.push_back(std::make_pair(0.131682, 0.609287));
  dim1.push_back(std::make_pair(0.138669, 0.15837));
  dim1.push_back(std::make_pair(0.138669, 0.212147));
  dim1.push_back(std::make_pair(0.138669, 0.299792));
  dim1.push_back(std::make_pair(0.146283, 0.375957));
  dim1.push_back(std::make_pair(0.150632, 0.354247));
  dim1.push_back(std::make_pair(0.150632, 0.597383));
  dim1.push_back(std::make_pair(0.152943, 0.283281));
  dim1.push_back(std::make_pair(0.15912, 0.642942));
  dim1.push_back(std::make_pair(0.160004, 0.683229));
  dim1.push_back(std::make_pair(0.165983, 0.287496));
  dim1.push_back(std::make_pair(0.165983, 0.307383));
  dim1.push_back(std::make_pair(0.171803, 0.190303));
  dim1.push_back(std::make_pair(0.177224, 0.559145));
  dim1.push_back(std::make_pair(0.177959, 0.219433));
  dim1.push_back(std::make_pair(0.177959, 0.481775));
  dim1.push_back(std::make_pair(0.179939, 0.426291));
  dim1.push_back(std::make_pair(0.180772, 0.600957));
  dim1.push_back(std::make_pair(0.182727, 0.593896));
  dim1.push_back(std::make_pair(0.184009, 0.504145));
  dim1.push_back(std::make_pair(0.193217, 0.714115));
  dim1.push_back(std::make_pair(0.20074, 0.26084));
  dim1.push_back(std::make_pair(0.203173, 0.493206));
  dim1.push_back(std::make_pair(0.203173, 0.538447));
  dim1.push_back(std::make_pair(0.203361, 0.618832));
  dim1.push_back(std::make_pair(0.203361, 0.656591));
  dim1.push_back(std::make_pair(0.205348, 0.427394));
  dim1.push_back(std::make_pair(0.205348, 0.556077));
  dim1.push_back(std::make_pair(0.213038, 0.642942));
  dim1.push_back(std::make_pair(0.21696, 0.397078));
  dim1.push_back(std::make_pair(0.21696, 0.493206));
  dim1.push_back(std::make_pair(0.218123, 0.812171));
  dim1.push_back(std::make_pair(0.222869, 0.249479));
  dim1.push_back(std::make_pair(0.2241, 0.692066));
  dim1.push_back(std::make_pair(0.224358, 0.667901));
  dim1.push_back(std::make_pair(0.224437, 0.254382));
  dim1.push_back(std::make_pair(0.224437, 0.393736));
  dim1.push_back(std::make_pair(0.224437, 0.562944));
  dim1.push_back(std::make_pair(0.226809, 0.342489));
  dim1.push_back(std::make_pair(0.227838, 0.58166));
  dim1.push_back(std::make_pair(0.227838, 0.671021));
  dim1.push_back(std::make_pair(0.229861, 0.259465));
  dim1.push_back(std::make_pair(0.233203, 0.504145));
  dim1.push_back(std::make_pair(0.2347, 0.253456));
  dim1.push_back(std::make_pair(0.2347, 0.401604));
  dim1.push_back(std::make_pair(0.23542, 0.509011));
  dim1.push_back(std::make_pair(0.23542, 0.656926));
  dim1.push_back(std::make_pair(0.238023, 0.306904));
  dim1.push_back(std::make_pair(0.238187, 0.560137));
  dim1.push_back(std::make_pair(0.242099, 0.375768));
  dim1.push_back(std::make_pair(0.251213, 0.419192));
  dim1.push_back(std::make_pair(0.251213, 0.470896));
  dim1.push_back(std::make_pair(0.25154, 0.441861));
  dim1.push_back(std::make_pair(0.252778, 0.683845));
  dim1.push_back(std::make_pair(0.253932, 0.303766));
  dim1.push_back(std::make_pair(0.25521, 0.318305));
  dim1.push_back(std::make_pair(0.25521, 0.366637));
  dim1.push_back(std::make_pair(0.256709, 0.552841));
  dim1.push_back(std::make_pair(0.256709, 0.566244));
  dim1.push_back(std::make_pair(0.256923, 0.662304));
  dim1.push_back(std::make_pair(0.261521, 0.634892));
  dim1.push_back(std::make_pair(0.270009, 0.642002));
  dim1.push_back(std::make_pair(0.271325, 0.329129));
  dim1.push_back(std::make_pair(0.27284, 0.801019));
  dim1.push_back(std::make_pair(0.273637, 0.283281));
  dim1.push_back(std::make_pair(0.273637, 0.505287));
  dim1.push_back(std::make_pair(0.276398, 0.667667));
  dim1.push_back(std::make_pair(0.279702, 0.365257));
  dim1.push_back(std::make_pair(0.280491, 0.579041));
  dim1.push_back(std::make_pair(0.284502, 0.405565));
  dim1.push_back(std::make_pair(0.284932, 0.377696));
  dim1.push_back(std::make_pair(0.291091, 0.624392));
  dim1.push_back(std::make_pair(0.291315, 0.443989));
  dim1.push_back(std::make_pair(0.293683, 0.522518));
  dim1.push_back(std::make_pair(0.293683, 0.741924));
  dim1.push_back(std::make_pair(0.307746, 0.588104));
  dim1.push_back(std::make_pair(0.318734, 0.498793));
  dim1.push_back(std::make_pair(0.32605, 0.52335));
  dim1.push_back(std::make_pair(0.335347, 0.58166));
  dim1.push_back(std::make_pair(0.336204, 0.71353));
  dim1.push_back(std::make_pair(0.342399, 0.595119));
  dim1.push_back(std::make_pair(0.348409, 0.664241));
  dim1.push_back(std::make_pair(0.357035, 0.70087));
  dim1.push_back(std::make_pair(0.357874, 0.476921));
  dim1.push_back(std::make_pair(0.36049, 0.640687));
  dim1.push_back(std::make_pair(0.363861, 0.494538));
  dim1.push_back(std::make_pair(0.363861, 0.494538));
  dim1.push_back(std::make_pair(0.364503, 0.640687));
  dim1.push_back(std::make_pair(0.364503, 0.640687));
  dim1.push_back(std::make_pair(0.364571, 0.661572));
  dim1.push_back(std::make_pair(0.367343, 0.597987));
  dim1.push_back(std::make_pair(0.387655, 0.569132));
  dim1.push_back(std::make_pair(0.397369, 0.644901));
  dim1.push_back(std::make_pair(0.397369, 0.659947));
  dim1.push_back(std::make_pair(0.397386, 0.431372));
  dim1.push_back(std::make_pair(0.39905, 0.440457));
  dim1.push_back(std::make_pair(0.399317, 0.566299));
  dim1.push_back(std::make_pair(0.4005, 0.406971));
  dim1.push_back(std::make_pair(0.400725, 0.484885));
  dim1.push_back(std::make_pair(0.405866, 0.426268));
  dim1.push_back(std::make_pair(0.407038, 0.478693));
  dim1.push_back(std::make_pair(0.409665, 0.429602));
  dim1.push_back(std::make_pair(0.418432, 0.476169));
  dim1.push_back(std::make_pair(0.420472, 0.478693));
  dim1.push_back(std::make_pair(0.424782, 0.626694));
  dim1.push_back(std::make_pair(0.424782, 0.629479));
  dim1.push_back(std::make_pair(0.438875, 0.504145));
  dim1.push_back(std::make_pair(0.439381, 0.539878));
  dim1.push_back(std::make_pair(0.439816, 0.527948));
  dim1.push_back(std::make_pair(0.442226, 0.529129));
  dim1.push_back(std::make_pair(0.446357, 0.454291));
  dim1.push_back(std::make_pair(0.452574, 0.527445));
  dim1.push_back(std::make_pair(0.457876, 0.504913));
  dim1.push_back(std::make_pair(0.476169, 0.529129));
  dim1.push_back(std::make_pair(0.479762, 0.51536));
  dim1.push_back(std::make_pair(0.479762, 0.528202));
  dim1.push_back(std::make_pair(0.479762, 0.612093));
  dim1.push_back(std::make_pair(0.499352, 0.552841));
  dim1.push_back(std::make_pair(0.515218, 0.579041));
  dim1.push_back(std::make_pair(0.515218, 0.696697));
  dim1.push_back(std::make_pair(0.51536, 0.661572));
  dim1.push_back(std::make_pair(0.516364, 0.579478));
  dim1.push_back(std::make_pair(0.524837, 0.747096));
  dim1.push_back(std::make_pair(0.530507, 0.611131));
  dim1.push_back(std::make_pair(0.540148, 0.560137));
  dim1.push_back(std::make_pair(0.544518, 0.547481));
  dim1.push_back(std::make_pair(0.547481, 0.595119));
  dim1.push_back(std::make_pair(0.552251, 0.635619));
  dim1.push_back(std::make_pair(0.567289, 0.595748));
  dim1.push_back(std::make_pair(0.571876, 0.593348));
  dim1.push_back(std::make_pair(0.597383, 0.698937));
  dim1.push_back(std::make_pair(0.629459, 0.65752));
  dim1.push_back(std::make_pair(0.634643, 0.667667));
  dim1.push_back(std::make_pair(0.634892, 0.667901));


  std::vector< std::pair<double, double> > dim2;
  dim2.push_back(std::make_pair(0.159766, 0.514062));
  dim2.push_back(std::make_pair(0.397435, 0.857648));
  dim2.push_back(std::make_pair(0.424396, 0.494389));
  dim2.push_back(std::make_pair(0.426268, 0.831535));
  dim2.push_back(std::make_pair(0.443989, 0.737309));
  dim2.push_back(std::make_pair(0.521027, 0.55268));
  dim2.push_back(std::make_pair(0.5294, 0.827161));
  dim2.push_back(std::make_pair(0.531269, 0.975176));
  dim2.push_back(std::make_pair(0.538447, 0.893562));
  dim2.push_back(std::make_pair(0.559145, 0.78499));
  dim2.push_back(std::make_pair(0.560137, 0.764799));
  dim2.push_back(std::make_pair(0.562944, 0.672131));
  dim2.push_back(std::make_pair(0.568065, 0.833694));
  dim2.push_back(std::make_pair(0.582169, 0.614194));
  dim2.push_back(std::make_pair(0.595119, 0.99503));
  dim2.push_back(std::make_pair(0.6178, 0.645186));
  dim2.push_back(std::make_pair(0.633861, 0.894727));
  dim2.push_back(std::make_pair(0.65058, 0.720305));
  dim2.push_back(std::make_pair(0.654454, 0.882628));
  dim2.push_back(std::make_pair(0.656591, 0.911192));
  dim2.push_back(std::make_pair(0.663122, 0.835692));
  dim2.push_back(std::make_pair(0.677699, 0.828483));
  dim2.push_back(std::make_pair(0.684151, 0.714499));
  dim2.push_back(std::make_pair(0.691498, 0.79696));
  dim2.push_back(std::make_pair(0.692066, 0.836093));
  dim2.push_back(std::make_pair(0.696697, 0.924678));
  dim2.push_back(std::make_pair(0.698937, 0.888513));
  dim2.push_back(std::make_pair(0.708468, 0.825471));
  dim2.push_back(std::make_pair(0.714853, 0.998385));
  dim2.push_back(std::make_pair(0.723869, 0.972154));
  dim2.push_back(std::make_pair(0.725867, 0.900344));
  dim2.push_back(std::make_pair(0.726513, 0.972148));
  dim2.push_back(std::make_pair(0.728115, 0.853808));
  dim2.push_back(std::make_pair(0.730914, 0.95898));
  dim2.push_back(std::make_pair(0.734646, 0.940408));
  dim2.push_back(std::make_pair(0.735457, 0.784873));
  dim2.push_back(std::make_pair(0.737998, 0.989994));
  dim2.push_back(std::make_pair(0.749057, 0.796175));
  dim2.push_back(std::make_pair(0.755644, 0.760529));
  dim2.push_back(std::make_pair(0.763753, 0.770437));
  dim2.push_back(std::make_pair(0.775107, 0.801097));
  dim2.push_back(std::make_pair(0.776465, 0.808596));
  dim2.push_back(std::make_pair(0.78092, 0.836048));
  dim2.push_back(std::make_pair(0.784286, 0.957797));
  dim2.push_back(std::make_pair(0.792679, 0.834564));
  dim2.push_back(std::make_pair(0.792679, 0.940152));
  dim2.push_back(std::make_pair(0.794893, 0.919911));
  dim2.push_back(std::make_pair(0.808823, 0.986678));
  dim2.push_back(std::make_pair(0.813632, 0.914683));
  dim2.push_back(std::make_pair(0.823974, 0.96578));
  dim2.push_back(std::make_pair(0.826181, 0.952838));
  dim2.push_back(std::make_pair(0.839664, 0.932946));
  dim2.push_back(std::make_pair(0.843357, 0.921802));
  dim2.push_back(std::make_pair(0.84446, 0.950299));
  dim2.push_back(std::make_pair(0.84446, 0.966199));
  dim2.push_back(std::make_pair(0.848269, 0.937031));
  dim2.push_back(std::make_pair(0.856414, 0.911481));
  dim2.push_back(std::make_pair(0.882767, 0.919131));
  dim2.push_back(std::make_pair(0.882767, 0.939881));
  dim2.push_back(std::make_pair(0.889066, 0.972132));
  dim2.push_back(std::make_pair(0.891205, 0.917646));
  dim2.push_back(std::make_pair(0.891205, 0.920893));
  dim2.push_back(std::make_pair(0.895655, 0.901842));
  dim2.push_back(std::make_pair(0.911343, 0.930285));
  dim2.push_back(std::make_pair(0.916196, 0.933545));
  dim2.push_back(std::make_pair(0.93033, 0.950769));
  dim2.push_back(std::make_pair(0.936075, 0.936608));


  std::vector< std::vector<double> > betti_numbers;
  std::vector<double> bn_dim_0;
  bn_dim_0.push_back(0.000103006);
  std::vector<double> bn_dim_1;
  std::vector<double> bn_dim_2;
  std::vector<double> bn_dim_3;
  betti_numbers.push_back(bn_dim_0);
  betti_numbers.push_back(bn_dim_1);
  betti_numbers.push_back(bn_dim_2);
  betti_numbers.push_back(bn_dim_3);



  Bitmap_cubical_complex< Bitmap_cubical_complex_base<double> > b("random_cubical_complex");
  Compute_persistence_with_phat< Bitmap_cubical_complex< Bitmap_cubical_complex_base<double> > > phat(&b);
  phat::persistence_pairs pairs = phat.compute_persistence_pairs_spectral_sequence_reduction();
  std::pair< std::vector< std::vector<double> >, std::vector< std::vector< std::pair<double, double> > > > persistence = phat.get_the_intervals(pairs);



  //compare Betti numbers:  
  BOOST_CHECK(betti_numbers.size() == persistence.first.size());
  for (size_t dim = 0; dim != persistence.first.size() - 1; ++dim) {
    BOOST_CHECK(betti_numbers[dim].size() == persistence.first[dim].size());

    for (size_t i = 0; i != persistence.first[dim].size(); ++i) {
      BOOST_CHECK(fabs(persistence.first[dim][i] - betti_numbers[dim][i]) < EPSILON);
    }

  }



  //compare persistence :
  //first check if we get persistence in the same dimensions:
  BOOST_CHECK(persistence.second.size() == 3);

  //check dimension 0:
  BOOST_CHECK(persistence.second[0].size() == dim0.size());
  for (size_t i = 0; i != persistence.second[0].size(); ++i) {
    BOOST_CHECK((fabs(persistence.second[0][i].first - (double) dim0[i].first) < EPSILON));
    BOOST_CHECK((fabs(persistence.second[0][i].second - (double) dim0[i].second) < EPSILON));
  }

  //check dimension 1:
  BOOST_CHECK(persistence.second[1].size() == dim1.size());
  for (size_t i = 0; i != persistence.second[1].size(); ++i) {
    BOOST_CHECK((fabs(persistence.second[1][i].first - (double) dim1[i].first) < EPSILON));
    BOOST_CHECK((fabs(persistence.second[1][i].second - (double) dim1[i].second) < EPSILON));
  }

  //check dimension 2:
  BOOST_CHECK(persistence.second[2].size() == dim2.size());
  for (size_t i = 0; i != persistence.second[2].size(); ++i) {
    BOOST_CHECK((fabs(persistence.second[2][i].first - (double) dim2[i].first) < EPSILON));
    BOOST_CHECK((fabs(persistence.second[2][i].second - (double) dim2[i].second) < EPSILON));
  }
}
