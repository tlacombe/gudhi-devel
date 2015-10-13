/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Vincent Rouvreau
 *
 *    Copyright (C) 2015  INRIA Saclay (France)
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

#include <cmath> // float comparison
#include <iostream>
#include <vector>
#include <algorithm>
#include <limits>  // for numeric_limits<>

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE "ann_graph"
#include <boost/test/unit_test.hpp>

#include <gudhi/ANN_graph.h>
#include <gudhi/Cluster.h>

using namespace std;
using namespace Gudhi::ANN_graph;
using namespace Gudhi::cluster;

//rename for brevity
typedef vector<ANN_point>::iterator Iterator;
typedef ANN_graph<Iterator> ANN_ngbh_graph;
typedef Cluster<ANN_ngbh_graph> ANN_Cluster;

BOOST_AUTO_TEST_CASE(cluster_when_ANN_graph_empty) {
  // create point cloud
  vector<ANN_point> point_cloud;
  // create distance structure
  ANN_ngbh_graph ngbh_graph(point_cloud.begin(), point_cloud.end(), 0, 0);
  
  ANN_Cluster cluster(ngbh_graph, 0.0);
  
  // The aim of this test is a non-crashing test on an empty graph
  
  cluster.output_intervals(std::cout);

  vector<Iterator> pperm;
  cluster.output_clusters(std::cout, pperm.begin(), pperm.end());
  cluster.output_clusters_coff(std::cout);

  // Check graph has not been modified
  BOOST_CHECK(ngbh_graph.get_num_points() == 0);
  BOOST_CHECK(ngbh_graph.get_start() == ngbh_graph.get_end());
}

BOOST_AUTO_TEST_CASE(simple_ann_graph) {
  double func_value = 1.0;
  // create point cloud
  vector<ANN_point> point_cloud;
  
  // 1st cluster
  vector<double> point_coord = {0, 0, 0};
  int point_dimension = point_coord.size();
  point_cloud.push_back(ANN_point(point_coord, point_dimension, func_value));
  point_coord = {1, 0, 0};
  point_cloud.push_back(ANN_point(point_coord, point_dimension, func_value));
  point_coord = {0, 1, 0};
  point_cloud.push_back(ANN_point(point_coord, point_dimension, func_value));
  point_coord = {1, 1, 0};
  point_cloud.push_back(ANN_point(point_coord, point_dimension, func_value));
  point_coord = {0, 0, 1};
  point_cloud.push_back(ANN_point(point_coord, point_dimension, func_value));
  point_coord = {1, 0, 1};
  point_cloud.push_back(ANN_point(point_coord, point_dimension, func_value));
  point_coord = {0, 1, 1};
  point_cloud.push_back(ANN_point(point_coord, point_dimension, func_value));
  point_coord = {1, 1, 1};
  point_cloud.push_back(ANN_point(point_coord, point_dimension, func_value));

  // 2nd cluster
  point_coord = {1, 1, 4};
  point_cloud.push_back(ANN_point(point_coord, point_dimension, func_value));
  point_coord = {1, 1, 4.5};
  point_cloud.push_back(ANN_point(point_coord, point_dimension, func_value));
  point_coord = {0.5, 1, 4};
  point_cloud.push_back(ANN_point(point_coord, point_dimension, func_value));

  // 3rd cluster
  point_coord = {3, 1, 1};
  point_cloud.push_back(ANN_point(point_coord, point_dimension, func_value));
  point_coord = {3.5, 1, 1};
  point_cloud.push_back(ANN_point(point_coord, point_dimension, func_value));
  point_coord = {3, 1.1, 1.5};
  point_cloud.push_back(ANN_point(point_coord, point_dimension, func_value));
  point_coord = {3.2, 0.9, 0.9};
  point_cloud.push_back(ANN_point(point_coord, point_dimension, func_value));
  point_coord = {3, 1.4, 1.2};
  point_cloud.push_back(ANN_point(point_coord, point_dimension, func_value));
  
  // create distance structure
  double mu = 2;
  ANN_ngbh_graph ngbh_graph(point_cloud.begin(), point_cloud.end(), point_dimension, mu);
  
  ANN_Cluster cluster(ngbh_graph, 1);
  cout << "### Number of clusters: " << cluster.get_nb_clusters() << endl;
  BOOST_CHECK(3 == cluster.get_nb_clusters());
  cout << "### Points with its cluster value:" << endl;
  for (Iterator pit = point_cloud.begin(); pit < point_cloud.end(); pit++) {
    cout << "Point[" << ngbh_graph.get_xyz(pit) << "] - cluster is :" << cluster.get_cluster(pit) << endl;
    if ((pit - point_cloud.begin()) < 8) {
      BOOST_CHECK(0 == cluster.get_cluster(pit));
    } else if ((pit - point_cloud.begin()) < 11) {
      BOOST_CHECK(1 == cluster.get_cluster(pit));
    } else {
      BOOST_CHECK(2 == cluster.get_cluster(pit));
    }
  }
}
