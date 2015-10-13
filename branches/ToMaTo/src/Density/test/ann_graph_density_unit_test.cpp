#include <cmath> // float comparison
#include <iostream>
#include <vector>
#include <algorithm>
#include <limits>  // for numeric_limits<>

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE "ann_graph"
#include <boost/test/unit_test.hpp>

#include <gudhi/ANN_graph.h>
#include <gudhi/Density.h>

using namespace std;
using namespace Gudhi::ANN_graph;
using namespace Gudhi::density;

//rename for brevity
typedef vector<ANN_point>::iterator Iterator;
typedef ANN_graph<Iterator> ANN_ngbh_graph;
typedef Density<ANN_ngbh_graph> ANN_Density;

BOOST_AUTO_TEST_CASE(density_when_ANN_graph_empty) {
  // create point cloud
  vector<ANN_point> point_cloud;
  // create distance structure
  ANN_ngbh_graph ngbh_graph(point_cloud.begin(), point_cloud.end(), 0, 0);
  
  ANN_Density density(ngbh_graph);
  
  // The aim of this test is a non-crashing test on an empty graph
  
  density.ball_density(0.5);  // in a ball of radius 0.5
  density.gaussian_NN(2, 0.5);  // Gaussian for the 2-closest neighbors, with a height of 0.5
  density.gaussian_mu(0.5, 2);  // Gaussian in a ball of radius 0.5, with a height of 2
  density.distance_to_density(3);  // density for the 3-closest neighbors

  // Check graph has not been modified
  BOOST_CHECK(ngbh_graph.get_num_points() == 0);
  BOOST_CHECK(ngbh_graph.get_start() == ngbh_graph.get_end());
}

bool AreAlmostTheSame(double a, double b) {
  return std::fabs(a - b) < std::numeric_limits<double>::epsilon();
}

void test_function(double func_value, vector<ANN_point>& point_cloud, ANN_ngbh_graph& ngbh_graph) {
  for (auto point_iterator = point_cloud.begin(); point_iterator < point_cloud.end(); point_iterator++) {
    cout << "Point [" << point_iterator->first_3_cordinates() << "] - function = ";
    cout << ngbh_graph.get_func(point_iterator) << endl;
    BOOST_CHECK(AreAlmostTheSame(ngbh_graph.get_func(point_iterator), func_value));
  }
}

BOOST_AUTO_TEST_CASE(simple_ann_graph) {
    // create point cloud
  vector<ANN_point> point_cloud;
  vector<double> point_coord = {0, 0, 0};
  int point_dimension = point_coord.size();
  point_cloud.push_back(ANN_point(point_coord, point_dimension));
  point_coord = {1, 0, 0};
  point_cloud.push_back(ANN_point(point_coord, point_dimension));
  point_coord = {0, 1, 0};
  point_cloud.push_back(ANN_point(point_coord, point_dimension));
  point_coord = {1, 1, 0};
  point_cloud.push_back(ANN_point(point_coord, point_dimension));
  point_coord = {0, 0, 1};
  point_cloud.push_back(ANN_point(point_coord, point_dimension));
  point_coord = {1, 0, 1};
  point_cloud.push_back(ANN_point(point_coord, point_dimension));
  point_coord = {0, 1, 1};
  point_cloud.push_back(ANN_point(point_coord, point_dimension));
  point_coord = {1, 1, 1};
  point_cloud.push_back(ANN_point(point_coord, point_dimension));
  
  // create distance structure
  double mu = 2;
  ANN_ngbh_graph ngbh_graph(point_cloud.begin(), point_cloud.end(), point_dimension, mu);
  
  ANN_Density density(ngbh_graph);
  density.ball_density(0.5);  // in a ball of radius 0.5
  cout << "--- BALL DENSITY\n";
  test_function(-0.25, point_cloud, ngbh_graph);

  density.gaussian_NN(1, 0.5);  // Gaussian for the 1-closest neighbors, with a height of 0.5
  cout << "--- GAUSSIAN NN\n";
  test_function(-1.0, point_cloud, ngbh_graph);

  density.gaussian_mu(0.5, 2);  // Gaussian in a ball of radius 0.5, with a height of 2
  cout << "--- GAUSSIAN MU\n";
  test_function(-1.0, point_cloud, ngbh_graph);

  density.distance_to_density(2);  // density for the 2-closest neighbors
  cout << "--- DISTANCE TO DENSITY\n";
  test_function(sqrt(2), point_cloud, ngbh_graph);
}
