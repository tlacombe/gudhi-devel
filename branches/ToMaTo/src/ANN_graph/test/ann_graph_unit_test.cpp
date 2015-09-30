#include <cmath> // float comparison
#include <iostream>
#include <vector>
#include <algorithm>
#include <limits>  // for numeric_limits<>

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE "ann_graph"
#include <boost/test/unit_test.hpp>

#include <gudhi/ANN_graph.h>

using namespace std;
using namespace Gudhi::ANN_graph;

//rename for brevity
typedef vector<ANN_point>::iterator Iterator;
typedef ANN_graph<Iterator> ANN_ngbh_graph;

BOOST_AUTO_TEST_CASE(ANN_graph_when_empty) {
  // create point cloud
  vector<ANN_point> point_cloud;
  // create distance structure
  ANN_ngbh_graph ngbh_graph(point_cloud.begin(), point_cloud.end(), 0, 0);
  BOOST_CHECK(ngbh_graph.get_num_points() == 0);
  BOOST_CHECK(ngbh_graph.get_start() == ngbh_graph.get_end());
}

bool AreAlmostTheSame(double a, double b) {
  return std::fabs(a - b) < std::numeric_limits<double>::epsilon();
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
  
  for (auto point_iterator = point_cloud.begin(); point_iterator < point_cloud.end(); point_iterator++) {
    vector<Iterator> vi;
    ngbh_graph.get_neighbors(point_iterator, vi);
    cout << "Point [" << point_iterator->first_3_cordinates() << "] - number of neighbors in a ball of squared radius ";
    cout << mu << " = " << vi.size() << endl;
    // get_neighbors inits a vector of iterator on neighbors points except the given point
    // In the cube of edge 1 case, when mu (squared radius) is 2, all the points are neighbors, except the opposite
    // point on the cube
    BOOST_CHECK(vi.size() == 6);
  }
  cout << endl;

  double radius_1 = 1;
  for (auto point_iterator = point_cloud.begin(); point_iterator < point_cloud.end(); point_iterator++) {
    int nb_neighbor = ngbh_graph.get_num_neighbors(point_iterator, radius_1);
    cout << "Point [" << point_iterator->first_3_cordinates() << "] - has ";
    cout << nb_neighbor << " neighbors in a ball of radius " << radius_1 << endl;
    // get_num_neighbors returns the number of neighbors including the given point
    // In the cube of edge 1 case, when radius is 1, all the points are neighbors, except the 4 opposite
    // points on the cube
    BOOST_CHECK(nb_neighbor == 4);
  }
  cout << endl;
  
  int k = 4;
  for (auto point_iterator = point_cloud.begin(); point_iterator < point_cloud.end(); point_iterator++) {
    double* ndist = new double[k];
    // get_neighbors_dist inits an array of distance with k-neighbors points including the given point
    // In the cube of edge 1 case, when k is 4, all the points are neighbors, except the 4 opposite
    // points on the cube
    ngbh_graph.get_neighbors_dist(point_iterator, k, ndist);
    cout << "Point [" << point_iterator->first_3_cordinates() << "] - " << k << "-closest neighbors distances are:";
    for (int i = 0; i < k; i++) {
      cout << " " << ndist[i];
      if (i == 0)
        BOOST_CHECK(AreAlmostTheSame(ndist[i], 0.0));  // The given point
      else
        BOOST_CHECK(AreAlmostTheSame(ndist[i], 1.0));
    }
    cout << endl;
    delete[] ndist;
  }
  cout << endl;
  
  double radius_2 = 1.415;  // A little bit more than 2^(1/2)
  for (auto point_iterator = point_cloud.begin(); point_iterator < point_cloud.end(); point_iterator++) {
    double* ndist = new double[point_cloud.size()];
    // get_neighbors_dist_r inits an array of distance with neighbors points in a ball of a given radius including the
    // given point
    // In the cube of edge 1 case, when radius is 2^(1/2), all the points are neighbors, except the opposite
    // point on the cube
    int nb_neighbor = ngbh_graph.get_neighbors_dist_r(point_iterator, radius_2, ndist);
    BOOST_CHECK(nb_neighbor == 7);
    cout << "Point [" << point_iterator->first_3_cordinates() << "] - closest neighbors, in a ball of radius ";
    cout << radius_2 << ", distances are:";
    for (int i = 0; i < nb_neighbor; i++) {
      cout << " " << ndist[i];
      if (i == 0)
        BOOST_CHECK(AreAlmostTheSame(ndist[i], 0.0));  // The given point
      else if ((i > 0) && (i < 4))
        BOOST_CHECK(AreAlmostTheSame(ndist[i], 1.0));
      else
        BOOST_CHECK(AreAlmostTheSame(ndist[i], 2.0));
    }
    cout << endl;
    delete[] ndist;
  }
  cout << endl;

}
