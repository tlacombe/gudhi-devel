#include <gudhi/ANN_graph.h>

#include <iostream>
#include <vector>
#include <algorithm>

using namespace std;
using namespace Gudhi::ANN_graph;

// rename for brevity
typedef vector<ANN_point>::iterator Iterator;
typedef ANN_graph<Iterator> ANN_ngbh_graph;

// main function
int main(int argc, char *argv[]) {
  // create point cloud
  vector<ANN_point> point_cloud;
  vector<double> point_coord = {0.1, 0.2};
  int point_dimension = point_coord.size();
  point_cloud.push_back(ANN_point(point_coord, point_dimension));
  point_coord = {0.2, 0.1};
  point_cloud.push_back(ANN_point(point_coord, point_dimension));
  point_coord = {0.15, 0.15};
  point_cloud.push_back(ANN_point(point_coord, point_dimension));
  point_coord = {12, 20};
  point_cloud.push_back(ANN_point(point_coord, point_dimension));

  // create distance structure
  double mu = 1e20;
  ANN_ngbh_graph ngbh_graph(point_cloud.begin(), point_cloud.end(), point_dimension, mu);

  for (auto point_iterator = point_cloud.begin(); point_iterator < point_cloud.end(); point_iterator++) {
    vector<Iterator> vi;
    // get_neighbors inits a vector of iterator on neighbors points except the given point
    ngbh_graph.get_neighbors(point_iterator, vi);
    cout << "Point [" << point_iterator->first_3_cordinates() << "] - number of neighbors in a ball of radius ";
    cout << mu << " = " << vi.size() << endl;
  }
  cout << endl;

  double radius_1 = 0.5;
  for (auto point_iterator = point_cloud.begin(); point_iterator < point_cloud.end(); point_iterator++) {
    // get_num_neighbors returns the number of neighbors including the given point
    cout << "Point [" << point_iterator->first_3_cordinates() << "] - has ";
    cout << ngbh_graph.get_num_neighbors(point_iterator, radius_1) << " neighbors in a ball of radius " << radius_1;
    cout << endl;
  }
  cout << endl;

  int k = 3;
  for (auto point_iterator = point_cloud.begin(); point_iterator < point_cloud.end(); point_iterator++) {
    double* ndist = new double[k];
    // get_neighbors_dist inits an array of distance with k-neighbors points including the given point
    ngbh_graph.get_neighbors_dist(point_iterator, k, ndist);
    cout << "Point [" << point_iterator->first_3_cordinates() << "] - " << k << "-closest neighbors distances are:";
    for (int i = 0; i < k; i++)
      cout << " " << ndist[i];
    cout << endl;
    delete[] ndist;
  }
  cout << endl;

  double radius_2 = 1.0;
  for (auto point_iterator = point_cloud.begin(); point_iterator < point_cloud.end(); point_iterator++) {
    double* ndist = new double[point_cloud.size()];
    // get_neighbors_dist_r inits an array of distance with neighbors points in a ball of a given radius including the
    // given point
    int nb_points = ngbh_graph.get_neighbors_dist_r(point_iterator, radius_2, ndist);
    cout << "Point [" << point_iterator->first_3_cordinates() << "] - closest neighbors, in a ball of radius ";
    cout << radius_2 << ", distances are:";
    for (int i = 0; i < nb_points; i++)
      cout << " " << ndist[i];
    cout << endl;
    delete[] ndist;
  }
  cout << endl;
}
