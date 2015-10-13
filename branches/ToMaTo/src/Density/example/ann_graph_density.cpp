#include <gudhi/ANN_graph.h>
#include <gudhi/Density.h>

#include <iostream>
#include <vector>
#include <algorithm>

using namespace std;
using namespace Gudhi::ANN_graph;
using namespace Gudhi::density;

// rename for brevity
typedef vector<ANN_point>::iterator Iterator;
typedef ANN_graph<Iterator> ANN_ngbh_graph;
typedef Density<ANN_ngbh_graph> ANN_Density;

void print_function(const char* func_str, vector<ANN_point>& point_cloud, ANN_ngbh_graph& ngbh_graph) {
  cout << func_str << endl;
  for (auto point_iterator = point_cloud.begin(); point_iterator < point_cloud.end(); point_iterator++) {
    cout << "Point [" << point_iterator->first_3_cordinates() << "] - function = ";
    cout << ngbh_graph.get_func(point_iterator) << endl;
  }
  cout << endl;
}

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

  ANN_Density density(ngbh_graph);
  density.ball_density(0.5);  // in a ball of radius 0.5
  print_function("--- BALL DENSITY", point_cloud, ngbh_graph);

  density.gaussian_NN(2, 0.5);  // Gaussian for the 2-closest neighbors, with a height of 0.5
  print_function("--- GAUSSIAN NN", point_cloud, ngbh_graph);

  density.gaussian_mu(0.5, 2);  // Gaussian in a ball of radius 0.5, with a height of 2
  print_function("--- GAUSSIAN MU", point_cloud, ngbh_graph);

  density.distance_to_density(3);  // density for the 3-closest neighbors
  print_function("--- DISTANCE TO DENSITY", point_cloud, ngbh_graph);
}
