#include <gudhi/ANN_graph.h>
#include <gudhi/Cluster.h>

#include <iostream>
#include <vector>
#include <algorithm>

using namespace std;
using namespace Gudhi::ANN_graph;
using namespace Gudhi::cluster;

// rename for brevity
typedef vector<ANN_point>::iterator Iterator;
typedef ANN_graph<Iterator> ANN_ngbh_graph;
typedef Cluster<ANN_ngbh_graph> ANN_Cluster;

// main function
int main(int argc, char *argv[]) {
  double func_value = 1.0;
  // create point cloud
  vector<ANN_point> point_cloud;

  vector<double> point_coord = {0, 0, 0};
  int point_dimension = point_coord.size();
  point_cloud.push_back(ANN_point(point_coord, point_dimension, func_value++));
  point_coord = {1, 0, 0};
  point_cloud.push_back(ANN_point(point_coord, point_dimension, func_value++));
  point_coord = {0, 1, 0};
  point_cloud.push_back(ANN_point(point_coord, point_dimension, func_value));
  point_coord = {1, 1, 0};
  point_cloud.push_back(ANN_point(point_coord, point_dimension, func_value++));
  point_coord = {0, 0, 1};
  point_cloud.push_back(ANN_point(point_coord, point_dimension, func_value++));
  point_coord = {1, 0, 1};
  point_cloud.push_back(ANN_point(point_coord, point_dimension, func_value--));
  point_coord = {0, 1, 1};
  point_cloud.push_back(ANN_point(point_coord, point_dimension, func_value--));
  point_coord = {1, 1, 1};
  point_cloud.push_back(ANN_point(point_coord, point_dimension, func_value));

  func_value = 0.5;
  point_coord = {1, 1, 4};
  point_cloud.push_back(ANN_point(point_coord, point_dimension, func_value++));
  point_coord = {1, 1, 4.5};
  point_cloud.push_back(ANN_point(point_coord, point_dimension, func_value--));
  point_coord = {0.5, 1, 4};
  point_cloud.push_back(ANN_point(point_coord, point_dimension, func_value));

  // create distance structure
  double mu = 1;
  ANN_ngbh_graph ngbh_graph(point_cloud.begin(), point_cloud.end(), point_dimension, mu);

  // Let's see the 0-persistence with a really high persistence threshold
  ANN_Cluster zero_persistence(ngbh_graph, 1e20);
  cout << "### 0-persistence intervals:" << endl;
  zero_persistence.output_intervals(std::cout);

  // -----------------------------------------------------------------------------
  // We can see from the traces, we have 2 interesting intervals starting from 0.5
  // -----------------------------------------------------------------------------

  // Let's construct a new cluster with a 0.5 persistence threshold
  ANN_Cluster cluster(ngbh_graph, 0.5);

  cout << "### Number of clusters: " << cluster.get_nb_clusters() << endl;
  cout << "### Points with its cluster value:" << endl;
  for (Iterator pit = point_cloud.begin(); pit < point_cloud.end(); pit++) {
    cout << "Point[" << ngbh_graph.get_xyz(pit) << "] - cluster is :" << cluster.get_cluster(pit) << endl;
  }

}
