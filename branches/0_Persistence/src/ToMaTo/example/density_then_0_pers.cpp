#include <iostream>

#include "gudhi/ToMaTo.h"
#include "gudhi/ToMaTo/ANN/ANN_point.h"
#include "gudhi/ToMaTo/ANN/ANN_graph_with_density.h"
#include "gudhi/ToMaTo/ANN/ANN_graph_off_io.h"

typedef Vertex<ANN_point> Point;
typedef ANN_graph_with_density< Point > Tree;

using namespace Gudhi;
using namespace ToMaTo;

int main(int argc, char *argv[]) {
  std::string offFileName("2pointscloud.off");
  
  Tree ANN_tree;
  // Set the ANN tree search in a squared radius of 0.25
  ANN_tree.set_sqrad(0.25*0.25);
  // Set the minimum persistence value to 2
  ANN_tree.set_persistence_threshold(2);

  ANN_graph_off_reader<Tree> off_reader(offFileName, ANN_tree);
  // Compute density from the 200th nearest neighbors
  ANN_tree.distance_to_density(200);
  
  ANN_tree.compute_persistence();

  //output barcode in a file
  ANN_tree.output_intervals("2points_diag.txt");
  
  //output colored clusters to COFF file (first 3 dimensions are selected)
  ANN_tree.output_clusters_off("2points_cluster.off");

  // output clusters (use permutation to preserve original point order)
  ANN_tree.output_clusters("2points_cluster.txt");

}
