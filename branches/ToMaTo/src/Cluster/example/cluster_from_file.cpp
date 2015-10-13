/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Primoz Skraba
 *
 *    Copyright (C) 2009 Primoz Skraba.  All Rights Reserved.
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

#include <gudhi/ANN_graph.h>
#include <gudhi/Cluster.h>

#include <iostream>
#include <sstream>
#include <cassert>
#include <vector>  // std::begin, std::end
#include <algorithm>
#include <string>

using namespace std;
using namespace Gudhi::ANN_graph;
using namespace Gudhi::cluster;

// rename for brevity
typedef ANN_graph< vector< ANN_point >::iterator > ANN_ngbh_graph;
typedef Cluster< ANN_ngbh_graph > ANN_Cluster;

// comparison function object for vector indices
template<class V> class Less_Than {
 protected:
  V& v;
 public:
  Less_Than(V& v_) : v(v_) { }

  bool operator()(const int a, const int b) const {
    return ANN_point::Less_Than()(v[a], v[b]);
  }
};


// main function
int main(int argc, char *argv[]) {
  if (argc != 4) {
    cout << "Usage:" << endl << argv[0] << " <input filename> <Rips radius> <prominence threshold>" << endl;
    exit(0);
  }

  int com = 1;
  vector< ANN_point > point_cloud;

  // read in data points with density values
  string input_file_name = argv[com++];
  ifstream input;
  input.open(input_file_name.c_str());
  assert(input.good());

  int dim = -1;
  int nb_points = 0;
  string lineData;
  while (getline(input, lineData)) {
    // read next point's coordinates
    double d;
    vector<double> row;
    stringstream lineStream(lineData);
    while (lineStream >> d)
      row.push_back(d);

    // set up dim if not already done
    if (dim < 0) {
      // last coordinate is reserved for function value
      dim = static_cast<int>(row.size()) - 1;
      cout << "Dimension: " << dim << endl;
    } else if (dim != static_cast<int>(row.size()) - 1) {
      // else check that dimension is preserved
      cerr << "Error: mismatched dimension in " << input_file_name << " at line " << (nb_points + 1) << endl;
      return -1;
    }

    // create vertex from point and function value
    ANN_point v(row, dim, row[dim]);
    point_cloud.push_back(v);
    nb_points++;
  }
  input.close();
  cout << "Number of input points: " << nb_points << endl;

  // sort point cloud and retrieve permutation (for pretty output)
  vector<int> perm;
  perm.reserve(nb_points);
  for (int i = 0; i < nb_points; i++)
    perm.push_back(i);
  std::sort(perm.begin(), perm.end(), Less_Than<vector<ANN_point> >(point_cloud));
  // store inverse permutation as array of iterators on initial point cloud
  vector< vector<ANN_point>::iterator> pperm;
  pperm.reserve(nb_points);
  for (int i = 0; i < nb_points; i++)
    pperm.push_back(point_cloud.begin());
  for (int i = 0; i < nb_points; i++)
    pperm[perm[i]] = (point_cloud.begin() + i);
  // operate permutation on initial point cloud
  vector<ANN_point> pc;
  pc.reserve(nb_points);
  for (int i = 0; i < nb_points; i++)
    pc.push_back(point_cloud[i]);
  for (int i = 0; i < nb_points; i++)
    point_cloud[i] = pc[perm[i]];

  // get rips parameter
  double r = atof(argv[com++]);

  // create distance structure
  ANN_ngbh_graph metric_information(point_cloud.begin(), point_cloud.end(), dim, r * r);

  // create cluster data structure
  ANN_Cluster output_clusters(metric_information, atof(argv[com++]));

  cout << "Nb clusters = " << output_clusters.get_nb_clusters() << endl;
  // output clusters (use permutation to preserve original point order)
  ofstream out;
  out.open("clusters.txt");
  output_clusters.output_clusters(out, pperm.begin(), pperm.end());
  out.close();

  // output barcode
  out.open("diagram.txt");
  output_clusters.output_intervals(out);
  out.close();

  // output colored clusters to COFF file (first 3 dimensions are selected)
  out.open("clusters_3d.coff");
  output_clusters.output_clusters_coff(out);
  out.close();
}
