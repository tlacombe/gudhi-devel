//----------------------------------------------------------------------
//----------------------------------------------------------------------
// File:	        main.C
// Programmer:		Primoz Skraba
// Description:		Example persistence clustering pipeline
// Last modified:	Sept. 8, 2009 (Version 0.2)
//----------------------------------------------------------------------
//  Copyright (c) 2009 Primoz Skraba.  All Rights Reserved.
//-----------------------------------------------------------------------
//
//
//    This program is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
//
//-----------------------------------------------------------------------


#include <iostream>
#include <sstream>
#include <cassert>

#include <algorithm>


#include "gudhi/ToMaTo.h"
#include "gudhi/ToMaTo/ANN/ANN_graph.h"
#include "gudhi/ToMaTo/ANN/ANN_point.h"

using namespace std;

//using namespace Gudhi::ToMaTo;

//rename for brevity
typedef Vertex<ANN_point> Point;


// main function

int main(int argc, char *argv[]) {

  if (argc != 4) {
    cout << "Usage:" << endl << argv[0] << " <input filename> <Rips radius> <prominence threshold>" << endl;
    exit(0);
  }


  int com = 1;
  vector< Point > point_cloud;

  //read in data points with density values
  string input_file_name = argv[com++];
  ifstream input;
  input.open(input_file_name.c_str());
  assert(input.good());


  // create distance structure
  ANN_graph< Point > metric_information;

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
      dim = row.size() - 1; // reserve last coordinate for function value
      cout << "Dimension: " << dim << endl;
    }// else check that dimension is preserved
    else if (dim != static_cast<int> (row.size()) - 1) {
      cerr << "Error: mismatched dimension in "
          << input_file_name << " at line " << (nb_points + 1) << endl;
      return -1;
    }

    // create new point and corresponding vertex
    ANN_point p(dim);
    p.coord = new double[dim];
    for (int i = 0; i < dim; i++)
      p.coord[i] = row[i];
    Point v(p);
    v.set_func(row[dim]);
    metric_information.push_back(v);
    nb_points++;
  }
  input.close();
  cout << "Number of input points: " << nb_points << " - dimension = " << dim << endl;

  metric_information.set_dimension(dim);
  //set rips parameter
  double r = atof(argv[com++]);
  metric_information.set_sqrad(r*r);
  
  metric_information.set_persistence_threshold(atof(argv[com++]));

  metric_information.compute_persistence();
  
  //output barcode
  metric_information.output_intervals("simple_diagram.txt");

  //output colored clusters to OFF file (first 3 dimensions are selected)
  metric_information.output_clusters_off("simple_clusters.off");

  // output clusters (use permutation to preserve original point order)
  metric_information.output_clusters("simple_clusters.txt");
}
