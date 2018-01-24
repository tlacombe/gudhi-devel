#include <iostream>
#include <vector>

#include <gudhi/Points_off_io.h>
#include <gudhi/Coxeter_system.h>
#include <gudhi/Coxeter_complex.h>
#include <gudhi/Coxeter_complex/Off_point_range.h>

#include <CGAL/Epick_d.h>

//#include <Eigen/Dense>

#include "cxx-prettyprint/prettyprint.hpp"

using K = CGAL::Epick_d<CGAL::Dynamic_dimension_tag>;
using FT = K::FT;
using Point_d = K::Point_d;
using Point_vector = std::vector< Point_d >;
using Coxeter_complex = Gudhi::Coxeter_complex<Point_vector, Coxeter_system>;

std::vector<FT> bounding_box_dimensions(Point_vector& points) {
  std::vector<FT> lower, upper, difference;
  for (auto x: points[0]) {
    lower.push_back(x);
    upper.push_back(x);
  }
  for (auto p: points)
    for (unsigned i = 0; i < p.size(); i++) {
      if (p[i] < lower[i])
        lower[i] = p[i];
      if (p[i] > upper[i])
        upper[i] = p[i];
    }
  for (unsigned i = 0; i < lower.size(); i++)
    difference.push_back(upper[i]-lower[i]);
  return difference;
}


/** Current state of the algorithm.
 *  Input: a point cloud 'point_vector'
 *  Output: a reconstruction (a simplicial complex?, a Czech-like complex?)
 */

int main(int argc, char * const argv[]) {
  std::cout << "Marching cube adaptation for Coxeter triangulations\n";
  if (argc > 3 || argc < 2) {
    std::cerr << "Usage: " << argv[0]
        << " path_to_off_point_file [initial level]\n";
    return 0;
  }
  Point_vector point_vector;
  double init_level = 1;
  if (argc == 3)
    init_level = atof(argv[2]);
  // Gudhi::Points_off_reader<Point_d> off_reader(argv[1]);
  // if (!off_reader.is_valid()) {
  //     std::cerr << "Coxeter triangulations - Unable to read file " << argv[1] << "\n";
  //     exit(-1);  // ----- >>
  //   }
  // point_vector = Point_vector(off_reader.get_point_cloud());
  // int N = point_vector.size();
  // unsigned short d = point_vector[0].size();
  // // short d = 2;
  // std::cout << "Successfully read " << N << " points in dimension " << d << std::endl;
  Gudhi::Off_point_range<Point_d> off_range(argv[1]);
  using Coxeter_complex_off = Gudhi::Coxeter_complex<Gudhi::Off_point_range<Point_d>, Coxeter_system>;
  int d = off_range.dimension();
  {
    Coxeter_system cs_A('A', d);
    Coxeter_complex_off cc(off_range, cs_A, init_level);  
    cc.write_mesh("sphere_coxeter_complex_A.mesh");
    cc.collapse();
  }
  // {
  //   Coxeter_system cs_B('B', d);
  //   Coxeter_complex cc(point_vector, cs_B, init_level);
  //   cc.write_mesh("sphere_coxeter_complex_B.mesh");
  // }
  // {
  //   Coxeter_system cs_C('C', d);
  //   Coxeter_complex cc(point_vector, cs_C, init_level); 
  //   cc.write_mesh("sphere_coxeter_complex_C.mesh");
  // }
 // Coxeter_system cs_D('D', d);
  // Coxeter_complex(point_vector, cs_D);  
  // Coxeter_system cs_E6('E', 6);
  // cs_E6.emplace_back('A', d-6);
  // Coxeter_complex(point_vector, cs_E6);  

  
}