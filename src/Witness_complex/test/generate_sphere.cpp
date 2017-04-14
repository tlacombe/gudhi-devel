#include <sys/types.h>
#include <sys/stat.h>

#include <gudhi/Points_off_io.h>
#include <gudhi/choose_n_farthest_points.h>
#include "../example/generators.h"

#include <CGAL/Epick_d.h>

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

typedef CGAL::Epick_d<CGAL::Dynamic_dimension_tag> Kernel;
typedef typename Kernel::Point_d Point_d;
typedef std::vector<Point_d> Point_range;

int main(int argc, char * const argv[]) {
  if (argc != 4) {
    std::cerr << "Usage: " << argv[0]
        << " dimension number_of_points output_file\n";
    return 0;
  }
  
  std::string out_name = argv[3];
  int dim = atoi(argv[1]), nbW = atoi(argv[2]);
  
  // Read the point file
  Point_range witnesses;
  
  // std::cout << "Successfully read " << witnesses.size() << " points.\n";
  // std::cout << "Ambient dimension is " << witnesses[0].dimension() << ".\n";

  // Generate uniformly random points on the (dim-1)-sphere
  generate_points_sphere(witnesses, nbW, dim);
  
  std::ofstream ofs(out_name, std::ofstream::out);
  ofs << "OFF\n" << nbW << " 0 0\n";
  for (auto w: witnesses) {
    for (auto x: w)
      ofs << x << " ";
    ofs << "\n";
  }
  ofs.close();

  std::cout << "Successfully generated " << nbW << " points.\n";
}
