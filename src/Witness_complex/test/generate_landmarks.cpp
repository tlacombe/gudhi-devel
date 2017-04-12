#include <sys/types.h>
#include <sys/stat.h>

#include <gudhi/Points_off_io.h>
#include <gudhi/choose_n_farthest_points.h>

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
        << " path_to_point_file number_of_landmarks output_file\n";
    return 0;
  }
  
  std::string file_name = argv[1], out_name = argv[3];
  int nbL = atoi(argv[2]);
  
  // Read the point file
  Point_range witnesses, landmarks;
  Gudhi::Points_off_reader<Point_d> off_reader(file_name);
  if (!off_reader.is_valid()) {
      std::cerr << "Witness complex - Unable to read file " << file_name << "\n";
      exit(-1);  // ----- >>
    }
  witnesses = Point_range(off_reader.get_point_cloud());
  
  // std::cout << "Successfully read " << witnesses.size() << " points.\n";
  // std::cout << "Ambient dimension is " << witnesses[0].dimension() << ".\n";

  // Choose landmarks
  Gudhi::subsampling::choose_n_farthest_points(Kernel(), witnesses, nbL, 0, std::back_inserter(landmarks));

  std::ofstream ofs(out_name, std::ofstream::out);
  ofs << "OFF\n" << nbL << " 0 0\n";
  for (auto l: landmarks) {
    for (auto x: l)
      ofs << x << " ";
    ofs << "\n";
  }
  ofs.close();

  std::cout << "Successfully generated " << nbL << " points.\n";
}
