#include <sys/types.h>
#include <sys/stat.h>

#include <gudhi/Simplex_tree.h>
#include <gudhi/Witness_complex.h>
#include <gudhi/Witness_complex_new.h>
#include <gudhi/choose_n_farthest_points.h>
#include <gudhi/Kd_tree_search.h>
#include <gudhi/Points_off_io.h>

#include <CGAL/Epick_d.h>

#include <iostream>
#include <fstream>
#include <ctime>
#include <string>
#include <vector>

typedef CGAL::Epick_d<CGAL::Dynamic_dimension_tag> Kernel;
typedef typename Kernel::Point_d Point_d;
typedef std::vector<Point_d> Point_range;
typedef Gudhi::spatial_searching::Kd_tree_search<Kernel, Point_range> Kd_tree;
typedef Kd_tree::INS_range Nearest_landmark_range; 
typedef std::vector<Nearest_landmark_range> Nearest_landmark_table;

typedef typename Gudhi::witness_complex::Witness_complex<Nearest_landmark_table> Witness_complex;
typedef typename Gudhi::witness_complex::Witness_complex_new<Nearest_landmark_table> Witness_complex_new;

int main(int argc, char * const argv[]) {
  if (argc != 5) {
    std::cerr << "Usage: " << argv[0]
        << " path_to_point_file number_of_landmarks max_squared_alpha limit_dimension\n";
    return 0;
  }

  #ifndef DEBUG_TRACES
  #define DEBUG_TRACES
  #endif
  
  std::string file_name = argv[1];
  int nbL = atoi(argv[2]), lim_dim = atoi(argv[4]);
  double alpha2 = atof(argv[3]);
  clock_t start, end;
  Gudhi::Simplex_tree<> simplex_tree, simplex_tree2;

  // Read the point file
  Point_range witnesses, landmarks;
  Gudhi::Points_off_reader<Point_d> off_reader(file_name);
  if (!off_reader.is_valid()) {
      std::cerr << "Witness complex - Unable to read file " << file_name << "\n";
      exit(-1);  // ----- >>
    }
  witnesses = Point_range(off_reader.get_point_cloud());
  
  std::cout << "Successfully read " << witnesses.size() << " points.\n";
  std::cout << "Ambient dimension is " << witnesses[0].dimension() << ".\n";

  // Choose landmarks
  Gudhi::subsampling::choose_n_farthest_points(Kernel(), witnesses, nbL, 0, std::back_inserter(landmarks));

  // Compute nearest neighbor table
  Kd_tree landmark_tree(landmarks);
  Nearest_landmark_table nearest_landmark_table;
  for (auto w: witnesses)
    nearest_landmark_table.push_back(landmark_tree.query_incremental_nearest_neighbors(w));

  
  // Compute witness complex
  start = clock();
  Witness_complex witness_complex(nearest_landmark_table);

  witness_complex.create_complex(simplex_tree, alpha2, lim_dim);
  end = clock();
  std::cout << "Witness complex took "
      << static_cast<double>(end - start) / CLOCKS_PER_SEC << " s. \n";
  std::cout << "Number of simplices is: " << simplex_tree.num_simplices() << "\n";

  // Compute witness complex - 2
  start = clock();
  Witness_complex_new witness_complex_new(nearest_landmark_table);

  witness_complex_new.create_complex(simplex_tree2, alpha2, lim_dim);
  end = clock();
  std::cout << "Witness complex 2 took "
      << static_cast<double>(end - start) / CLOCKS_PER_SEC << " s. \n";
  std::cout << "Number of simplices is: " << simplex_tree2.num_simplices() << "\n";
  
}
