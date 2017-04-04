#include <sys/types.h>
#include <sys/stat.h>

#include <gudhi/Simplex_tree.h>
#include <gudhi/SAL.h>
// #include <gudhi/SALW.h>

#include <gudhi/Witness_complex.h>
#include <gudhi/Witness_complex_new.h>
#include <gudhi/Witness_complex_cof.h>
#include <gudhi/Witness_complex_wmap.h>
//#include <gudhi/Witness_complex_sal2.h> 
#include <gudhi/Witness_complex_sal4.h>
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
typedef typename Gudhi::witness_complex::Witness_complex_cof<Nearest_landmark_table> Witness_complex_cof;
typedef typename Gudhi::witness_complex::Witness_complex_wmap<Nearest_landmark_table> Witness_complex_wmap;

typedef typename Gudhi::witness_complex::Witness_complex_sal4<Nearest_landmark_table> Witness_complex_sal4;


int main(int argc, char * const argv[]) {
  if (argc != 5) {
    std::cerr << "Usage: " << argv[0]
        << " path_to_point_file number_of_landmarks max_squared_alpha limit_dimension\n";
    return 0;
  }
  
  std::string file_name = argv[1];
  int nbL = atoi(argv[2]), lim_dim = atoi(argv[4]);
  double alpha2 = atof(argv[3]);
  clock_t start, end;
  Gudhi::Simplex_tree<> simplex_tree, simplex_tree2, simplex_tree3, simplex_tree4;
  Gudhi::SAL sal1, sal2, sal3, sal4;
  
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

  
  // Compute witness complex - 1
  start = clock();
  Witness_complex witness_complex(nearest_landmark_table);

  witness_complex.create_complex(simplex_tree, alpha2, lim_dim);
  end = clock();
  std::cout << "Witness complex 1 (no cofaces, no witlists) took "
      << static_cast<double>(end - start) / CLOCKS_PER_SEC << " s. \n";
  std::cout << "Number of simplices is: " << simplex_tree.num_simplices() << "\n";
                                                                             
  // std::cout << simplex_tree << std::endl;
  
  // // Compute witness complex - 2
  // start = clock();
  // Witness_complex_new witness_complex_new(nearest_landmark_table);

  // witness_complex_new.create_complex(simplex_tree2, alpha2, lim_dim);
  // end = clock();
  // std::cout << "Witness complex 2 (cofaces and witlists) took "
  //     << static_cast<double>(end - start) / CLOCKS_PER_SEC << " s. \n";
  // std::cout << "Number of simplices is: " << simplex_tree2.num_simplices() << "\n";

  // // Compute witness complex - 3
  // start = clock();
  // Witness_complex_cof witness_complex_cof(nearest_landmark_table);

  // witness_complex_cof.create_complex(simplex_tree3, alpha2, lim_dim);
  // end = clock();
  // std::cout << "Witness complex 3 (cofaces, no witlists) took "
  //     << static_cast<double>(end - start) / CLOCKS_PER_SEC << " s. \n";
  // std::cout << "Number of simplices is: " << simplex_tree3.num_simplices() << "\n";
  
  // // Compute witness complex - 4
  // start = clock();
  // Witness_complex_wmap witness_complex_wmap(nearest_landmark_table);

  // witness_complex_wmap.create_complex(simplex_tree4, alpha2, lim_dim);
  // end = clock();
  // std::cout << "Witness complex 4 (no cofaces, but witlists) took "
  //     << static_cast<double>(end - start) / CLOCKS_PER_SEC << " s. \n";
  // std::cout << "Number of simplices is: " << simplex_tree4.num_simplices() << "\n";

  // Compute witness complex - SAL 1
  start = clock();
  // Witness_complex witness_complex(nearest_landmark_table);

  witness_complex.create_complex(sal1, alpha2, lim_dim);
  end = clock();
  std::cout << "Witness complex SAL 1 (no cofaces, no witlists) took "
      << static_cast<double>(end - start) / CLOCKS_PER_SEC << " s. \n";
  std::cout << "Number of critical simplices is: " << sal1.num_simplices() << "\n";
  // Compute the total number of simplices with ST
  Gudhi::Simplex_tree<> st_temp1;
  for (auto sh: sal1.critical_cofaces(Gudhi::Simplex()) ) {
    st_temp1.insert_simplex_and_subfaces(*sh);
  }
  std::cout << "Number of simplices is: " << st_temp1.num_simplices() << "\n";  
  
  // Compute witness complex - SAL 4
  start = clock();
  Witness_complex_sal4 witness_complex_sal4(nearest_landmark_table);

  witness_complex_sal4.create_complex(sal4, alpha2, lim_dim);
  end = clock();
  std::cout << "Witness complex SAL 4 (no cofaces, but witlists) took "
      << static_cast<double>(end - start) / CLOCKS_PER_SEC << " s. \n";
  std::cout << "Number of critical simplices is: " << sal4.num_simplices() << "\n";

  Gudhi::Simplex_tree<> st_temp4;
  for (auto sh: sal1.critical_cofaces(Gudhi::Simplex()) ) {
    st_temp4.insert_simplex_and_subfaces(*sh);
  }
  std::cout << "Number of simplices is: " << st_temp4.num_simplices() << "\n";  

  
}
