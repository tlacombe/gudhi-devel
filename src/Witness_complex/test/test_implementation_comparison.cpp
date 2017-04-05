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

std::string file_name;

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

/* Returns true if and only if the simplex is critical.
 */
template <class SimplexTree>
bool is_critical(typename SimplexTree::Simplex_handle sh,
                 SimplexTree& st)
{
  assert(st.num_simplices() != 0);
  // if (st.cofaces_simplex_range(sh,1).empty())
  //   std::cout << "! ";
  for (auto sh_cof: st.cofaces_simplex_range(sh, 1))
    if (st.filtration(sh_cof) == st.filtration(sh))
      return false;
  return true;
}

/* Returns the number of critical simplices.
 */
template <class SimplexTree>
int num_crit_simplices(SimplexTree& st)
{
  std::ofstream ofs("num_crit.out", std::ofstream::out);
  int count = 0;
  for (auto sh: st.complex_simplex_range())
    if (is_critical(sh, st)) {
      count++;
      ofs << st.dimension(sh);
      for (auto v: st.simplex_vertex_range(sh))
        ofs << " " << v;
      ofs << " " << st.filtration(sh) << "\n";
    }
  ofs.close();
  return count;
}

/* Returns the number of critical simplices. Might be faster.
 */
template <class SimplexTree>
int num_crit_simplices2(SimplexTree& st)
{
  int count = 0;
  double curr_filtr = 0;
  st.initialize_filtration();
  SimplexTree* st_temp = new SimplexTree;
  std::ofstream ofs(file_name, std::ofstream::out);
  for (auto sh: st.filtration_simplex_range()) {
    if (curr_filtr == st.filtration(sh))
      st_temp->insert_simplex(st.simplex_vertex_range(sh));
    else {
      st_temp->set_dimension(st.dimension());
      for (auto sh_temp: st_temp->complex_simplex_range())
        if ((!st_temp->has_children(sh_temp)) && is_critical(sh_temp, *st_temp)) {
          count++;
          ofs << st_temp->dimension(sh_temp);
          for (auto v: st_temp->simplex_vertex_range(sh_temp))
            ofs << " " << v;
          ofs << " " << curr_filtr << "\n";
        }
      curr_filtr = st.filtration(sh);
      delete st_temp;
      st_temp = new SimplexTree;
      st_temp->insert_simplex(st.simplex_vertex_range(sh));
    }
  }
  st_temp->set_dimension(st.dimension());
  for (auto sh_temp: st_temp->complex_simplex_range())
    if ((!st_temp->has_children(sh_temp)) && is_critical(sh_temp, *st_temp)) {
      // for (auto v: st_temp->simplex_vertex_range(sh_temp))
      //   std::cout << v << " ";
      // std::cout << "\n";
      count++;
      ofs << st_temp->dimension(sh_temp);
      for (auto v: st_temp->simplex_vertex_range(sh_temp))
        ofs << " " << v;
      ofs << " " << curr_filtr << "\n";
    }
  delete st_temp;
  ofs.close();
  return count;
}

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
  std::ofstream ofs("st1.out", std::ofstream::out);
  ofs << simplex_tree << "\n";
  ofs.close();
  std::cout << "Witness complex 1 (no cofaces, no witlists) took "
      << static_cast<double>(end - start) / CLOCKS_PER_SEC << " s. \n";
  // start = clock();
  // int crit1_st1 = num_crit_simplices(simplex_tree);
  // end = clock();
  // std::cout << "Number of critical simplices: " << crit1_st1 << ". Time = " << static_cast<double>(end - start) / CLOCKS_PER_SEC << "s.\n";
  // start = clock();
  file_name = "num_crit1.out";
  int crit2_st1 = num_crit_simplices2(simplex_tree);
  end = clock();
  std::cout << "Number of critical simplices: " << crit2_st1 << ". Time = " << static_cast<double>(end - start) / CLOCKS_PER_SEC << "s.\n";

  
  std::cout << "Number of simplices is: " << simplex_tree.num_simplices() << "\n";
  
  // std::cout << simplex_tree << std::endl;
  
  // Compute witness complex - 2
  start = clock();
  Witness_complex_new witness_complex_new(nearest_landmark_table);

  witness_complex_new.create_complex(simplex_tree2, alpha2, lim_dim);
  end = clock();
  std::cout << "Witness complex 2 (cofaces and witlists) took "
      << static_cast<double>(end - start) / CLOCKS_PER_SEC << " s. \n";
  file_name = "num_crit2.out";
  std::cout << "Number of critical simplices: " << num_crit_simplices2(simplex_tree2) << "\n";
  std::cout << "Number of simplices is: " << simplex_tree2.num_simplices() << "\n";
  assert(simplex_tree == simplex_tree2);

  // // Compute witness complex - 3
  start = clock();
  Witness_complex_cof witness_complex_cof(nearest_landmark_table);
  ofs = std::ofstream("st2.out", std::ofstream::out);
  ofs << simplex_tree3 << "\n";
  ofs.close();
  
  witness_complex_cof.create_complex(simplex_tree3, alpha2, lim_dim);
  end = clock();
  std::cout << "Witness complex 3 (cofaces, no witlists) took "
      << static_cast<double>(end - start) / CLOCKS_PER_SEC << " s. \n";
  file_name = "num_crit3.out";
  std::cout << "Number of critical simplices: " << num_crit_simplices2(simplex_tree3) << "\n";
  std::cout << "Number of simplices is: " << simplex_tree3.num_simplices() << "\n";
  assert(simplex_tree == simplex_tree3);
  
  // // Compute witness complex - 4
  // start = clock();
  // Witness_complex_wmap witness_complex_wmap(nearest_landmark_table);

  // witness_complex_wmap.create_complex(simplex_tree4, alpha2, lim_dim);
  // end = clock();
  // std::cout << "Witness complex 4 (no cofaces, but witlists) took "
  //     << static_cast<double>(end - start) / CLOCKS_PER_SEC << " s. \n";
  // std::cout << "Number of critical simplices: " << num_crit_simplices2(simplex_tree4) << "\n";
  // std::cout << "Number of simplices is: " << simplex_tree4.num_simplices() << "\n";
  // assert(simplex_tree == simplex_tree4);

  // // Compute witness complex - SAL 1
  // start = clock();
  // // Witness_complex witness_complex(nearest_landmark_table);

  // witness_complex.create_complex(sal1, alpha2, lim_dim);
  // end = clock();
  // std::cout << "Witness complex SAL 1 (no cofaces, no witlists) took "
  //     << static_cast<double>(end - start) / CLOCKS_PER_SEC << " s. \n";
  // std::cout << "Number of critical simplices is: " << sal1.num_simplices() << "\n";
  // // Compute the total number of simplices with ST
  // Gudhi::Simplex_tree<> st_temp1;
  // for (auto sh: sal1.critical_cofaces(Gudhi::Simplex()) ) {
  //   st_temp1.insert_simplex_and_subfaces(*sh);
  // }
  // std::cout << "Number of simplices is: " << st_temp1.num_simplices() << "\n";  
  
  // // Compute witness complex - SAL 4
  // start = clock();
  // Witness_complex_sal4 witness_complex_sal4(nearest_landmark_table);

  // witness_complex_sal4.create_complex(sal4, alpha2, lim_dim);
  // end = clock();
  // std::cout << "Witness complex SAL 4 (no cofaces, but witlists) took "
  //     << static_cast<double>(end - start) / CLOCKS_PER_SEC << " s. \n";
  // std::cout << "Number of critical simplices is: " << sal4.num_simplices() << "\n";

  // Gudhi::Simplex_tree<> st_temp4;
  // for (auto sh: sal1.critical_cofaces(Gudhi::Simplex()) ) {
  //   st_temp4.insert_simplex_and_subfaces(*sh);
  // }
  // std::cout << "Number of simplices is: " << st_temp4.num_simplices() << "\n";  

  
}
