#include <iostream>
#include <fstream>
#include <chrono>
#include <gudhi/Simplex_tree.h>
#include <gudhi/Zigzag_persistence.h>
#include "gudhi/reader_utils.h"
#include <gudhi/distance_functions.h>
#include <gudhi/Zigzag_filtration.h>
#include <gudhi/Points_off_io.h>
#include <boost/program_options.hpp>
#include <CGAL/Epick_d.h>
#include <gudhi/choose_n_farthest_points.h>
#include <gudhi/pick_n_random_points.h>

// Types definition
using Simplex_tree = Gudhi::Simplex_tree<Gudhi::Simplex_tree_options_zigzag_persistence>;
using Zz_persistence = Gudhi::zigzag_persistence::Zigzag_persistence<Simplex_tree> ;
using Zz_edge = Zigzag_edge<Simplex_tree>;
using Filtration_value = Simplex_tree::Filtration_value;
// using Zigzag_persistence = Gudhi::zigzag_persistence::Zigzag_persistence<Simplex_tree >;
// using Point = std::vector<double>;
using K = CGAL::Epick_d<CGAL::Dynamic_dimension_tag>;
using Point_d = typename K::Point_d;
using Points_off_reader = Gudhi::Points_off_reader<Point_d>;

void program_options( int argc, char* argv[]
                    , std::string& off_file_points
                    , Filtration_value& nu
                    , Filtration_value& mu
                    , int& dim_max);

bool lex_cmp( Point_d p, Point_d q ) {
  auto itp = p.begin(); auto itq = q.begin();
  while(itp != p.end() && itq != q.end()) {
    if(*itp != *itq) { return *itp < *itq; }
    ++itp; ++itq;
  }
  return false;
}

int main(int argc, char* argv[])
{
  std::chrono::time_point<std::chrono::system_clock> start, end;
  int enlapsed_sec;

  std::string off_file_points;
  Filtration_value nu, mu;
  int dim_max;

  program_options(argc, argv, off_file_points, nu, mu, dim_max);
//sequence of insertion and deletions of vertices and edges
  std::vector< Zz_edge >        edge_filtration;
//epsilon_i values, size() == #points
  std::vector<Filtration_value> filtration_values;
  //extract points from file
  Points_off_reader off_reader(off_file_points); //read points
//kernel
  K k_d; 

//check whether there are duplicates
  bool(*tmp_cmp)(Point_d,Point_d) = lex_cmp;
  std::set<Point_d, bool(*)(Point_d,Point_d) > no_dup(tmp_cmp);
  for(auto p : off_reader.get_point_cloud()) { no_dup.insert(p); }
  if(no_dup.size() != off_reader.get_point_cloud().size()) { 
    std::cout << "Duplicates " << no_dup.size() << " vs. " << off_reader.get_point_cloud().size() << "\n"; 
    return 0; 
  }

  //sort points
  start = std::chrono::system_clock::now();
  std::vector<Point_d> sorted_points;
  Gudhi::subsampling::choose_n_farthest_points( k_d, off_reader.get_point_cloud() 
    , off_reader.get_point_cloud().size() //all points
    , 0//start with point [0]//Gudhi::subsampling::random_starting_point
    , std::back_inserter(sorted_points));

  end = std::chrono::system_clock::now();
  enlapsed_sec =std::chrono::duration_cast<std::chrono::seconds>(end-start).count();

  std::cout << "Furthest point sort: " << enlapsed_sec << " sec.\n";

  auto sqdist = k_d.squared_distance_d_object();

  //Compute edge filtration
  start = std::chrono::system_clock::now();
  points_to_edge_filtration( sorted_points, 
                             sqdist, //Gudhi::Euclidean_distance(),
                             nu*nu, mu*mu, filtration_values, edge_filtration );
  end = std::chrono::system_clock::now();
  enlapsed_sec =std::chrono::duration_cast<std::chrono::seconds>(end-start).count();

  std::cout << "Edge filtration computation: " << enlapsed_sec << " sec.\n";
  
  //apply sqrt to correct the use of sqrt distance
  for(auto & f : filtration_values) { f = std::sqrt(f); }
  for(auto & e : edge_filtration) { e.assign_fil(std::sqrt(e.fil())); }

  Simplex_tree st;
  st.initialize_filtration(edge_filtration, dim_max); 
  // auto zz_rg = st.filtration_simplex_range();

  start = std::chrono::system_clock::now();
  Zz_persistence zz(st);
  zz.compute_zigzag_persistence();
  end = std::chrono::system_clock::now();
  enlapsed_sec =std::chrono::duration_cast<std::chrono::seconds>(end-start).count();

  std::cout << "Compute zigzag persistence in: " << enlapsed_sec << " sec.\n";
  std::cout << std::endl;
  std::cout << "Persistence diagram (log2): \n";
  zz.output_log2_diagram();

  std::cout << "\n\n\n";
  std::cout << "Persistence diagram: \n";
  zz.output_diagram();

  return 0;
}

 
//program options 
void program_options(int argc, char* argv[], std::string& off_file_points, Filtration_value& nu, Filtration_value &mu, int& dim_max) {
  namespace po = boost::program_options;
  po::options_description hidden("Hidden options");
  hidden.add_options()("input-file", po::value<std::string>(&off_file_points),
                       "Name of an OFF file containing a point set.\n");

  po::options_description visible("Allowed options", 100);
  visible.add_options()("help,h", "produce help message")
  ( "nu-factor,nu",
    po::value<Filtration_value>(&nu)->default_value(3.0),
    "Lower multiplicative factor in the oscillating Rips zigzag filtration.")
  ( "mu-factor,mu",
    po::value<Filtration_value>(&mu)->default_value(3.2),
    "Upper multiplicative factor in the oscillating Rips zigzag filtration.")
  ( "cpx-dimension,d", po::value<int>(&dim_max)->default_value(1),
    "Maximal dimension of the oscillating Rips complexes in the filtration.")
  // ( "field-charac,p", po::value<int>(&p)->default_value(11),
  //   "Characteristic p of the coefficient field Z/pZ for computing homology.")
  // ( "min-persistence,m", po::value<Filtration_value>(&min_persistence),
  //   "Minimal lifetime of homology feature to be recorded. Default is 0. Enter a negative value to see zero length "
  //     "intervals")
  ;

  po::positional_options_description pos;
  pos.add("input-file", 1);

  po::options_description all;
  all.add(visible).add(hidden);

  po::variables_map vm;
  po::store(po::command_line_parser(argc, argv).options(all).positional(pos).run(), vm);
  po::notify(vm);

  if (vm.count("help") || !vm.count("input-file")) {
    std::cout << std::endl;
    std::cout << "Compute the oscillating Rips zigzag filtration based on a point cloud, with Euclidean metric.\n\n";
    std::cout << "Usage: " << argv[0] << " [options] input-file" << std::endl << std::endl;
    std::cout << visible << std::endl;
    std::abort();
  }
}
