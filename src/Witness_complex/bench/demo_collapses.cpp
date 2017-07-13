/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Siargey Kachanovich
 *
 *    Copyright (C) 2016  INRIA Sophia Antipolis-Méditerranée (France)
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

#include <gudhi/Simplex_tree.h>
#include <gudhi/Flag_simplex_tree.h>
#include <gudhi/Euclidean_strong_witness_complex.h>
#include <gudhi/Euclidean_witness_complex.h>
#include <gudhi/Dim_lists.h>
#include <gudhi/reader_utils.h>
#include <gudhi/Persistent_cohomology.h>
#include <gudhi/Good_links.h>
#include <gudhi/pick_n_random_points.h>
#include <gudhi/choose_n_farthest_points.h>
#include <gudhi/sparsify_point_set.h>
#include <gudhi/Points_off_io.h>

#include <iostream>
#include <fstream>
#include <ctime>
#include <utility>
#include <algorithm>
#include <set>
#include <queue>
#include <iterator>
#include <string>

#include <CGAL/Epick_d.h>

#include <boost/tuple/tuple.hpp>
#include <boost/iterator/zip_iterator.hpp>
#include <boost/iterator/counting_iterator.hpp>
#include <boost/range/iterator_range.hpp>

#include <boost/program_options.hpp>

#include "../example/generators.h"
#include "output.h"
#include "output_tikz.h"
#include "tetrahedron_properties.h"


using namespace Gudhi;
using namespace Gudhi::witness_complex;
using namespace Gudhi::persistent_cohomology;
using namespace Gudhi::subsampling;

typedef CGAL::Epick_d<CGAL::Dynamic_dimension_tag> K;
typedef typename K::Point_d Point_d;
typedef std::vector<Point_d> Point_Vector;
//typedef Simplex_tree<Simplex_tree_options_fast_persistence> STree;
typedef Simplex_tree<> STree;
typedef Flag_simplex_tree<> FTree;
typedef STree::Simplex_handle Simplex_handle;

typedef Euclidean_strong_witness_complex<K> SRWit;
typedef Euclidean_witness_complex<K> WRWit;

/** Program options ***************************************************************
***********************************************************************************
***********************************************************************************
***********************************************************************************
**********************************************************************************/


void program_options(int argc, char * const argv[]
                     , int & experiment_number
                     , std::string & filepoints
                     , std::string & landmark_file
                     , std::string & experiment_name
                     , int & nbL
                     , double & alpha2_s
                     , double & alpha2_w
                     , double & mu_epsilon
                     , int & dim_max
                     , std::vector<int> & desired_homology
                     , double & min_persistence) {
  namespace po = boost::program_options;
  po::options_description hidden("Hidden options");
  hidden.add_options()
    ("option", po::value<int>(& experiment_number),
     "Experiment id.")
    ("input-file", po::value<std::string>(&filepoints),
       "Name of file containing a point set. Format is one point per line:   X1 ... Xd ");
  
  po::options_description visible("Allowed options", 100);
  visible.add_options()
      ("help,h", "produce help message")
      ("output-file,o", po::value<std::string>(&experiment_name)->default_value("witness"),
       "The prefix of all the output files. Default is 'witness'")
      ("landmarks,L", po::value<int>(&nbL)->default_value(0),
       "Number of landmarks.")
      ( "landmark-file,l", po::value<std::string>(&landmark_file),
        "Name of a file containing landmarks")
      ("alpha2_s,A", po::value<double>(&alpha2_s)->default_value(0),
       "Relaxation parameter for the strong complex.")
      ("alpha2_w,a", po::value<double>(&alpha2_w)->default_value(0),
       "Relaxation parameter for the weak complex.")
      ("mu_epsilon,e", po::value<double>(&mu_epsilon)->default_value(0),
       "Sparsification parameter.")
      ("cpx-dimension,d", po::value<int>(&dim_max)->default_value(1),
       "Maximal dimension of the Witness complex we want to compute.")
    
    ("homology,H", po::value<std::vector<int>>(&desired_homology)->multitoken(),
    "The desired Betti numbers.")
      ("min-persistence,m", po::value<Filtration_value>(&min_persistence),
       "Minimal lifetime of homology feature to be recorded. Default is 0. Enter a negative value to see zero length intervals");

  po::positional_options_description pos;
  pos.add("option", 1);
  pos.add("input-file", 2);

  po::options_description all;
  all.add(visible).add(hidden);

  po::variables_map vm;
  po::store(po::command_line_parser(argc, argv).
            options(all).positional(pos).run(), vm);
  po::notify(vm);

  if (vm.count("help") || !vm.count("input-file")) {
    std::cout << std::endl;
    std::cout << "Compute the persistent homology with coefficient field Z/3Z \n";
    std::cout << "of a Strong relaxed witness complex defined on a set of input points.\n \n";
    std::cout << "The output diagram contains one bar per line, written with the convention: \n";
    std::cout << "   p   dim b d \n";
    std::cout << "where dim is the dimension of the homological feature,\n";
    std::cout << "b and d are respectively the birth and death of the feature and \n";
    std::cout << "p is the characteristic of the field Z/pZ used for homology coefficients." << std::endl << std::endl;

    std::cout << "Usage: " << argv[0] << " [options] input-file" << std::endl << std::endl;
    std::cout << visible << std::endl;
    std::abort();
  }
}



/********************************************************************************************
 * Length of the good interval experiment
 *******************************************************************************************/

struct Pers_endpoint {
  double alpha;
  bool start;
  int dim;
  Pers_endpoint(double alpha_, bool start_, int dim_)
    : alpha(alpha_), start(start_), dim(dim_)
  {}
};

/*
struct less_than_key {
  inline bool operator() (const MyStruct& struct1, const MyStruct& struct2) {
    return (struct1.key < struct2.key);
  }
};
*/

template <class SimplicialComplex>
double good_interval_length(const std::vector<int> & desired_homology, SimplicialComplex & simplex_tree, double alpha2, std::string out_name = "barcodes.tikz.tex", std::string out_name2 = "pers_diag.tikz.tex")
{
  int nbL = simplex_tree.num_vertices();
  int p = 3;
  persistent_cohomology::Persistent_cohomology< SimplicialComplex, Field_Zp > pcoh(simplex_tree, true);
  pcoh.init_coefficients( p ); //initilizes the coefficient field for homology
  pcoh.compute_persistent_cohomology( 0.000001 );
  std::ofstream out_stream("pers_diag.tmp");
  pcoh.output_diagram(out_stream);
  out_stream.close();
  system("echo Hallo!");
  system("sort -k 2 -t ' ' pers_diag.tmp > pers_diag.tmp2");
  std::ifstream in_stream("pers_diag.tmp2", std::ios::in);
  std::string line;
  std::vector<Pers_endpoint> pers_endpoints;
  while (getline(in_stream, line)) {
    unsigned p, dim;
    double alpha_start, alpha_end;
    std::istringstream iss(line);
    iss >> p >> dim >> alpha_start >> alpha_end;
    if (iss.fail())
      alpha_end = alpha2;
    //std::cout << p << " " << dim << " " << alpha_start << " " << alpha_end << "\n";
    //if (dim < desired_homology.size()+1)
    if (alpha_start != alpha_end) {
    //   if (alpha_end < alpha_start)
    //     alpha_end = alpha2;
      pers_endpoints.push_back(Pers_endpoint(alpha_start, true, dim));
      pers_endpoints.push_back(Pers_endpoint(alpha_end, false, dim));
      //std::cout << p << " " << dim << " " << alpha_start << " " << alpha_end << "\n";
    }
  }
  std::cout << "desired_homology.size() = " << desired_homology.size() << "\n";
  for (auto nd: desired_homology)
    std::cout << nd << std::endl;
  std::cout << "Pers_endpoints.size = " << pers_endpoints.size() << std::endl;
  in_stream.close();
  std::sort(pers_endpoints.begin(),
            pers_endpoints.end(),
            [](const Pers_endpoint & p1, const Pers_endpoint & p2){
              return p1.alpha < p2.alpha;}
            );
  write_barcodes("pers_diag.tmp2", alpha2, out_name);
  write_persistence_diagram("pers_diag.tmp2", alpha2, out_name2);
  /*
  for (auto p: pers_endpoints) {
    std::cout << p.alpha << " " << p.dim << " " << p.start << "\n";
  }
  */
  std::vector<int> current_homology(desired_homology.size(),0);
  //current_homology[0] = 1; // for the compulsary "0 0 inf" entry
  double good_start = 0, good_end = 0;
  double sum_intervals = 0;
  int num_pieces = 0;
  bool interval_in_process = (desired_homology == current_homology);
  for (auto p: pers_endpoints) {
    /*
      std::cout << p.alpha << " " << p.dim << " ";
    if (p.start)
      std::cout << "s\n";
    else
      std::cout << "e\n";
    */
    /*
    std::cout << "Treating " << p.alpha << " " << p.dim << " " << p.start
              << " [";
    for (int v: current_homology)
      std::cout << v << " ";
    std::cout << "]\n";
    */
    if (p.start)
      current_homology[p.dim]++;
    else
      current_homology[p.dim]--;
    if (interval_in_process) {
      good_end = p.alpha;
      sum_intervals += good_end - good_start;
      std::cout << "good_start = " << good_start
                << ", good_end = " << good_end << "\n";
                
      Dim_lists<SimplicialComplex> simplices(simplex_tree, nbL-1, (good_end - good_start)/2);
      //simplices.collapse();
      //simplices.output_simplices();
      interval_in_process = false;
      //break;
    }
    else if (desired_homology == current_homology) {
      interval_in_process = true;
      good_start = p.alpha;
      num_pieces++;
    }
  }
  std::cout << "Number of good homology intervals: " << num_pieces << "\n";
  return sum_intervals;
}

/////////////////////////////////////////////
// Other stuff
////////////////////////////////////////////

/**
 * \brief Customized version of read_points
 * which takes into account a possible nbP first line
 *
 */
inline void
read_points_cust(std::string file_name, std::vector< std::vector< double > > & points) {
  std::ifstream in_file(file_name.c_str(), std::ios::in);
  if (!in_file.is_open()) {
    std::cerr << "Unable to open file " << file_name << std::endl;
    return;
  }
  std::string line;
  double x;
  while (getline(in_file, line)) {
    std::vector< double > point;
    std::istringstream iss(line);
    while (iss >> x) {
      point.push_back(x);
    }
    if (point.size() != 1)
      points.push_back(point);
  }
  in_file.close();
}

void rips(Point_Vector & points, double alpha2, int dim_max, STree& st)
{
  Graph_t prox_graph = compute_proximity_graph(points, sqrt(alpha2)
                                               , euclidean_distance<std::vector<FT> >);
  // Construct the Rips complex in a Simplex Tree
  // insert the proximity graph in the simplex tree
  st.insert_graph(prox_graph);
  // expand the graph until dimension dim_max
  st.expansion(dim_max);
}

void output_experiment_information(char * const file_name)
{
    std::cout << "Enter a valid experiment number. Usage: "
              << file_name << " exp_no options\n";
    std::cout << "Experiment description:\n"
              << "0 nbP nbL dim alpha limD mu_epsilon: "
              << "Build persistence diagram on relaxed witness complex "
              << "built from a point cloud on (dim-1)-dimensional sphere "
              << "consisting of nbP witnesses and nbL landmarks. "
              << "The maximal relaxation is alpha and the limit on simplicial complex "
              << "dimension is limD.\n";
    std::cout << "1 file_name nbL alpha limD: "
              << "Build persistence diagram on relaxed witness complex "
              << "build from a point cloud stored in a file and nbL landmarks. "
              << "The maximal relaxation is alpha and the limit on simplicial complex dimension is limD\n";
}

void sparsify_until(Point_Vector& point_vector,
                    double mu_epsilon,
                    unsigned nbL,
                    Point_Vector& landmarks)
{
  Point_Vector landmarks_sparsified;
  sparsify_point_set(K(), point_vector, mu_epsilon, std::back_inserter(landmarks_sparsified)); 
  if (landmarks_sparsified.size() > nbL)
    pick_n_random_points(landmarks_sparsified, nbL, std::back_inserter(landmarks));
  else
    landmarks = landmarks_sparsified;
}

template <class Simplicial_complex >
void collapse_experiment(Simplicial_complex& simplex_tree, std::vector<Point_d>& landmarks, std::vector<int>& landmarks_ind,
                         bool is2d, bool l_is_v, int nbL, FT alpha2, int limD, std::string mesh_filename = "witness")
{
  clock_t start, end;
  Dim_lists<Simplicial_complex> simplices(simplex_tree, limD, alpha2); 
  start = clock();
  simplices.collapse();
  end = clock();
  double time = static_cast<double>(end - start) / CLOCKS_PER_SEC;
  std::cout << "Collapses took "
            << time << " s. \n";
  // simplices.output_simplices();
  Simplicial_complex collapsed_tree;
  for (auto sh: simplices) {
    std::vector<int> vertices;
    for (int v: collapsed_tree.simplex_vertex_range(sh))
      vertices.push_back(v);
    collapsed_tree.insert_simplex(vertices, simplex_tree.filtration(sh));
  }
  // for (unsigned i = 0; i != nbL; ++i) {
  //   if (distances[i][0] == 0)
  //     landmarks_ind[knn[i][0]] = i;
  // }
  write_witness_mesh(landmarks, landmarks_ind, collapsed_tree, collapsed_tree.complex_simplex_range(), is2d, l_is_v, mesh_filename+"_after_collapse.mesh");

  std::cout << "Before collapse:\n";
  std::vector<int> simplices_by_dimension(limD, 0);
  for (auto sh: simplex_tree.complex_simplex_range()) {
    while ((unsigned)simplex_tree.dimension(sh) > simplices_by_dimension.size()+1)
      simplices_by_dimension.push_back(0);
    simplices_by_dimension[simplex_tree.dimension(sh)]++;
  }
  for (unsigned i = 0; i < simplices_by_dimension.size() && simplices_by_dimension[i] != 0; i++)
    std::cout << "Dimension " << i << ": " << simplices_by_dimension[i] << " simplices.\n";

  int alt_sum = 0;
  for (unsigned i = 0; i < simplices_by_dimension.size(); i++)
    alt_sum += (1-2*(i%2))*simplices_by_dimension[i];
  std::cout << "Alt_sum = " << alt_sum << std::endl;

  std::cout << "After collapse:\n";
  int chi = 0;
  for (auto sh: simplices)
    chi += 1-2*(simplex_tree.dimension(sh)%2);
  std::cout << "Euler characteristic is " << chi << std::endl;

  simplices_by_dimension = std::vector<int>(limD, 0);
  for (auto sh: simplices) {
    while ((unsigned)simplex_tree.dimension(sh) > simplices_by_dimension.size()+1)
      simplices_by_dimension.push_back(0);
    simplices_by_dimension[simplex_tree.dimension(sh)]++;
  }
  unsigned i = 0;
  for (i = 0; i < simplices_by_dimension.size() && simplices_by_dimension[i] != 0; i++)
    std::cout << "Dimension " << i << ": " << simplices_by_dimension[i] << " simplices.\n";
  std::cout << "Total simplices before: " << simplex_tree.num_simplices() << " after: " << collapsed_tree.num_simplices() << "\n";

  double sum_angles = 0;

  if (i == 3) // If the collapsed complex is 2d
    for (auto sh: simplices)
      if (simplex_tree.dimension(sh) == 2) {
        double tetra[12];
        tetra[0] = 0; tetra[1] = 0; tetra[2] = 0;
        unsigned j = 3;
        for (auto v: simplex_tree.simplex_vertex_range(sh))
          for (auto pci = landmarks[v].begin(); pci != landmarks[v].end(); ++pci, ++j)
            tetra[j] = *pci;
        double* angles = tetrahedron_solid_angles(tetra);
        //std::cout << angles[0] << " " << angles[1] << " " << angles[2] << " " << angles[3] << "\n";
        sum_angles += angles[0];
      }
  std::cout << "Sum of angles = " << sum_angles << ", sum/4pi = " << sum_angles/(4*3.1415926536) << "\n";
  Good_links<STree> gl(collapsed_tree);
  if (gl.complex_is_pseudomanifold())
    std::cout << "Collapsed complex is a pseudomanifold.\n";
  else
    std::cout << "Collapsed complex is NOT a pseudomanifold.\n";
  bool good = true;
  for (auto v: collapsed_tree.complex_vertex_range())
    if (!gl.has_good_link(v)) {
      std::cout << "Bad link around " << v << std::endl;
      good = false;
    }
  if (good)
    std::cout << "All the links are good.\n";
  else
    std::cout << "There are bad links.\n";
  std::ofstream out_stream(mesh_filename+"_before_collapse.stree");
  out_stream << simplex_tree;
  out_stream.close();
  out_stream = std::ofstream(mesh_filename+"_after_collapse.stree");
  out_stream << collapsed_tree;
  out_stream.close();
  
}

void rw_experiment(Point_Vector & point_vector, int nbL, FT alpha2, int limD, FT mu_epsilon = 0.1,
                   std::string mesh_filename = "witness", const std::vector<int>& desired_homology = std::vector<int>())
{
  clock_t start, end;
  STree simplex_tree, flag_tree, wflag_tree;

  // Choose landmarks
  //std::vector<std::vector< int > > knn;
  //std::vector<std::vector< FT > > distances;
  start = clock();
  //Gudhi::witness_complex::landmark_choice_by_random_knn(point_vector, nbL, alpha, limD, knn, distances);

  std::vector<Point_d> landmarks;
  sparsify_until(point_vector, mu_epsilon, nbL, landmarks);
  // choose_n_farthest_points(K(), point_vector, nbL, 0, std::back_inserter(landmarks));
  end = clock();
  double time = static_cast<double>(end - start) / CLOCKS_PER_SEC;
  std::cout << "Choice of " << nbL << " landmarks took "
            << time << " s. \n";
  // Compute witness complex
  start = clock();
  SRWit swit(landmarks, point_vector);
  swit.create_complex(simplex_tree, alpha2, limD);
  end = clock();
  time = static_cast<double>(end - start) / CLOCKS_PER_SEC;
  std::cout << "Witness complex for " << nbL << " landmarks took "
            << time << " s. \n";
  std::cout << "The complex contains " << simplex_tree.num_simplices() << " simplices \n";
  // Compute flag complex from witness graph
  start = clock();
  swit.create_complex(flag_tree, alpha2, 1);
  flag_tree.expansion(limD);
  end = clock();
  time = static_cast<double>(end - start) / CLOCKS_PER_SEC;
  std::cout << "Flag complex from witness graph for " << nbL << " landmarks took "
            << time << " s. \n";
  std::cout << "The complex contains " << flag_tree.num_simplices() << " simplices \n";

  // Compute flag witness complex
  start = clock();
  WRWit wwit(landmarks, point_vector);
  wwit.create_complex(wflag_tree, alpha2, 1);
  wflag_tree.expansion(limD);
  end = clock();
  time = static_cast<double>(end - start) / CLOCKS_PER_SEC;
  std::cout << "Flag witness complex for " << nbL << " landmarks took "
            << time << " s. \n";
  std::cout << "The complex contains " << wflag_tree.num_simplices() << " simplices \n";

  
  //std::cout << simplex_tree << "\n";
  
  // Compute the persistence diagram of the complex
  simplex_tree.set_dimension(limD);
  persistent_cohomology::Persistent_cohomology< STree, Field_Zp > pcoh(simplex_tree, true);
  int p = 3;
  pcoh.init_coefficients( p ); //initilizes the coefficient field for homology
  start = clock();
  pcoh.compute_persistent_cohomology( alpha2/10 );
  end = clock();
  time = static_cast<double>(end - start) / CLOCKS_PER_SEC;
  std::cout << "Persistence diagram took "
            << time << " s. \n";
  // pcoh.output_diagram();
  good_interval_length(desired_homology, simplex_tree, alpha2, mesh_filename + "_barcode_before.tikz.tex", mesh_filename + "_persdiag_before.tikz.tex");
  int chi = 0;
  for (auto sh: simplex_tree.complex_simplex_range())
    chi += 1-2*(simplex_tree.dimension(sh)%2);
  std::cout << "Euler characteristic is " << chi << std::endl;

  good_interval_length(desired_homology, flag_tree, alpha2, mesh_filename + "_flag_barcode_before.tikz.tex", mesh_filename + "_flag_persdiag_before.tikz.tex");
  chi = 0;
  for (auto sh: flag_tree.complex_simplex_range())
    chi += 1-2*(flag_tree.dimension(sh)%2);
  std::cout << "Euler characteristic is " << chi << std::endl;

  good_interval_length(desired_homology, wflag_tree, alpha2, mesh_filename + "_weak_flag_barcode_before.tikz.tex", mesh_filename + "_weak_flag_persdiag_before.tikz.tex");
    chi = 0;
  for (auto sh: wflag_tree.complex_simplex_range())
    chi += 1-2*(wflag_tree.dimension(sh)%2);
  std::cout << "Euler characteristic is " << chi << std::endl;
  
  std::vector<int> landmarks_ind(nbL); 
  bool is2d = point_vector.begin()->size() == 2, l_is_v = false;
  write_witness_mesh(landmarks, landmarks_ind, wflag_tree, wflag_tree.complex_simplex_range(), is2d, l_is_v, mesh_filename+"_before_collapse.mesh");

  // Collapse part
  std::cout << "Starting collapses for SWit...\n";
  collapse_experiment(simplex_tree, landmarks, landmarks_ind, is2d, l_is_v, nbL, alpha2, limD, mesh_filename+"_swit");

  std::cout << "Starting collapses for FlagSWit...\n";  
  collapse_experiment(flag_tree, landmarks, landmarks_ind, is2d, l_is_v, nbL, alpha2, limD, mesh_filename+"_flagswit");

  std::cout << "Starting collapses for FlagWWit...\n";  
  collapse_experiment(wflag_tree, landmarks, landmarks_ind, is2d, l_is_v, nbL, alpha2, limD, mesh_filename+"_wflag");

  
  // collapsed_tree.set_dimension(limD);
  // Persistent_cohomology< STree, Field_Zp > pcoh2(collapsed_tree, true);
  // pcoh2.init_coefficients( p ); //initilizes the coefficient field for homology
  // pcoh2.compute_persistent_cohomology( alpha2/10 );
  // pcoh2.output_diagram();
  // good_interval_length(desired_homology, simplex_tree, alpha2, mesh_filename + "barcode_after.tikz.tex");
  
  // write_witness_mesh(landmarks, landmarks_ind, collapsed_tree, collapsed_tree.complex_simplex_range(), is2d, l_is_v, mesh_filename+"_after_collapse.mesh");
}

void rips_experiment(Point_Vector & points, double threshold, int dim_max)
{
  typedef STree ST;
  clock_t start, end;
  ST st;

  // Compute the proximity graph of the points
  start = clock();
  rips(points, threshold, dim_max, st);
  end = clock();
  
  double time = static_cast<double>(end - start) / CLOCKS_PER_SEC;
  std::cout << "Rips complex took "
            << time << " s. \n";
  std::cout << "The complex contains " << st.num_simplices() << " simplices \n";
  //std::cout << "   and has dimension " << st.dimension() << " \n";

  // Sort the simplices in the order of the filtration
  st.initialize_filtration();

  // Compute the persistence diagram of the complex
  persistent_cohomology::Persistent_cohomology<ST, Field_Zp > pcoh(st);
  // initializes the coefficient field for homology
  int p = 3;
  double min_persistence = -1; //threshold/5;
  pcoh.init_coefficients(p);
  pcoh.compute_persistent_cohomology(min_persistence);
  pcoh.output_diagram();
}


int experiment0 (int argc, char * const argv[])
{
  if (argc != 8) {
    std::cerr << "Usage: " << argv[0]
              << " 0 nbP nbL dim alpha limD mu_epsilon\n";
    return 0;
  }
  /*
    boost::filesystem::path p;
    for (; argc > 2; --argc, ++argv)
    p /= argv[1];
  */
  
  int nbP       = atoi(argv[2]);
  int nbL       = atoi(argv[3]);
  int dim       = atoi(argv[4]);
  double alpha  = atof(argv[5]);
  int limD      = atoi(argv[6]);
  double mu_epsilon = atof(argv[7]);
  std::vector<int> desired_homology = {1};
  for (int i = 1; i < dim-1; ++i)
    desired_homology.push_back(0);
  desired_homology.push_back(1);
  
  // Read the point file
  Point_Vector point_vector;
  generate_points_sphere(point_vector, nbP, dim);
  std::cout << "Successfully generated " << point_vector.size() << " points.\n";
  std::cout << "Ambient dimension is " << point_vector[0].size() << ".\n";

  std::stringstream ss;
  rw_experiment(point_vector, nbL, alpha, limD, mu_epsilon, "sphere_w"+std::to_string(nbP)+"_l"+std::to_string(nbL), desired_homology);
  return 0;
}


void run_comparison(Point_Vector & landmarks,
                    Point_Vector & points,
                    unsigned limD,
                    double alpha2_s,
                    double alpha2_w,
                    std::vector<int>& desired_homology)
{
  clock_t start, end;
  STree simplex_tree;

  //std::cout << "alpha2 = " << alpha2_s << "\n";
  start = clock();
  SRWit srwit(landmarks, points);
  srwit.create_complex(simplex_tree, alpha2_s, limD);
  end = clock();
  std::cout << "SRWit.size = " << simplex_tree.num_simplices() << std::endl;
  simplex_tree.set_dimension(desired_homology.size());

  std::cout << "Good homology interval length for SRWit is "
                << good_interval_length(desired_homology, simplex_tree, alpha2_s) << "\n";
  std::cout << "Time: " << static_cast<double>(end - start) / CLOCKS_PER_SEC << " s. \n";
  int chi = 0;
  for (auto sh: simplex_tree.complex_simplex_range())
    chi += 1-2*(simplex_tree.dimension(sh)%2);
  std::cout << "Euler characteristic is " << chi << std::endl;


  STree simplex_tree2;
  std::cout << "alpha2 = " << alpha2_w << "\n";
  start = clock();
  WRWit wrwit(landmarks, points);
  wrwit.create_complex(simplex_tree, alpha2_s, limD);
  end = clock();
  std::cout << "WRWit.size = " << simplex_tree2.num_simplices() << std::endl;
  simplex_tree2.set_dimension(landmarks.size()-1);
  
  std::cout << "Good homology interval length for WRWit is "
            << good_interval_length(desired_homology, simplex_tree2, alpha2_w) << "\n";
  std::cout << "Time: " << static_cast<double>(end - start) / CLOCKS_PER_SEC << " s. \n";
  chi = 0;
  for (auto sh: simplex_tree2.complex_simplex_range())
    chi += 1-2*(simplex_tree2.dimension(sh)%2);
  std::cout << "Euler characteristic is " << chi << std::endl;

  //write_witness_mesh(points, landmarks_ind, simplex_tree2, simplex_tree2.complex_simplex_range(), false, true, "wrwit.mesh");

  
}

int experiment1 (int argc, char * const argv[])
{
  /*
    boost::filesystem::path p;
    for (; argc > 2; --argc, ++argv)
    p /= argv[1];
  */
  
  // std::string file_name = argv[2];
  // int nbL = atoi(argv[3]), limD = atoi(argv[6]);
  // double alpha2 = atof(argv[4]), mu_epsilon = atof(argv[5]);
  // std::string experiment_name = argv[7];

  int option = 1;
  std::string file_name, landmark_file;
  int nbL = 0, limD;
  double alpha2_s, alpha2_w, mu_epsilon, min_pers;
  std::string experiment_name;
  std::vector<int> desired_homology = {1};
  std::vector<Point_d> landmarks;
  
  program_options(argc, argv, option, file_name, landmark_file, experiment_name, nbL, alpha2_s, alpha2_w, mu_epsilon, limD, desired_homology, min_pers);

  // Read the point file
  Point_Vector point_vector;
  // read_points_cust(file_name, point_vector);
  Gudhi::Points_off_reader<Point_d> off_reader(file_name);
  point_vector = std::vector<Point_d>(off_reader.get_point_cloud());

  //std::cout << "The file contains " << point_vector.size() << " points.\n";
  //std::cout << "Ambient dimension is " << point_vector[0].size() << ".\n";
  //std::cout << "Limit dimension for the complexes is " << limD << ".\n";

  // if (landmark_file == "")
  //   sparsify_until(point_vector,  mu_epsilon, nbL, landmarks);
  // else
  //   read_points_cust(landmark_file, landmarks);
  // nbL = landmarks.size();
  STree simplex_tree;
  // std::vector<std::vector< int > > knn;
  // std::vector<std::vector< FT > > distances;

  rw_experiment(point_vector, nbL, alpha2_s, limD, mu_epsilon, experiment_name, desired_homology);
  //run_comparison(landmarks, point_vector, limD, alpha2_s, alpha2_w, desired_homology);
  return 0;
}


int experiment2(int argc, char * const argv[])
{
  for (unsigned d = 3; d < 4; d++) {
    // Sphere S^d 
    Point_Vector point_vector;
    unsigned N = 1;
    double alpha2 = 2.4 - 0.4*d; 
    switch (d) {
    case 1: alpha2 = 2.2; break;
    case 2: alpha2 = 1.7; break;
    case 3: alpha2 = 1.5; break;
    case 4: alpha2 = 1.4; break;
    default: alpha2 = 1.4; break;
    }
    unsigned nbL = 20;
    std::vector<int> desired_homology(nbL-1,0);
    desired_homology[0] = 1; desired_homology[d] = 1;
    
    
    for (unsigned i = 1; i <= N; ++i) {
      unsigned nbW = 1000*i;//, nbL = 20;
      double mu_epsilon = 1/sqrt(nbL);
      std::cout << "Running test S"<< d <<", |W|=" << nbW << ", |L|=" << nbL << std::endl;
      generate_points_sphere(point_vector, i*1000, d+1);
      std::vector<Point_d> landmarks;

      sparsify_until(point_vector, mu_epsilon, nbL, landmarks);
            
      std::cout << "|L| after sparsification: " << landmarks.size() << "\n";
      
      run_comparison(landmarks, point_vector, nbL-1, alpha2, alpha2, desired_homology);
    }
  }
  /*
  {
    // SO(3)
    Point_Vector point_vector;
    double alpha2 = 0.6; 
    std::cout << "alpha2 = " << alpha2 << "\n";
    unsigned nbL = 150;
    std::vector<int> desired_homology(nbL-1,0);
    desired_homology[0] = 1; desired_homology[1] = 1; desired_homology[2] = 1; //Kl
    // desired_homology[0] = 1; desired_homology[3] = 1; //SO3

    double mu_epsilon = 1/sqrt(nbL);
    if (argc < 3) std::cerr << "No file name indicated!\n";
    read_points_cust(argv[2], point_vector);
    int nbW = point_vector.size();
    std::cout << "Running test SO(3), |W|=" << nbW << ", |L|=" << nbL << std::endl;
    std::vector<Point_d> landmarks;
    Gudhi::witness_complex::landmark_choice_by_sparsification(point_vector, nbL, mu_epsilon, landmarks);
    
    std::vector<std::vector< int > > knn;
    std::vector<std::vector< FT > > distances;
    
    std::cout << "|L| after sparsification: " << landmarks.size() << "\n";
    Gudhi::witness_complex::build_distance_matrix(point_vector,   // aka witnesses
                                                  landmarks,  // aka landmarks
                                                  alpha2,
                                                  nbL-1,
                                                  knn,
                                                  distances);
    run_comparison(knn, distances, point_vector, nbL, alpha2, desired_homology);
  }
  */
  return 0;
}

int experiment3(int argc, char * const argv[])
{
  // Both witnesses and landmarks are given as input
  
  
  return 0;
}

int main (int argc, char * const argv[])
{
  if (argc == 1) {
    output_experiment_information(argv[0]);
    return 1;
  }
  switch (atoi(argv[1])) {
  case 0 :
    return experiment0(argc, argv);
    break;
  case 1 :
    return experiment1(argc, argv);
    break;
  case 2 :
    return experiment2(argc, argv);
    break;
  case 3 :
    return experiment3(argc, argv);
    break;
  default :
    output_experiment_information(argv[0]);
    return 1;
  }
}
