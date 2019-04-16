/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Vincent Rouvreau
 *
 *    Copyright (C) 2019 Inria
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

#include <gudhi/Points_off_io.h>
#include <gudhi/reader_utils.h>
#include <gudhi/distance_functions.h>
#include <gudhi/graph_simplicial_complex.h>
#include <gudhi/Flag_complex_strong_collapse.h>
#include <gudhi/Simplex_tree.h>
#include <gudhi/Rips_complex.h>
#include <gudhi/Persistent_cohomology.h>

#include <boost/program_options.hpp>

#include <string>
#include <vector>
#include <limits>  // infinity
#include <algorithm>  // for std::sort

// Types definition
using Simplex_tree = Gudhi::Simplex_tree<Gudhi::Simplex_tree_options_full_featured>;
using Filtration_value = Simplex_tree::Filtration_value;
using Point = std::vector<Filtration_value>;
using Points_off_reader = Gudhi::Points_off_reader<Point>;
using Flag_complex_strong_collapse = Gudhi::strong_collapse::Flag_complex_strong_collapse;
using Filtration_set = std::vector<Filtration_value>;
using Rips_complex = Gudhi::rips_complex::Rips_complex<Filtration_value>;
using Field_Zp = Gudhi::persistent_cohomology::Field_Zp;
using Persistent_cohomology = Gudhi::persistent_cohomology::Persistent_cohomology<Simplex_tree, Field_Zp>;

void program_options(int argc, char* argv[], std::string& off_file_points, std::string& output_file,
                     bool& dry_run, Filtration_set& edge_length_set, int& dim_max, int& p,
                     Filtration_value& min_persistence);

int main(int argc, char* argv[]) {
  std::string off_file_points;
  std::string output_file;
  bool dry_run = false;
  Filtration_set edge_length_set;
  int dim_max;
  int p;
  Filtration_value min_persistence;

  program_options(argc, argv, off_file_points, output_file, dry_run, edge_length_set, dim_max, p, min_persistence);
  std::sort (edge_length_set.begin(), edge_length_set.end());

  Points_off_reader off_reader(off_file_points);
  if (!off_reader.is_valid()) {
    std::cout << "OFF file " << off_file_points << "is not valid." << std::endl;
    exit(-1);
  }

  std::cout << "Number of points=" << off_reader.get_point_cloud().size() << std::endl;

  Gudhi::Filtered_edges_vector<Simplex_tree> edge_graph =
      Gudhi::compute_edge_graph<Gudhi::Filtered_edges_vector, Simplex_tree>(
          off_reader.get_point_cloud(),
          edge_length_set.back(),
          Gudhi::Euclidean_distance());

  std::cout << "Edge graph contains " << edge_graph.size() << " edges. ";
  std::cout << "Minimal edge length is '" << edge_graph.get_filtration_min() << "'. ";
  std::cout << "Maximal edge length is '" << edge_graph.get_filtration_max() << "'." << std::endl;

  if (dry_run) {
    exit(0);
  }

  Flag_complex_strong_collapse flag_complex(off_reader.get_point_cloud().size());
  if (edge_length_set.size() == 1)
    flag_complex.initialize_exact_version(edge_graph);
  else
    flag_complex.initialize_approximate_version(edge_graph, edge_length_set);

  Rips_complex rips_complex_after_collapse(flag_complex.get_distance_matrix(), edge_graph.get_filtration_max());

  Simplex_tree simplex_tree_after_collapse;
  rips_complex_after_collapse.create_complex(simplex_tree_after_collapse, dim_max);

  std::cout << "The complex contains " << simplex_tree_after_collapse.num_simplices() << " simplices after collapse.\n"
            << "   and has dimension " << simplex_tree_after_collapse.dimension() << " \n";

  // Sort the simplices in the order of the filtration
  simplex_tree_after_collapse.initialize_filtration();

  // Compute the persistence diagram of the complex
  Persistent_cohomology pcoh(simplex_tree_after_collapse);
  // initializes the coefficient field for homology
  pcoh.init_coefficients(p);

  pcoh.compute_persistent_cohomology(min_persistence);

  // Output the diagram in output_file
  if (output_file.empty()) {
    pcoh.output_diagram();
  } else {
    std::ofstream out(output_file);
    pcoh.output_diagram(out);
    out.close();
  }

  return 0;
}

void program_options(int argc, char* argv[], std::string& off_file_points, std::string& output_file,
                     bool& dry_run, Filtration_set& edge_length_set, int& dim_max, int& p,
                     Filtration_value& min_persistence) {
  namespace po = boost::program_options;
  po::options_description visible("Allowed options", 100);

  visible.add_options()("help,h", "produce help message")(
      "input-file,i", po::value<std::string>(&off_file_points),
      "Name of an OFF file containing a point set.")(
      "output-file,o", po::value<std::string>(&output_file)->default_value(std::string()),
      "Name of file in which the distance matrix is written. Default print in std::cout")(
      "dry-run",
      "If the dry-run is set, only a summary of the edge graph will be displayed. This can help to define edge-length-set")(
      "edge-length-set,s", po::value<Filtration_set>(&edge_length_set)->multitoken()->default_value(
          Filtration_set{std::numeric_limits<Filtration_value>::infinity()}, "{inf}"),
      "Edge length set. Default computes exact version. i.e. '-s 0.2 0.3 0.4'")(
      "cpx-dimension,d", po::value<int>(&dim_max)->default_value(1),
      "Maximal dimension of the Rips complex we want to compute.")(
      "field-charac,p", po::value<int>(&p)->default_value(11),
      "Characteristic p of the coefficient field Z/pZ for computing homology.")(
      "min-persistence,m", po::value<Filtration_value>(&min_persistence),
      "Minimal lifetime of homology feature to be recorded. Default is 0. Enter a negative value to see zero length "
      "intervals");

  po::options_description all;
  all.add(visible);

  po::variables_map vm;
  po::store(po::command_line_parser(argc, argv).options(all).run(), vm);
  po::notify(vm);

  dry_run = vm.count("dry-run");

  if (vm.count("help") || !vm.count("input-file")) {
    std::cout << std::endl;
    std::cout << "Compute the strong collapsed distance matrix \n";
    std::cout << "of a set of input points." << std::endl << std::endl;

    std::cout << "Usage: " << argv[0] << " [options] -i input-off-file" << std::endl << std::endl;
    std::cout << visible << std::endl;
    exit(-1);
  }
}
