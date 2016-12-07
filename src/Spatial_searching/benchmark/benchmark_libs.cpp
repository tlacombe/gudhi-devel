/******************************************************************************
This benchmark allows to compare different libraries.

It reads the benchmark_script.txt file (located in the same folder as this 
file) and performs neighboring operations on the point set.
An XML file is created at each run of the benchmark. 
It contains statistics about the tests performed. This XML file 
can be processed in Excel, for example.
 ******************************************************************************/

// Without TBB_USE_THREADING_TOOL Intel Inspector XE will report false positives in Intel TBB
// (http://software.intel.com/en-us/articles/compiler-settings-for-threading-error-analysis-in-intel-inspector-xe/)
#ifdef _DEBUG
# define TBB_USE_THREADING_TOOL
#endif

#ifdef _DEBUG
# define PRINT_FOUND_NEIGHBORS
#endif

#include <cstddef>

//#define INPUT_STRIDES 3 // only take one point every INPUT_STRIDES points

#include "utilities.h"
#include "functor_GUDHI_Kd_tree_search.h"

#ifdef GUDHI_NMSLIB_IS_AVAILABLE
#include <init.h>
#include "functor_NMSLIB_hnsw.h"
#include "functor_NMSLIB_swgraph.h"
#endif

#ifdef GUDHI_NANOFLANN_IS_AVAILABLE
#include "functor_NANOFLANN.h"
#endif

#ifdef GUDHI_FLANN_IS_AVAILABLE
#include "functor_FLANN.h"
#endif

#include <gudhi/console_color.h>
#include <gudhi/Points_off_io.h>
#include <gudhi/Debug_utils.h>
#include <gudhi/Clock.h>
#include <gudhi/random_point_generators.h>

#include <CGAL/IO/Triangulation_off_ostream.h>
#include <CGAL/assertions_behaviour.h>
#include <CGAL/Epick_d.h>
#include <CGAL/Random.h>

#include <boost/algorithm/string/replace.hpp>
#include <boost/algorithm/string/trim_all.hpp>
#include <boost/range/adaptor/strided.hpp>

#include <utility>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <algorithm>
#include <map>
#include <string>
#include <iomanip>

#ifdef GUDHI_USE_TBB
#include <tbb/task_scheduler_init.h>
#endif

const char * const BENCHMARK_SCRIPT_FILENAME = "benchmark_script.txt";

typedef CGAL::Epick_d<CGAL::Dynamic_dimension_tag> Kernel;
typedef Kernel::FT FT;
typedef Kernel::Point_d Point;
typedef Kernel::Vector_d Vector;

// Returns the timing for the building and the queries
// (or -1 if time_limit reached)
template <typename ANN_Functor, typename Point, typename... Targs>
std::pair<double, double> test__ANN_queries(
  std::vector<Point> const& points,
  std::vector<Point> const& queries,
  int K,
  double epsilon,
  double time_limit = 0.,
  Targs... functor_additional_args)
{
  Gudhi::Clock t;

  // Build the structure
  t.begin();
  ANN_Functor functor(points, epsilon, functor_additional_args...);
  t.end();
  double build_time = t.num_seconds();

  // Perform the queries
  t.begin();
  for (Point const& q : queries)
  {
    if (t.num_seconds() > time_limit)
      return std::make_pair(build_time, -1.);

    functor.query_k_nearest_neighbors(q, K, epsilon);
  }
  t.end();

  return std::make_pair(build_time, t.num_seconds());
}

void run_tests(
  int ambient_dim,
  std::vector<Point> &points,
  int num_queries,
  double epsilon,
  double time_limit = 0.,
  const char *input_name = "") {

  Kernel k;

  //===========================================================================
  // Init
  //===========================================================================
  Gudhi::Clock t;

  // Get input_name_stripped
  std::string input_name_stripped(input_name);
  size_t slash_index = input_name_stripped.find_last_of('/');
  if (slash_index == std::string::npos)
    slash_index = input_name_stripped.find_last_of('\\');
  if (slash_index == std::string::npos)
    slash_index = 0;
  else
    ++slash_index;
  input_name_stripped = input_name_stripped.substr(
    slash_index, input_name_stripped.find_last_of('.') - slash_index);

  SET_PERFORMANCE_DATA("Num_points", points.size());

#ifdef GUDHI_TC_USE_ANOTHER_POINT_SET_FOR_TANGENT_SPACE_ESTIM
  std::vector<Point> points_not_sparse = points;
#endif

  SET_PERFORMANCE_DATA("Num_points", points.size());

  //===========================================================================
  // Run the test
  //===========================================================================

  std::map<std::string, std::pair<double, double>> timings;

  std::vector<Point> queries(points.begin(), points.begin() + std::min(points.size(), (std::size_t)num_queries));
  timings["GUDHI"] = test__ANN_queries<GUDHI_Kd_tree_search>(
    points, queries, 10, epsilon, time_limit);

#ifdef GUDHI_NMSLIB_IS_AVAILABLE
  similarity::initLibrary(LIB_LOGNONE, NULL); // No logging
  timings["NMSLIB HNSW"] = test__ANN_queries<NMSLIB_hnsw>(
    points, queries, 10, epsilon, time_limit);
  timings["NMSLIB SWgraph"] = test__ANN_queries<NMSLIB_swgraph>(
    points, queries, 10, epsilon, time_limit);
#endif

#ifdef GUDHI_NANOFLANN_IS_AVAILABLE
  timings["Nanoflann"] = test__ANN_queries<Nanoflann>(
    points, queries, 10, epsilon, time_limit);
#endif

#ifdef GUDHI_FLANN_IS_AVAILABLE
  timings["Flann - linear bruteforce"] = test__ANN_queries<Flann>(
    points, queries, 10, epsilon, time_limit, flann::LinearIndexParams());

  timings["Flann - randomized kd-trees"] = test__ANN_queries<Flann>(
    points, queries, 10, epsilon, time_limit, flann::KDTreeIndexParams());

  timings["Flann - hierarchical k-means"] = test__ANN_queries<Flann>(
    points, queries, 10, epsilon, time_limit, flann::KMeansIndexParams());

  timings["Flann - composite (kd-trees + k-means)"] = test__ANN_queries<Flann>(
    points, queries, 10, epsilon, time_limit, flann::CompositeIndexParams());

  timings["Flann - single kd-tree"] = test__ANN_queries<Flann>(
    points, queries, 10, epsilon, time_limit, flann::KDTreeSingleIndexParams());

  timings["Flann - hierarchical clustering"] = test__ANN_queries<Flann>(
    points, queries, 10, epsilon, time_limit, flann::HierarchicalClusteringIndexParams());
#endif

  //===========================================================================
  // Display info
  //===========================================================================

  std::cerr
    << "\n================================================\n"
    << "Dimension: " << ambient_dim << "\n"
    << "Number of points: " << points.size() << "\n"
    << "Number of queries: " << queries.size() << "\n"
    << "Computation times in seconds (build + queries): \n\n";
  
  int c = 0;
  for (auto timing : timings)
  {
    //std::cerr << (c % 2 ? white : black_on_white);
    std::cerr << std::left << "  "
      << std::setw(50) << std::setfill(' ') << timing.first
      << std::setw(10) << std::setfill(' ') << timing.second.first
      << std::setw(10) << std::setfill(' ') << timing.second.second << "\n\n";
    ++c;
  }
  std::cerr << white;
  std::cerr << "================================================\n";

  //===========================================================================
  // Export info
  //===========================================================================
  /*SET_PERFORMANCE_DATA("Init_time", init_time);
  SET_PERFORMANCE_DATA("Comput_time", computation_time);
  SET_PERFORMANCE_DATA("Perturb_successful",
                       (perturb_success ? 1 : 0));
  SET_PERFORMANCE_DATA("Perturb_time", perturb_time);
  SET_PERFORMANCE_DATA("Perturb_steps", num_perturb_steps);
  SET_PERFORMANCE_DATA("Info", "");
  */
}

int main() {
  CGAL::set_error_behaviour(CGAL::ABORT);

#ifdef GUDHI_USE_TBB
#ifdef _DEBUG
  int num_threads = 1;
#else
  int num_threads = tbb::task_scheduler_init::default_num_threads() - 4;
#endif
#endif

  unsigned int seed = static_cast<unsigned int> (time(NULL));
  CGAL::default_random = CGAL::Random(seed);  // TODO(CJ): use set_default_random
  std::cerr << "Random seed = " << seed << "\n";

  std::ifstream script_file;
  script_file.open(BENCHMARK_SCRIPT_FILENAME);
  // Script?
  // Script file format: each line gives
  //    - Filename (point set) or "generate_XXX" (point set generation)
  //    - Ambient dim
  //    - Intrinsic dim
  //    - Number of iterations with these parameters
  if (script_file.is_open()) {
    int i = 1;

#ifdef GUDHI_USE_TBB
    tbb::task_scheduler_init init(num_threads > 0 ? num_threads : tbb::task_scheduler_init::automatic);
#endif

    std::cerr << "Script file '" << BENCHMARK_SCRIPT_FILENAME << "' found.\n";
    script_file.seekg(0);
    while (script_file.good()) {
      std::string line;
      std::getline(script_file, line);
      if (line.size() > 1 && line[0] != '#') {
        boost::replace_all(line, "\t", " ");
        boost::trim_all(line);
        std::cerr << "\n\n";
        std::cerr << "*****************************************\n";
        std::cerr << "******* " << line << "\n";
        std::cerr << "*****************************************\n";
        std::stringstream sstr(line);

        std::string input;
        std::string param1;
        std::string param2;
        std::string param3;
        std::size_t num_points;
        int ambient_dim;
        double epsilon;
        double time_limit;
        int num_queries;
        int num_iteration;
        sstr >> input;
        sstr >> param1;
        sstr >> param2;
        sstr >> param3;
        sstr >> num_points;
        sstr >> ambient_dim;
        sstr >> epsilon;
        sstr >> time_limit;
        sstr >> num_queries;
        sstr >> num_iteration;

        for (int j = 0; j < num_iteration; ++j) {
          std::string input_stripped = input;
          size_t slash_index = input_stripped.find_last_of('/');
          if (slash_index == std::string::npos)
            slash_index = input_stripped.find_last_of('\\');
          if (slash_index == std::string::npos)
            slash_index = 0;
          else
            ++slash_index;
          input_stripped = input_stripped.substr(
            slash_index, input_stripped.find_last_of('.') - slash_index);

          SET_PERFORMANCE_DATA("Input", input_stripped);
          SET_PERFORMANCE_DATA("Param1", param1);
          SET_PERFORMANCE_DATA("Param2", param2);
          SET_PERFORMANCE_DATA("Param3", param3);
          SET_PERFORMANCE_DATA("Ambient_dim", ambient_dim);

#ifdef GUDHI_USE_TBB
          SET_PERFORMANCE_DATA(
            "Num_threads",
            (num_threads == -1 ? tbb::task_scheduler_init::default_num_threads() : num_threads));
#else
          SET_PERFORMANCE_DATA("Num_threads", "N/A");
#endif

          std::cerr << "\nRun #" << i << "...\n";

#ifdef GUDHI_TC_PROFILING
          Gudhi::Clock t_gen;
#endif

          std::vector<Point> points;

          if (input == "generate_moment_curve") {
            points = Gudhi::generate_points_on_moment_curve<Kernel>(num_points, ambient_dim,
                                                                    std::atof(param1.c_str()), std::atof(param2.c_str()));
          } else if (input == "generate_plane") {
            points = Gudhi::generate_points_on_plane<Kernel>(num_points, std::atoi(param1.c_str()), ambient_dim);
          } else if (input == "generate_sphere_d") {
            points = Gudhi::generate_points_on_sphere_d<Kernel>(num_points, ambient_dim,
                                                                std::atof(param1.c_str()),  // radius
                                                                std::atof(param2.c_str()));  // radius_noise_percentage
          } else if (input == "generate_two_spheres_d") {
            points = Gudhi::generate_points_on_two_spheres_d<Kernel>(num_points, ambient_dim,
                                                                     std::atof(param1.c_str()),
                                                                     std::atof(param2.c_str()),
                                                                     std::atof(param3.c_str()));
          } else if (input == "generate_3sphere_and_circle_d") {
            GUDHI_CHECK(ambient_dim == 5,
                        std::logic_error("Ambient dim should be 5"));
            points = Gudhi::generate_points_on_3sphere_and_circle<Kernel>(num_points,
                                                                          std::atof(param1.c_str()));
          } else if (input == "generate_torus_3D") {
            points = Gudhi::generate_points_on_torus_3D<Kernel>(num_points,
                                                                std::atof(param1.c_str()),
                                                                std::atof(param2.c_str()),
                                                                param3 == "Y");
          } else if (input == "generate_torus_d") {
            points = Gudhi::generate_points_on_torus_d<Kernel>(num_points,
                                                               ambient_dim / 2,
                                                               param2 == "Y",  // uniform
                                                               std::atof(param2.c_str()));  // radius_noise_percentage
          } else if (input == "generate_klein_bottle_3D") {
            points = Gudhi::generate_points_on_klein_bottle_3D<Kernel>(num_points,
                                                                       std::atof(param1.c_str()), std::atof(param2.c_str()));
          } else if (input == "generate_klein_bottle_4D") {
            points = Gudhi::generate_points_on_klein_bottle_4D<Kernel>(num_points,
                                                                       std::atof(param1.c_str()), std::atof(param2.c_str()),
                                                                       std::atof(param3.c_str()));  // noise
          } else if (input == "generate_klein_bottle_variant_5D") {
            points = Gudhi::generate_points_on_klein_bottle_variant_5D<Kernel>(num_points,
                                                                               std::atof(param1.c_str()), std::atof(param2.c_str()));
          } else {
            Gudhi::Points_off_reader<Point> off_reader(input);
            if (!off_reader.is_valid())
              std::cerr << "Unable to read file " << input << "\n";
            else
              points = off_reader.get_point_cloud();
          }

#ifdef GUDHI_TC_PROFILING
          t_gen.end();
          std::cerr << "Point set generated/loaded in " << t_gen.num_seconds()
              << " seconds.\n";
#endif

          if (!points.empty()) {
#if defined(INPUT_STRIDES) && INPUT_STRIDES > 1
            auto p = points | boost::adaptors::strided(INPUT_STRIDES);
            std::vector<Point> points(p.begin(), p.end());
            std::cerr << "****************************************\n"
                << "WARNING: taking 1 point every " << INPUT_STRIDES
                << " points.\n"
                << "****************************************\n";
#endif

            run_tests(ambient_dim, points, num_queries, epsilon, time_limit, input.c_str());

            std::cerr << "TC #" << i++ << " done.\n";
            std::cerr << "\n---------------------------------\n";
          } else {
            std::cerr << "TC #" << i++ << ": no points loaded.\n";
          }

          XML_perf_data::commit();
        }
      }
    }
    script_file.seekg(0);
    script_file.clear();

    script_file.close();
  }    // Or not script?
  else {
    std::cerr << "Script file '" << BENCHMARK_SCRIPT_FILENAME << "' NOT found.\n";
  }

#ifdef _MSC_VER
  system("pause");
#endif
  return 0;
}
