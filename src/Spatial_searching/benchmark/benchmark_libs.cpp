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

#else // RELEASE
//# define PRINT_FOUND_NEIGHBORS
# define CHECK_ACTUAL_EPSILON
//# define LOOP_ON_VARIOUS_PRECISIONS
#endif

#define GUDHI_DO_NOT_TEST_KD_TREE_SEARCH
#undef GUDHI_NMSLIB_IS_AVAILABLE
#undef GUDHI_NANOFLANN_IS_AVAILABLE
//#undef GUDHI_ANN_IS_AVAILABLE
#undef GUDHI_FLANN_IS_AVAILABLE
#undef GUDHI_COVERTREE_DNCRANE_IS_AVAILABLE  // DNCrane and Manzil cannot be used at the same time
#undef GUDHI_COVERTREE_MANZIL_IS_AVAILABLE

const int ONLY_THE_FIRST_N_POINTS = 100000; // 0 = no limit

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

#ifdef GUDHI_ANN_IS_AVAILABLE
#include "functor_ANN.h"
#endif

#ifdef GUDHI_COVERTREE_MANZIL_IS_AVAILABLE
#include "functor_Cover_tree_Manzil.h"
#endif

#ifdef GUDHI_COVERTREE_DNCRANE_IS_AVAILABLE
#include "functor_Cover_tree_DNCrane.h"
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
#include <CGAL/Memory_sizer.h>
#include <CGAL/Mesh_3/Profiling_tools.h>

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
#include <tuple>
#include <iomanip>
#include <typeinfo>

#ifdef GUDHI_USE_TBB
#include <tbb/task_scheduler_init.h>
#endif

const char * const BENCHMARK_SCRIPT_FILENAME = "benchmark_script.txt";

typedef CGAL::Epick_d<CGAL::Dynamic_dimension_tag> Kernel;
typedef Kernel::FT FT;
typedef Kernel::Point_d Point;
typedef Kernel::Vector_d Vector;

// Returns the timing for the building and the queries
template <typename ANN_Functor, typename Point, typename Range_of_query_results, typename... Targs>
std::tuple<std::string, double, double, std::size_t> test__ANN_queries(
  std::vector<Point> const& points,
  std::vector<Point> const& queries,
  int k,
  double epsilon,
  std::string const& algorith_name,
  std::string const& algorith_params,
  Range_of_query_results const *ground_truth,
  Targs... functor_additional_args)
{
  WallClockTimer t;
  
  std::cerr << black_on_white << "Testing " << algorith_name << "...\n" << white;

  SET_PERFORMANCE_DATA("Num_points", points.size());
  SET_PERFORMANCE_DATA("Epsilon", epsilon);
  SET_PERFORMANCE_DATA("Algorithm", algorith_name);
  SET_PERFORMANCE_DATA("Algo_params", algorith_params);
  SET_PERFORMANCE_DATA("Type_of_test", "ANN queries");
  SET_PERFORMANCE_DATA("Num_queries", queries.size());
  SET_PERFORMANCE_DATA("K", k);

  std::size_t mem_before = CGAL::Memory_sizer().virtual_size();

  // Build the structure
  t.reset();
  ANN_Functor functor(points, epsilon, functor_additional_args...);
  double build_time = t.elapsed();

  // Perform the queries
  std::size_t checksum = 0;
  std::vector<std::vector<std::pair<std::size_t, double>>> results(
    queries.size(), std::vector<std::pair<std::size_t, double>>());
  t.reset();
  auto it_res = results.begin();
  for (Point const& q : queries)
    checksum += functor.query_k_nearest_neighbors(q, k, epsilon, &(*it_res++));

  double q_time = t.elapsed();
  if (ground_truth) {
    auto actual_eps_and_recall = compute_actual_precision(results, *ground_truth);
    if (actual_eps_and_recall.first > epsilon)
      std::cerr << red << "WARNING: Actual epsilon = " << actual_eps_and_recall.first << " > " << epsilon << white << "\n";
    else
      std::cerr << green << "OK: Actual epsilon = " << actual_eps_and_recall.first << " <= " << epsilon << white << "\n";

    if (actual_eps_and_recall.second < 1. - epsilon)
      std::cerr << red << "WARNING: Actual recall = " << actual_eps_and_recall.second << " > " << 1. - epsilon << white << "\n";
    else
      std::cerr << green << "OK: Actual recall = " << actual_eps_and_recall.second << " >= " << 1. - epsilon << white << "\n";

    SET_PERFORMANCE_DATA("Actual_eps", actual_eps_and_recall.first);
    SET_PERFORMANCE_DATA("Actual_recall", actual_eps_and_recall.second);
  }
  else {
    SET_PERFORMANCE_DATA("Actual_eps", -1.);
    SET_PERFORMANCE_DATA("Actual_recall", -1.);
  }

  std::cerr << black_on_white << "DONE" << white << " testing " << algorith_name << " ("
    << "checksum = " << checksum
    << ", build = " << build_time
    << ", queries = " << q_time
    << ").\n\n";

  double queries_per_sec = static_cast<double>(queries.size()) / q_time;
  std::size_t mem_after = CGAL::Memory_sizer().virtual_size();

  SET_PERFORMANCE_DATA("Time1_label", "Build");
  SET_PERFORMANCE_DATA("Time1", build_time);
  SET_PERFORMANCE_DATA("Time2_label", "Queries/s");
  SET_PERFORMANCE_DATA("Time2", queries_per_sec);
  SET_PERFORMANCE_DATA("Mem_MB", (mem_after - mem_before)/(1024*1024));
  SET_PERFORMANCE_DATA("Checksum", checksum);
  XML_perf_data::commit(false);

  return std::make_tuple(algorith_name, build_time, queries_per_sec, (mem_after - mem_before) / (1024 * 1024));
}

void run_tests(
  int ambient_dim,
  std::vector<Point> const& points,
  std::vector<Point> const& queries,
  int k,
  double epsilon,
  const char *input_name = "") {

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

  //===========================================================================
  // Compute ground truth for all queries
  //===========================================================================
#ifdef CHECK_ACTUAL_EPSILON
  typedef std::vector<std::pair<std::size_t, double>> Query_res;

  GUDHI_Kd_tree_search gkts(points);

  std::vector<std::vector<std::pair<std::size_t, double>>> ground_truth(
    queries.size(), std::vector<std::pair<std::size_t, double>>());

  auto ground_truth_it = ground_truth.begin();
  for (auto const& q : queries)
    gkts.query_k_nearest_neighbors(q, k, 0., &(*ground_truth_it++));

  std::vector<std::vector<std::pair<std::size_t, double>>> const* p_ground_truth = &ground_truth;

#else
  std::vector<std::vector<std::pair<std::size_t, double>>> const* p_ground_truth = NULL;
#endif


  //===========================================================================
  // Run the tests
  //===========================================================================

  std::vector<std::tuple<std::string, double, double, std::size_t>> perfs;

  //---------------------------------------------------------------------------
  // NMSLIB
  //---------------------------------------------------------------------------

#ifdef GUDHI_NMSLIB_IS_AVAILABLE
  similarity::initLibrary(LIB_LOGNONE, NULL); // No logging

# ifdef LOOP_ON_VARIOUS_PRECISIONS
  const int M_min               = 8; // 8
  const int M_max               = 32; // 32
  const int efConstruction_min  = 100;
  const int efConstruction_max  = 400;
  const int efSearch_min        = 100;
  const int efSearch_max        = 400;
# else
  const int M_min = 16;
  const int efConstruction_min = 200;
  const int efSearch_min = 128;
  const int M_max = M_min;
  const int efConstruction_max = efConstruction_min;
  const int efSearch_max = efSearch_min;
# endif

  for (int M = M_min; M <= M_max; M *= 2) {
    std::string algo_params = "M=";
    algo_params += std::to_string(M);
    algo_params += ", post = 0";

    for (int efConstruction = efConstruction_min; efConstruction <= efConstruction_max; efConstruction *= 2) {
      std::string algo_params_2 = algo_params + ", efConstruction=";
      algo_params_2 += std::to_string(efConstruction);

      for (int efSearch = efSearch_min; efSearch <= efSearch_max; efSearch *= 2) {
        std::string algo_params_3 = algo_params_2 + ", efSearch=";
        algo_params_3 += std::to_string(efSearch);

        perfs.push_back(test__ANN_queries<NMSLIB_hnsw>(
          points, queries, k, epsilon, "NMSLIB HNSW", algo_params, p_ground_truth, M, 0, efConstruction, efSearch));
      }
    }
  }
  
  //perfs.push_back(test__ANN_queries<NMSLIB_swgraph>(
  //  points, queries, k, epsilon, "NMSLIB SWgraph", "", p_ground_truth));
#endif


  //---------------------------------------------------------------------------
  // CGAL/GUDHI
  //---------------------------------------------------------------------------

#ifdef LOOP_ON_VARIOUS_PRECISIONS
  for (double epsilon = 0.; epsilon <= 1.0; epsilon += 0.1)
#endif

#ifndef GUDHI_DO_NOT_TEST_KD_TREE_SEARCH
  perfs.push_back(test__ANN_queries<GUDHI_Kd_tree_search>(
    points, queries, k, epsilon, "GUDHI Kd_tree_search", "", p_ground_truth));
#endif

  //---------------------------------------------------------------------------
  // NANOFLANN
  //---------------------------------------------------------------------------

#ifdef GUDHI_NANOFLANN_IS_AVAILABLE
# ifdef LOOP_ON_VARIOUS_PRECISIONS
  for (double epsilon = 0.; epsilon <= 1.0; epsilon += 0.1)
# endif
  perfs.push_back(test__ANN_queries<Nanoflann>(
    points, queries, k, epsilon, "Nanoflann", "", p_ground_truth));
#endif

  //---------------------------------------------------------------------------
  // ANN
  //---------------------------------------------------------------------------

#ifdef GUDHI_ANN_IS_AVAILABLE
# ifdef LOOP_ON_VARIOUS_PRECISIONS
  for (double epsilon = 0.; epsilon <= 1.0; epsilon += 0.1)
# endif
    perfs.push_back(test__ANN_queries<Ann>(
      points, queries, k, epsilon, "ANN", "", p_ground_truth));
#endif

  //---------------------------------------------------------------------------
  // FLANN
  //---------------------------------------------------------------------------

#ifdef GUDHI_FLANN_IS_AVAILABLE

# ifdef LOOP_ON_VARIOUS_PRECISIONS

  std::vector<std::vector<std::pair<std::size_t, double>>> const* p_gt_for_flann = NULL;
  std::vector<Point> const* p_gt_queries = NULL;

  for (int flann_checks = 32; flann_checks <= 8096; flann_checks *= 2)
  {
# else

  //**********************************************************
  // Compute ground truth for N random queries
  // DO NOT CONFUSE WITH THE "all queries" GROUND TRUTH

  const int GROUND_TRUTH_FOR_FLANN_NUM_QUERIES = 200;
  typedef std::vector<std::pair<std::size_t, double>> Query_res;

  GUDHI_Kd_tree_search gkts_for_flann(points);

  std::vector<Point> gt_queries;
  gt_queries.reserve(GROUND_TRUTH_FOR_FLANN_NUM_QUERIES);
  std::vector<std::vector<std::pair<std::size_t, double>>> ground_truth_for_flann(
    GROUND_TRUTH_FOR_FLANN_NUM_QUERIES, std::vector<std::pair<std::size_t, double>>());

  auto ground_truth_for_flann_it = ground_truth_for_flann.begin();
  for (int i = 0; i < GROUND_TRUTH_FOR_FLANN_NUM_QUERIES; ++i)
  {
    // Randomly draw query point (might draw the same one sometimes, but it's ok)
    int q_index = (rand() % points.size());
    gt_queries.push_back(points[q_index]);

    gkts_for_flann.query_k_nearest_neighbors(points[q_index], k, 0., &(*ground_truth_for_flann_it++));
  }
  //**********************************************************

  const int flann_checks = 32;
  auto p_gt_for_flann = &ground_truth_for_flann;
  auto p_gt_queries = &gt_queries;
# endif

  perfs.push_back(test__ANN_queries<Flann>(
    points, queries, k, epsilon, 
    "Flann - randomized kd-trees - 4 trees",
    std::string("checks=") + std::to_string(flann_checks), 
    p_ground_truth, flann::KDTreeIndexParams(4), flann_checks, p_gt_for_flann, p_gt_queries));

  perfs.push_back(test__ANN_queries<Flann>(
    points, queries, k, epsilon,
    "Flann - randomized kd-trees - 16 trees",
    std::string("checks=") + std::to_string(flann_checks),
    p_ground_truth, flann::KDTreeIndexParams(16), flann_checks, p_gt_for_flann, p_gt_queries));

  perfs.push_back(test__ANN_queries<Flann>(
    points, queries, k, epsilon,
    "Flann - hierarchical k-means",
    std::string("checks=") + std::to_string(flann_checks),
    p_ground_truth, flann::KMeansIndexParams(), flann_checks, p_gt_for_flann, p_gt_queries));

  perfs.push_back(test__ANN_queries<Flann>(
    points, queries, k, epsilon,
    "Flann - composite (4 kd-trees + k-means)",
    std::string("checks=") + std::to_string(flann_checks),
    p_ground_truth, flann::CompositeIndexParams(4), flann_checks, p_gt_for_flann, p_gt_queries));

  perfs.push_back(test__ANN_queries<Flann>(
    points, queries, k, epsilon,
    "Flann - composite (16 kd-trees + k-means)",
    std::string("checks=") + std::to_string(flann_checks),
    p_ground_truth, flann::CompositeIndexParams(16), flann_checks, p_gt_for_flann, p_gt_queries));

  perfs.push_back(test__ANN_queries<Flann>(
    points, queries, k, epsilon,
    "Flann - hierarchical clustering",
    std::string("checks=") + std::to_string(flann_checks),
    p_ground_truth, flann::HierarchicalClusteringIndexParams(), flann_checks, p_gt_for_flann, p_gt_queries));

# ifdef LOOP_ON_VARIOUS_PRECISIONS
  }
# endif

  perfs.push_back(test__ANN_queries<Flann>(
    points, queries, k, epsilon, "Flann - linear bruteforce", "", p_ground_truth, flann::LinearIndexParams()));

  perfs.push_back(test__ANN_queries<Flann>(
    points, queries, k, epsilon, "Flann - single kd-tree", "", p_ground_truth, flann::KDTreeSingleIndexParams()));

#endif

  //---------------------------------------------------------------------------
  // Cover-tree
  //---------------------------------------------------------------------------

#ifdef GUDHI_COVERTREE_DNCRANE_IS_AVAILABLE
  perfs.push_back(test__ANN_queries<Cover_tree_DNCrane>(
    points, queries, k, epsilon, "Cover-tree DNCrane", "", p_ground_truth));
#endif

  //---------------------------------------------------------------------------
  // Cover-tree parallel
  //---------------------------------------------------------------------------

#ifdef GUDHI_COVERTREE_MANZIL_IS_AVAILABLE
  perfs.push_back(test__ANN_queries<Cover_tree_Manzil>(
    points, queries, k, epsilon, "Cover-tree Manzil", "", p_ground_truth));
#endif

  //===========================================================================
  // Display info
  //===========================================================================

  std::cerr
    << "\n================================================\n"
    << "Dimension: " << ambient_dim << "\n"
    << "Number of points: " << points.size() << "\n"
    << "Number of queries: " << queries.size() << "\n"
    << "Computation times in seconds: build + queries/s + memory (MB)\n\n";
  
  int c = 0;
  for (auto perf : perfs)
  {
    //std::cerr << (c % 2 ? white : black_on_white);
    std::cerr << std::left << "  "
      << std::setw(50) << std::setfill(' ') << std::get<0>(perf)
      << std::setw(20) << std::setfill(' ') << std::get<1>(perf)
      << std::setw(20) << std::setfill(' ') << std::get<2>(perf)
      << std::setw(20) << std::setfill(' ') << std::get<3>(perf) << "\n\n";
    ++c;
  }
  std::cerr << white;
  std::cerr << "================================================\n";
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
        int num_queries;
        int k;
        int num_iteration;
        sstr >> input;
        sstr >> param1;
        sstr >> param2;
        sstr >> param3;
        sstr >> num_points;
        sstr >> ambient_dim;
        sstr >> epsilon;
        sstr >> num_queries;
        sstr >> k;
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
          } else if (input == "generate_points_in_cube_d") {
            points = Gudhi::generate_points_in_cube_d<Kernel>(num_points, ambient_dim,
                                                              std::atof(param1.c_str()));  // side length
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
            // Contains tangent space basis
            if (input.substr(input.size() - 5) == "fvecs") {
              if (ONLY_THE_FIRST_N_POINTS > 0)
                std::cerr << yellow << "***************************************\n"
                  << "WARNING: only the first " << ONLY_THE_FIRST_N_POINTS
                  << " points have been loaded.\n"
                  << "****************************************\n\n" << white;

              load_points_from_fvecs_file<Kernel>(input, std::back_inserter(points), ONLY_THE_FIRST_N_POINTS);
            }
            else {
              Gudhi::Points_off_reader<Point> off_reader(input);
              if (!off_reader.is_valid())
                std::cerr << "Unable to read file " << input << "\n";
              else
              {
                points = off_reader.get_point_cloud();
                if (ONLY_THE_FIRST_N_POINTS > 0 && points.size() > ONLY_THE_FIRST_N_POINTS)
                {
                  std::cerr << yellow << "***************************************\n"
                    << "WARNING: only the first " << ONLY_THE_FIRST_N_POINTS
                    << " points have been loaded.\n"
                    << "****************************************\n\n" << white;

                  points.resize(ONLY_THE_FIRST_N_POINTS);
                  points.shrink_to_fit();
                }
              }
            }
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
                << "****************************************\n\n";
#endif

            std::vector<Point> queries;
            queries.reserve(num_queries);
            // Randomly pick points to fill "queries" (might get some duplicates but it is ok)
            for (int i = 0; i < num_queries; ++i) {
              queries.push_back(points[rand() % points.size()]);
            }

            run_tests(ambient_dim, points, queries, k, epsilon, input.c_str());

            std::cerr << "Run #" << i++ << " done.\n";
            std::cerr << "\n---------------------------------\n";
          } else {
            std::cerr << "Run #" << i++ << ": no points loaded.\n";
          }
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
