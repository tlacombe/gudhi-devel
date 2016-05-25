//#undef GUDHI_USE_TBB // CJTODO TEMP

// Without TBB_USE_THREADING_TOOL Intel Inspector XE will report false positives in Intel TBB
// (http://software.intel.com/en-us/articles/compiler-settings-for-threading-error-analysis-in-intel-inspector-xe/)
#ifdef _DEBUG
# define TBB_USE_THREADING_TOOL
#endif

#include <cstddef>

//#define GUDHI_TC_USE_ANOTHER_POINT_SET_FOR_TANGENT_SPACE_ESTIM
//#define TC_PROTECT_POINT_SET_DELTA  0.003
//#define JUST_BENCHMARK_SPATIAL_SEARCH // CJTODO: test
//#define CHECK_IF_ALL_SIMPLICES_ARE_IN_THE_AMBIENT_DELAUNAY
//#define TC_INPUT_STRIDES 3 // only take one point every TC_INPUT_STRIDES points
#define TC_NO_EXPORT
//#define TC_EXPORT_TO_RIB
//#define GUDHI_TC_ALVAREZ_SURFACE_WINDOW 0.95 // 1.9 - 0.95
//#define GUDHI_TC_EXPORT_SPARSIFIED_POINT_SET
//#define GUDHI_TC_EXPORT_ALL_COORDS_IN_OFF
//#define GUDHI_TC_RECOMPUTE_TANGENT_SPACE_EVERYTIME

const std::size_t ONLY_LOAD_THE_FIRST_N_POINTS = 1000000;

#include <gudhi/Tangential_complex/RIB_exporter.h>
#include <gudhi/Debug_utils.h>
#include <gudhi/Clock.h>
#include <gudhi/Tangential_complex.h>
#include <gudhi/sparsify_point_set.h>

#include <CGAL/assertions_behaviour.h>
#include <CGAL/Epick_d.h>
#include <CGAL/Random.h>

#include "../test/testing_utilities.h" // CJTODO: won't work?

#include <boost/algorithm/string/replace.hpp>
#include <boost/algorithm/string/trim_all.hpp>
#include <boost/range/adaptor/strided.hpp>

#include <cstdlib>
#include <ctime>
#include <fstream>
#include <cmath>  // for std::sqrt

#ifdef GUDHI_USE_TBB
# include <tbb/task_scheduler_init.h>
#endif
#include "XML_exporter.h"
#define GUDHI_TC_EXPORT_PERFORMANCE_DATA
#define GUDHI_TC_SET_PERFORMANCE_DATA(value_name, value) \
        XML_perf_data::set(value_name, value);

const char * const BENCHMARK_SCRIPT_FILENAME = "benchmark_script.txt";

typedef CGAL::Epick_d<CGAL::Dynamic_dimension_tag>              Kernel;
typedef Kernel::FT                                              FT;
typedef Kernel::Point_d                                         Point;
typedef Kernel::Vector_d                                        Vector;
typedef Gudhi::Tangential_complex<
  Kernel, CGAL::Dynamic_dimension_tag,
  CGAL::Parallel_tag>                                           TC;
typedef TC::Simplex                                             Simplex;
typedef TC::Simplex_set                                         Simplex_set;


#ifdef TC_PROTECT_POINT_SET_DELTA
# include <gudhi/Tangential_complex/protected_sets.h> // CJTODO TEST
#endif

#ifdef JUST_BENCHMARK_SPATIAL_SEARCH
std::ofstream spatial_search_csv_file("benchmark_spatial_search.csv");
#endif

using namespace Gudhi::Tangential_complex_;

class XML_perf_data
{
public:
  typedef Streaming_XML_exporter<std::string> XML_exporter;

  XML_perf_data(const std::string &filename)
    : m_xml(filename, "ContainerPerformance", "Perf",
            construct_subelements_names())
  {}

  virtual ~XML_perf_data()
  {
  }

  static XML_perf_data &get()
  {
    static XML_perf_data singleton(build_filename());
    return singleton;
  }

  template <typename Value_type>
  static void set(const std::string &name, Value_type value)
  {
    get().set_data(name, value);
  }

  static void commit()
  {
    get().commit_current_element();
  }

protected:
  static std::string build_filename()
  {
    std::stringstream sstr;
    sstr << "perf_logs/Performance_log_" << time(0) << ".xml";
    return sstr.str();
  }

  static std::vector<std::string> construct_subelements_names()
  {
    std::vector<std::string> subelements;
    subelements.push_back("Input");
    subelements.push_back("Param1");
    subelements.push_back("Param2");
    subelements.push_back("Param3");
    subelements.push_back("Intrinsic_dim");
    subelements.push_back("Ambient_dim");
    subelements.push_back("Num_threads");
    subelements.push_back("Sparsity");
    subelements.push_back("Max_perturb");
    subelements.push_back("Num_points_in_input");
    subelements.push_back("Num_points");
    subelements.push_back("Perturb_technique");
    subelements.push_back("Perturb_which_points");
    subelements.push_back("Initial_num_inconsistent_local_tr");
    subelements.push_back("Best_num_inconsistent_local_tr");
    subelements.push_back("Final_num_inconsistent_local_tr");
    subelements.push_back("Init_time");
    subelements.push_back("Comput_time");
    subelements.push_back("Perturb_successful");
    subelements.push_back("Perturb_time");
    subelements.push_back("Perturb_steps");
    subelements.push_back("Add_higher_dim_simpl_time");
    subelements.push_back("Result_pure_pseudomanifold");
    subelements.push_back("Result_num_wrong_dim_simplices");
    subelements.push_back("Result_num_wrong_number_of_cofaces");
    subelements.push_back("Result_num_unconnected_stars");
    subelements.push_back("Info");

    return subelements;
  }

  void set_data(const std::string &name, const std::string &value)
  {
    m_current_element[name] = value;
  }

  template <typename Value_type>
  void set_data(const std::string &name, Value_type value)
  {
    std::stringstream sstr;
    sstr << value;
    set_data(name, sstr.str());
  }

  void commit_current_element()
  {
    m_xml.add_element(m_current_element);
    m_current_element.clear();
  }

  XML_exporter m_xml;
  XML_exporter::Element_with_map m_current_element;
};

class Test_dim
{
public:
  Test_dim(
    int min_allowed_dim = 0, 
    int max_allowed_dim = std::numeric_limits<int>::max())
    : m_min_allowed_dim(min_allowed_dim), m_max_allowed_dim(max_allowed_dim)
  {}

  template <typename Simplex>
  bool operator()(Simplex const& s)
  {
    return s.size() - 1 >= m_min_allowed_dim
      && s.size() - 1 <= m_max_allowed_dim;
  }

private:
  int m_min_allowed_dim;
  int m_max_allowed_dim;
};

// color_inconsistencies: only works if p_complex = NULL
template <typename TC>
bool export_to_off(
  TC const& tc, 
  std::string const& input_name_stripped,
  std::string const& suffix,
  bool color_inconsistencies = false,
  typename TC::Simplicial_complex const* p_complex = NULL,
  Simplex_set const *p_simpl_to_color_in_red = NULL,
  Simplex_set const *p_simpl_to_color_in_green = NULL,
  Simplex_set const *p_simpl_to_color_in_blue = NULL)
{
#ifdef TC_NO_EXPORT
  return true;
#endif

#if 0
  Kernel k;
  FT center_pt[] = { -0.5, -std::sqrt(3.) / 2, -0.5, std::sqrt(3.) / 2 };
  FT proj_pt[] = { 0., 0., 0., 0.2 };
  S3_to_R3_stereographic_projection<Kernel>
    proj_functor(0.2, 
                 Point(4, &center_pt[0], &center_pt[4]),
                 k);
#else
  CGAL::Identity<Point> proj_functor;
  //Kernel k;
  //std::array<int, 3> sel = { 1, 3, 5 };
  //Orthogonal_projection<Kernel> proj_functor(sel, k);
#endif

  if (tc.intrinsic_dimension() <= 3)
  {
    std::stringstream output_filename;
    output_filename << "output/" << input_name_stripped << "_" 
      << tc.intrinsic_dimension() << "_in_R" 
      << tc.ambient_dimension() << "_"
      << tc.number_of_vertices() << "v"
      << suffix << ".off";
    std::ofstream off_stream(output_filename.str().c_str());

    if (p_complex)
    {
#ifndef TC_NO_EXPORT
      tc.export_to_off(
        *p_complex, off_stream, 
        p_simpl_to_color_in_red,
        p_simpl_to_color_in_green, 
        p_simpl_to_color_in_blue,
        proj_functor);
#endif
    }
    else
    {
      tc.export_to_off(
        off_stream, color_inconsistencies, 
        p_simpl_to_color_in_red,
        p_simpl_to_color_in_green, 
        p_simpl_to_color_in_blue,
        NULL,
        proj_functor);
    }
    return true;
  }
  return false;
}

void make_tc(std::vector<Point> &points, 
             TC::TS_container const& tangent_spaces, // can be empty
             int intrinsic_dim,
             double sparsity = 0.01,
             double max_perturb = 0.005,
             bool perturb = true, 
             bool add_high_dim_simpl = false, 
             bool collapse = false,
             double time_limit_for_perturb = 0.,
             const char *input_name = "tc")
{
  Kernel k;
  
  if (sparsity > 0. && !tangent_spaces.empty())
  {
    std::cerr << "Error: cannot sparsify point set with pre-computed normals.\n";
    return;
  }

  // CJTODO TEMP TEST
  //TC::Simplicial_complex compl;
  //{std::size_t ss[] = {0, 1, 2}; compl.add_simplex(Simplex(ss, ss + 3)); }
  //{std::size_t ss[] = {0, 2, 3}; compl.add_simplex(Simplex(ss, ss + 3)); }
  //{std::size_t ss[] = {0, 3, 4}; compl.add_simplex(Simplex(ss, ss + 3)); }
  //{std::size_t ss[] = {0, 4, 1}; compl.add_simplex(Simplex(ss, ss + 3)); }
  //{std::size_t ss[] = {0, 5, 6}; compl.add_simplex(Simplex(ss, ss + 3)); }
  //compl.is_pure_pseudomanifold(2, 7, false, 10);

  //TC::Simplicial_complex compl;
  //{std::size_t ss[] = {0, 1, 2, 5}; compl.add_simplex(Simplex(ss, ss + 4)); }
  //{std::size_t ss[] = {0, 2, 3, 5}; compl.add_simplex(Simplex(ss, ss + 4)); }
  //{std::size_t ss[] = {0, 3, 4, 5}; compl.add_simplex(Simplex(ss, ss + 4)); }
  //{std::size_t ss[] = {0, 4, 1, 5}; compl.add_simplex(Simplex(ss, ss + 4)); }
  //{std::size_t ss[] = {0, 1, 2, 6}; compl.add_simplex(Simplex(ss, ss + 4)); }
  //{std::size_t ss[] = {0, 2, 3, 6}; compl.add_simplex(Simplex(ss, ss + 4)); }
  //{std::size_t ss[] = {0, 3, 4, 6}; compl.add_simplex(Simplex(ss, ss + 4)); }
  //{std::size_t ss[] = {0, 4, 1, 6}; compl.add_simplex(Simplex(ss, ss + 4)); }
  //{std::size_t ss[] = {0, 4, 7, 8}; compl.add_simplex(Simplex(ss, ss + 4)); }
  //compl.is_pure_pseudomanifold(3, 9, false, 10);
  // /CJTODO TEMP TEST

#ifdef JUST_BENCHMARK_SPATIAL_SEARCH
  benchmark_spatial_search(points, k, spatial_search_csv_file);
  return;
#endif

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

  int ambient_dim = k.point_dimension_d_object()(*points.begin());

  GUDHI_TC_SET_PERFORMANCE_DATA("Num_points_in_input", points.size());

#ifdef GUDHI_TC_USE_ANOTHER_POINT_SET_FOR_TANGENT_SPACE_ESTIM
  std::vector<Point> points_not_sparse = points;
#endif

  //===========================================================================
  // Sparsify point set if requested
  //===========================================================================
  if (sparsity > 0.)
  {
    std::size_t num_points_before = points.size();
    std::vector<Point> sparsified_points;
    Gudhi::sparsify_point_set(k, points, sparsity*sparsity, 
      std::back_inserter(sparsified_points));
    sparsified_points.swap(points);
    std::cerr << "Number of points before/after sparsification: "
      << num_points_before << " / " << points.size() << "\n";

#ifdef GUDHI_TC_EXPORT_SPARSIFIED_POINT_SET
    std::ofstream ps_stream("output/sparsified_point_set.txt");
    export_point_set(k, points, ps_stream);
#endif
  }

#ifdef TC_PROTECT_POINT_SET_DELTA
  // CJTODO TEST    
# ifdef GUDHI_TC_PROFILING
  Wall_clock_timer t_protection;
# endif

  std::vector<Point> points2;
  std::vector<int> dummy;
  std::vector<std::vector<int> > dummy2;
  landmark_choice_protected_delaunay(
    points, points.size(), points2, dummy, TC_PROTECT_POINT_SET_DELTA, dummy2, false, true);
  points = points2;

# ifdef GUDHI_TC_PROFILING
  t_protection.end();
  std::cerr << "Point set protected in " << t_protection.num_seconds()
    << " seconds.\n";
# endif

  std::cerr << "Number of points after PROTECTION: " << points.size() << "\n";
#endif

  GUDHI_TC_SET_PERFORMANCE_DATA("Sparsity", sparsity);
  GUDHI_TC_SET_PERFORMANCE_DATA("Max_perturb", max_perturb);
  GUDHI_TC_SET_PERFORMANCE_DATA("Num_points", points.size());

  //===========================================================================
  // Compute Tangential Complex
  //===========================================================================

#ifdef GUDHI_TC_USE_ANOTHER_POINT_SET_FOR_TANGENT_SPACE_ESTIM
  TC tc(points.begin(), points.end(), sparsity, intrinsic_dim,
    points_not_sparse.begin(), points_not_sparse.end(), max_perturb, k);
#else
  TC tc(points.begin(), points.end(), sparsity, intrinsic_dim, max_perturb, k);
#endif

  // CJTODO TEMP
  tc.estimate_intrinsic_dimension();
  //return;

  if (!tangent_spaces.empty())
  {
    tc.set_tangent_planes(tangent_spaces);
  }

  t.end();
  double init_time = t.num_seconds();
  
  t.begin();
  tc.compute_tangential_complex();
  t.end();
  double computation_time = t.num_seconds();
  
#ifdef CHECK_IF_ALL_SIMPLICES_ARE_IN_THE_AMBIENT_DELAUNAY
  if (ambient_dim <= 4)
    tc.check_if_all_simplices_are_in_the_ambient_delaunay();
#endif

  //tc.check_correlation_between_inconsistencies_and_fatness();

  //===========================================================================
  // Export to OFF
  //===========================================================================

  // Create complex
  int max_dim = -1;
  TC::Simplicial_complex complex;
  Simplex_set inconsistent_simplices;
  max_dim = tc.export_TC(complex, false, 2, &inconsistent_simplices);

  t.begin();
  bool ret = export_to_off(
    tc, input_name_stripped, "_INITIAL_TC", true, 
    &complex, &inconsistent_simplices);
  t.end();
  double export_before_time = (ret ? t.num_seconds() : -1);

  unsigned int num_perturb_steps = 0;
  double perturb_time = -1;
    double export_after_perturb_time = -1.;
  Gudhi::Fix_inconsistencies_status perturb_ret = Gudhi::FIX_NOT_PERFORMED;
  if (perturb)
  {
    //=========================================================================
    // Try to fix inconsistencies by perturbing points
    //=========================================================================
    t.begin();
    std::size_t initial_num_inconsistent_local_tr;
    std::size_t best_num_inconsistent_local_tr;
    std::size_t final_num_inconsistent_local_tr;
    perturb_ret = tc.fix_inconsistencies_using_perturbation(
      num_perturb_steps, initial_num_inconsistent_local_tr,
      best_num_inconsistent_local_tr, final_num_inconsistent_local_tr,
      time_limit_for_perturb);
    t.end();
    perturb_time = t.num_seconds();

    GUDHI_TC_SET_PERFORMANCE_DATA("Initial_num_inconsistent_local_tr", 
                                 initial_num_inconsistent_local_tr);
    GUDHI_TC_SET_PERFORMANCE_DATA("Best_num_inconsistent_local_tr", 
                                 best_num_inconsistent_local_tr);
    GUDHI_TC_SET_PERFORMANCE_DATA("Final_num_inconsistent_local_tr", 
                                 final_num_inconsistent_local_tr);

    //tc.check_correlation_between_inconsistencies_and_fatness();

    // DEBUGGING: confirm that all stars were actually refreshed
    //std::cerr << yellow << "FINAL CHECK...\n" << white;
    //std::size_t num_inc = tc.number_of_inconsistent_simplices(true).second;
    //tc.refresh_tangential_complex();
    //if (CGAL::cpp11::get<1>(tc.number_of_inconsistent_simplices(true)) != num_inc)
    //  std::cerr << red << "FINAL CHECK: FAILED.\n" << white;
    //else
    //  std::cerr << green << "FINAL CHECK: PASSED.\n" << white;


    //=========================================================================
    // Export to OFF
    //=========================================================================

    // Re-build the complex
    Simplex_set inconsistent_simplices;
    max_dim = tc.export_TC(complex, false, 2, &inconsistent_simplices);

    t.begin();
    bool exported = export_to_off(
      tc, input_name_stripped, "_AFTER_FIX", true, &complex, 
      &inconsistent_simplices);
    t.end();
    double export_after_perturb_time = (exported ? t.num_seconds() : -1);

    //std::string fn = "output/inc_stars/";
    //fn += input_name_stripped;
    //tc.export_inconsistent_stars_to_OFF_files(fn);

#if !defined(TC_NO_EXPORT) && defined(TC_EXPORT_TO_RIB)
    std::ofstream rib(std::string("output/") + input_name_stripped + ".rib");
    RIB_exporter<TC::Points, TC::Simplicial_complex::Simplex_set> rib_exporter(
      tc.points(),
      complex.simplex_range(),
      rib,
      input_name_stripped + ".tif",
      false, // is_preview
      std::make_tuple(2,4,6),
      1600, 503 // resolution
      );
    rib_exporter.write_file();

    std::ofstream rib_LQ(std::string("output/") + input_name_stripped + "_LQ.rib");
    RIB_exporter<TC::Points, TC::Simplicial_complex::Simplex_set> rib_exporter_LQ(
      tc.points(),
      complex.simplex_range(),
      rib_LQ,
      input_name_stripped + "_LQ.tif",
      true, // is_preview
      std::make_tuple(0, 4, 5)
      );
    rib_exporter_LQ.write_file();
#endif
  }
  else
  {
    GUDHI_TC_SET_PERFORMANCE_DATA("Initial_num_inconsistent_local_tr", "N/A");
    GUDHI_TC_SET_PERFORMANCE_DATA("Best_num_inconsistent_local_tr", "N/A");
    GUDHI_TC_SET_PERFORMANCE_DATA("Final_num_inconsistent_local_tr", "N/A");
  }

  // CJTODO TEST
  //tc.check_and_solve_inconsistencies_by_filtering_simplices_out();

  double fix2_time = -1;
  double export_after_fix2_time = -1.;
  if (add_high_dim_simpl)
  {
    //=========================================================================
    // Try to fix inconsistencies by adding higher-dimension simplices
    //=========================================================================
    t.begin();
    // Try to solve the remaining inconstencies
    tc.check_and_solve_inconsistencies_by_adding_higher_dim_simplices();
    t.end();
    fix2_time = t.num_seconds();

    /*Simplex_set not_delaunay_simplices;
    if (ambient_dim <= 4)
    {
      tc.check_if_all_simplices_are_in_the_ambient_delaunay(
        &complex, true, &not_delaunay_simplices);
    }*/
  
    //=========================================================================
    // Export to OFF
    //=========================================================================

    // Re-build the complex
    Simplex_set inconsistent_simplices;
    max_dim = tc.export_TC(complex, false, 2, &inconsistent_simplices);

    t.begin();
    bool exported = export_to_off(
      tc, input_name_stripped, "_AFTER_FIX2", false, &complex, 
      &inconsistent_simplices);
    t.end();
    double export_after_fix2_time = (exported ? t.num_seconds() : -1);
  }
  else
  {
    Simplex_set inconsistent_simplices;
    max_dim = tc.export_TC(complex, false, 2, &inconsistent_simplices);
  }

  complex.display_stats();

  if (intrinsic_dim == 2)
    complex.euler_characteristic(true);

  // CJTODO TEMP: Export to OFF with higher-dim simplices colored
  /*Simplex_set higher_dim_simplices;
  complex.get_simplices_matching_test(
    Test_dim(intrinsic_dim + 1),
    std::inserter(higher_dim_simplices, higher_dim_simplices.begin()));
  export_to_off(
    tc, input_name_stripped, "_BEFORE_COLLAPSE", false, &complex, 
    &inconsistent_simplices, &higher_dim_simplices);*/
  
  //===========================================================================
  // Collapse
  //===========================================================================
  if (collapse)
  {
    complex.collapse(max_dim);
    complex.display_stats();
  }

  //===========================================================================
  // Is the result a pure pseudomanifold?
  //===========================================================================
  std::size_t num_wrong_dim_simplices, 
              num_wrong_number_of_cofaces, 
              num_unconnected_stars;
  Simplex_set wrong_dim_simplices;
  Simplex_set wrong_number_of_cofaces_simplices;
  Simplex_set unconnected_stars_simplices;
  bool is_pure_pseudomanifold = complex.is_pure_pseudomanifold(
    intrinsic_dim, tc.number_of_vertices(), 
#ifdef GUDHI_TC_ALVAREZ_SURFACE_WINDOW
    true, // allow borders
#else
    false, // do NOT allow borders
#endif
    false, 1,
    &num_wrong_dim_simplices, &num_wrong_number_of_cofaces, 
    &num_unconnected_stars,
    &wrong_dim_simplices, &wrong_number_of_cofaces_simplices, 
    &unconnected_stars_simplices);

  //===========================================================================
  // Export to OFF
  //===========================================================================
  
  double export_after_collapse_time = -1.;
  if (collapse)
  {
    t.begin();
    bool exported = export_to_off(
      tc, input_name_stripped, "_AFTER_COLLAPSE", false, &complex,
      &wrong_dim_simplices, &wrong_number_of_cofaces_simplices,
      &unconnected_stars_simplices);
    t.end();
    std::cerr
      << " OFF colors:\n"
      << "   * Red: wrong dim simplices\n"
      << "   * Green: wrong number of cofaces simplices\n"
      << "   * Blue: not-connected stars\n";
    double export_after_collapse_time = (exported ? t.num_seconds() : -1.);
  }

  //===========================================================================
  // Display info
  //===========================================================================

  std::cerr
    << "\n================================================\n"
    << "Number of vertices: " << tc.number_of_vertices() << "\n"
    << "Computation times (seconds): \n"
    << "  * Tangential complex: " << init_time + computation_time << "\n"
    << "    - Init + kd-tree = " << init_time << "\n"
    << "    - TC computation = " << computation_time << "\n"
    << "  * Export to OFF (before perturb): " << export_before_time << "\n"
    << "  * Fix inconsistencies 1: " << perturb_time
    <<      " (" << num_perturb_steps << " steps) ==> "
    <<      (perturb_ret == Gudhi::TC_FIXED ? "FIXED" : "NOT fixed") << "\n"
    << "  * Fix inconsistencies 2: " << fix2_time << "\n"
    << "  * Export to OFF (after perturb): " << export_after_perturb_time << "\n"
    << "  * Export to OFF (after fix2): "<< export_after_fix2_time << "\n"
    << "  * Export to OFF (after collapse): "
    <<      export_after_collapse_time << "\n"
    << "================================================\n";
  
  //===========================================================================
  // Export info
  //===========================================================================
  GUDHI_TC_SET_PERFORMANCE_DATA("Init_time", init_time);
  GUDHI_TC_SET_PERFORMANCE_DATA("Comput_time", computation_time);
  GUDHI_TC_SET_PERFORMANCE_DATA("Perturb_successful",
                                (perturb_ret == Gudhi::TC_FIXED ? 1 : 0));
  GUDHI_TC_SET_PERFORMANCE_DATA("Perturb_time", perturb_time);
  GUDHI_TC_SET_PERFORMANCE_DATA("Perturb_steps", num_perturb_steps);
  GUDHI_TC_SET_PERFORMANCE_DATA("Add_higher_dim_simpl_time", fix2_time);
  GUDHI_TC_SET_PERFORMANCE_DATA("Result_pure_pseudomanifold",
                                (is_pure_pseudomanifold ? 1 : 0));
  GUDHI_TC_SET_PERFORMANCE_DATA("Result_num_wrong_dim_simplices",
                                num_wrong_dim_simplices);
  GUDHI_TC_SET_PERFORMANCE_DATA("Result_num_wrong_number_of_cofaces", 
                                num_wrong_number_of_cofaces);
  GUDHI_TC_SET_PERFORMANCE_DATA("Result_num_unconnected_stars", 
                                num_unconnected_stars);
  GUDHI_TC_SET_PERFORMANCE_DATA("Info", "");
}

int main()
{
  CGAL::set_error_behaviour(CGAL::ABORT);

#ifdef GUDHI_USE_TBB
# ifdef _DEBUG
  int num_threads = 1;
# else
  int num_threads = 8;
# endif
#endif

  unsigned int seed = static_cast<unsigned int>(time(NULL));
  CGAL::default_random = CGAL::Random(seed); // CJTODO: use set_default_random
  std::cerr << "Random seed = " << seed << "\n";

  std::ifstream script_file;
  script_file.open(BENCHMARK_SCRIPT_FILENAME);
  // Script?
  // Script file format: each line gives
  //    - Filename (point set) or "generate_XXX" (point set generation)
  //    - Ambient dim
  //    - Intrinsic dim
  //    - Number of iterations with these parameters
  if (script_file.is_open())
  {
    int i = 1;
#ifdef GUDHI_USE_TBB
# ifdef BENCHMARK_WITH_1_TO_MAX_THREADS
    for(num_threads = 1 ;
          num_threads <= tbb::task_scheduler_init::default_num_threads() ;
          ++num_threads)
# endif
#endif
    /*for (Concurrent_mesher_config::get().num_work_items_per_batch = 5 ;
      Concurrent_mesher_config::get().num_work_items_per_batch < 100 ;
      Concurrent_mesher_config::get().num_work_items_per_batch += 5)*/
    {
#ifdef GUDHI_USE_TBB
      tbb::task_scheduler_init init(
        num_threads > 0 ? num_threads : tbb::task_scheduler_init::automatic);
#endif

      std::cerr << "Script file '" << BENCHMARK_SCRIPT_FILENAME << "' found.\n";
      script_file.seekg(0);
      while (script_file.good())
      {
        std::string line;
        std::getline(script_file, line);
        if (line.size() > 1 && line[0] != '#')
        {
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
          int intrinsic_dim;
          double sparsity;
          double max_perturb;
          char perturb, add_high_dim_simpl, collapse;
          double time_limit_for_perturb;
          int num_iteration;
          sstr >> input;
          sstr >> param1;
          sstr >> param2;
          sstr >> param3;
          sstr >> num_points;
          sstr >> ambient_dim;
          sstr >> intrinsic_dim;
          sstr >> sparsity;
          sstr >> max_perturb;
          sstr >> perturb;
          sstr >> add_high_dim_simpl;
          sstr >> collapse;
          sstr >> time_limit_for_perturb;
          sstr >> num_iteration;

          // CJTODO TEMP
          for (int ii = 0 ; ii < 1 ; ++ii)
          for (int jj = 0 ; jj < 1 ; ++jj)
          {
          std::cerr << red
            << "******************************************************\n"
            << "*** Sparsity = " << sparsity * (1 << ii) << "\n"
            << "*** Max_perturb = " << max_perturb * (1 << jj) << "\n"
            << "******************************************************\n"
            << white;

          for (int j = 0 ; j < num_iteration ; ++j)
          {
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

            GUDHI_TC_SET_PERFORMANCE_DATA("Input", input_stripped);
            GUDHI_TC_SET_PERFORMANCE_DATA("Param1", param1);
            GUDHI_TC_SET_PERFORMANCE_DATA("Param2", param2);
            GUDHI_TC_SET_PERFORMANCE_DATA("Param3", param3);
            GUDHI_TC_SET_PERFORMANCE_DATA("Ambient_dim", ambient_dim);
            GUDHI_TC_SET_PERFORMANCE_DATA("Intrinsic_dim", intrinsic_dim);
#ifdef GUDHI_TC_PERTURB_POSITION
# ifdef GUDHI_TC_PERTURB_POSITION_TANGENTIAL
            GUDHI_TC_SET_PERFORMANCE_DATA("Perturb_technique", "Tangential_translation");
# else
            GUDHI_TC_SET_PERFORMANCE_DATA("Perturb_technique", "Ambient_translation");
# endif
#elif defined(GUDHI_TC_PERTURB_WEIGHT)
            GUDHI_TC_SET_PERFORMANCE_DATA("Perturb_technique", "Weight");
#elif defined(GUDHI_TC_PERTURB_TANGENT_SPACE)
            GUDHI_TC_SET_PERFORMANCE_DATA("Perturb_technique", "Tangent_space");
#else
            GUDHI_TC_SET_PERFORMANCE_DATA("Perturb_technique", "Undefined_value");
#endif
#ifdef GUDHI_TC_PERTURB_THE_CENTER_VERTEX_ONLY
            GUDHI_TC_SET_PERFORMANCE_DATA("Perturb_which_points", "Center_vertex");
#elif defined(GUDHI_TC_PERTURB_THE_SIMPLEX_ONLY)
            GUDHI_TC_SET_PERFORMANCE_DATA("Perturb_which_points", "Simplex");
#elif defined(GUDHI_TC_PERTURB_THE_1_STAR)
            GUDHI_TC_SET_PERFORMANCE_DATA("Perturb_which_points", "1_star");
#elif defined(GUDHI_TC_PERTURB_N_CLOSEST_POINTS)
            std::stringstream sstr;
            sstr << GUDHI_TC_NUMBER_OF_PERTURBED_POINTS(intrinsic_dim) <<
              "_closest_points";
            GUDHI_TC_SET_PERFORMANCE_DATA("Perturb_which_points", sstr.str());
#else
            GUDHI_TC_SET_PERFORMANCE_DATA("Perturb_which_points", "Undefined_value");
#endif

#ifdef GUDHI_USE_TBB
            GUDHI_TC_SET_PERFORMANCE_DATA(
              "Num_threads",
              (num_threads == -1 ? tbb::task_scheduler_init::default_num_threads() : num_threads));
#else
            GUDHI_TC_SET_PERFORMANCE_DATA("Num_threads", "N/A");
#endif

            std::cerr << "\nTC #" << i << "...\n";
          
#ifdef GUDHI_TC_PROFILING
            Gudhi::Clock t_gen;
#endif

            std::vector<Point> points;
            TC::TS_container tangent_spaces;

            if (input == "generate_moment_curve")
            {
              points = generate_points_on_moment_curve<Kernel>(
                num_points, ambient_dim,
                std::atof(param1.c_str()), std::atof(param2.c_str()));
            }
            else if (input == "generate_plane")
            {
              points = generate_points_on_plane<Kernel>(
                num_points, intrinsic_dim, ambient_dim);
            }
            else if (input == "generate_sphere_d")
            {
              points = generate_points_on_sphere_d<Kernel>(
                num_points, ambient_dim,
                std::atof(param1.c_str()),  // radius
                std::atof(param2.c_str())); // radius_noise_percentage
            }
            else if (input == "generate_two_spheres_d")
            {
              points = generate_points_on_two_spheres_d<Kernel>(
                num_points, ambient_dim,
                std::atof(param1.c_str()),
                std::atof(param2.c_str()),
                std::atof(param3.c_str()));
            }
            else if (input == "generate_3sphere_and_circle_d")
            {
              GUDHI_CHECK(intrinsic_dim == 3, 
                std::logic_error("Intrinsic dim should be 3"));
              GUDHI_CHECK(ambient_dim == 5,
                std::logic_error("Ambient dim should be 5"));
              points = generate_points_on_3sphere_and_circle<Kernel>(
                num_points,
                std::atof(param1.c_str()));
            }
            else if (input == "generate_torus_3D")
            {
              points = generate_points_on_torus_3D<Kernel>(
                num_points,
                std::atof(param1.c_str()),
                std::atof(param2.c_str()),
                param3 == "Y");
            }
            else if (input == "generate_torus_d")
            {
              points = generate_points_on_torus_d<Kernel>(
                num_points, 
                intrinsic_dim,
                param1 == "Y", // uniform
                std::atof(param2.c_str())); // radius_noise_percentage
            }
            else if (input == "generate_klein_bottle_3D")
            {
              points = generate_points_on_klein_bottle_3D<Kernel>(
                num_points,
                std::atof(param1.c_str()), std::atof(param2.c_str()));
            }
            else if (input == "generate_klein_bottle_4D")
            {
              points = generate_points_on_klein_bottle_4D<Kernel>(
                num_points,
                std::atof(param1.c_str()), std::atof(param2.c_str()),
                std::atof(param3.c_str())); // noise
            }
            else if (input == "generate_klein_bottle_variant_5D")
            {
              points = generate_points_on_klein_bottle_variant_5D<Kernel>(
                num_points,
                std::atof(param1.c_str()), std::atof(param2.c_str()));
            }
            else
            {
              load_points_from_file<Kernel, typename TC::Tangent_space_basis>(
                input, std::back_inserter(points), 
                std::back_inserter(tangent_spaces), 
                ONLY_LOAD_THE_FIRST_N_POINTS);
            }

#ifdef GUDHI_TC_PROFILING
            t_gen.end();
            std::cerr << "Point set generated/loaded in " << t_gen.num_seconds()
                      << " seconds.\n";
#endif

            if (!points.empty())
            {
#if defined(TC_INPUT_STRIDES) && TC_INPUT_STRIDES > 1
              auto p = points | boost::adaptors::strided(TC_INPUT_STRIDES); // CJTODO C++11 (auto)
              std::vector<Point> points(p.begin(), p.end());
              std::cerr << "****************************************\n"
                << "WARNING: taking 1 point every " << TC_INPUT_STRIDES
                << " points.\n"
                << "****************************************\n";
#endif

              make_tc(points, tangent_spaces, intrinsic_dim, 
                sparsity * (1 << ii), max_perturb * (1 << jj), 
                perturb == 'Y', add_high_dim_simpl == 'Y',  collapse == 'Y', 
                time_limit_for_perturb, input.c_str());

              std::cerr << "TC #" << i++ << " done.\n";
              std::cerr << "\n---------------------------------\n";
            }
            else
            {
              std::cerr << "TC #" << i++ << ": no points loaded.\n";
            }

            XML_perf_data::commit();
          }
          } // CJTODO TEMP
        }
      }
      script_file.seekg(0);
      script_file.clear();
    }

    script_file.close();
  }
  // Or not script?
  else
  {
    std::cerr << "Script file '" << BENCHMARK_SCRIPT_FILENAME << "' NOT found.\n";
  }

  system("pause");
  return 0;
}
