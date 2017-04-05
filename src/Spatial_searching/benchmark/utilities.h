/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Clement Jamin
 *
 *    Copyright (C) 2016 INRIA
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

#ifndef UTILITIES_264161_
#define UTILITIES_264161_

#include "XML_exporter.h"
#define SET_PERFORMANCE_DATA(value_name, value) \
        XML_perf_data::set(value_name, value);

#include <fstream>
#include <sstream>

template<typename Pair>
struct Second_of_pair
{
  typedef typename Pair::second_type result_type;
  result_type operator()(Pair p) const { return p.second; }
};


class XML_perf_data {
public:
  typedef Streaming_XML_exporter<std::string> XML_exporter;

  XML_perf_data(const std::string &filename)
    : m_xml(filename, "ContainerPerformance", "Perf",
      construct_subelements_names()) { }

  virtual ~XML_perf_data() { }

  static XML_perf_data &get() {
    static XML_perf_data singleton(build_filename());
    return singleton;
  }

  template <typename Value_type>
  static void set(const std::string &name, Value_type value) {
    get().set_data(name, value);
  }

  static void commit(bool clear_current_element = true) {
    get().commit_current_element(clear_current_element);
  }

protected:

  static std::string build_filename() {
    std::stringstream sstr;
    sstr << "perf_logs/Performance_log_" << time(0) << ".xml";
    return sstr.str();
  }

  static std::vector<std::string> construct_subelements_names() {
    std::vector<std::string> subelements;
    subelements.push_back("Input");
    subelements.push_back("Param1");
    subelements.push_back("Param2");
    subelements.push_back("Param3");
    subelements.push_back("Ambient_dim");
    subelements.push_back("Num_threads");
    subelements.push_back("Num_points");
    subelements.push_back("Random_rotation");
    subelements.push_back("Ambient_noise");
    subelements.push_back("Epsilon");
    subelements.push_back("Algorithm");
    subelements.push_back("Algo_params");
    subelements.push_back("Type_of_test");
    subelements.push_back("Num_queries");
    subelements.push_back("K");
    subelements.push_back("Time1_label");
    subelements.push_back("Time1");
    subelements.push_back("Time2_label");
    subelements.push_back("Time2");
    subelements.push_back("Time3_label");
    subelements.push_back("Time3");
    subelements.push_back("Mem_MB");
    subelements.push_back("Actual_eps");
    subelements.push_back("Actual_recall");
    subelements.push_back("Tree_depth");
    subelements.push_back("Info");

    return subelements;
  }

  void set_data(const std::string &name, const std::string &value) {
    m_current_element[name] = value;
  }

  template <typename Value_type>
  void set_data(const std::string &name, Value_type value) {
    std::stringstream sstr;
    sstr << value;
    set_data(name, sstr.str());
  }

  void commit_current_element(bool clear_current_element = true) {
    m_xml.add_element(m_current_element);
    if (clear_current_element)
      m_current_element.clear();
  }

  XML_exporter m_xml;
  XML_exporter::Element_with_map m_current_element;
};

template<typename Kernel, typename OutputIteratorPoints>
bool load_points_from_fvecs_file(const std::string &filename, OutputIteratorPoints points, int only_the_first_n_points = -1)
{
  typedef typename Kernel::Point_d Point;

  std::ifstream in(filename, std::ios::binary);
  if (!in.is_open()) {
    std::cerr << "Could not open '" << filename << "'" << std::endl;
    return false;
  }

  Kernel k;
  unsigned long pt_dim = 0;

  in.read(reinterpret_cast<char*>(&pt_dim), 4);
  std::vector<float> current_pt;
  current_pt.reserve(pt_dim);
  for (int c = 0; !in.fail() && c != only_the_first_n_points; ++c) {

    for (int j = 0; j < pt_dim; ++j)
    {
      float coord = 0.f;
      in.read(reinterpret_cast<char*>(&coord), 4);
      current_pt.push_back(coord);
    }

    *points++ = Point(current_pt.begin(), current_pt.end());
    current_pt.clear();
    in.read(reinterpret_cast<char*>(&pt_dim), 4);
  }

#ifdef DEBUG_TRACES
  std::cerr << "'" << filename << "' loaded." << std::endl;
#endif

  return true;
}

// Range_of_query_results(2) is a range of ranges of std::pair<std::size_t, double> (index, squared distance)
// Returns a pair<epsilon, recall>
template <typename Range_of_query_results, typename Range_of_query_results2>
std::pair<double, double> compute_actual_precision(
  Range_of_query_results const& results, 
  Range_of_query_results2 const& ground_truth)
{
  if (results.size() == 0)
  {
    std::cerr << red << "WARNING: compute_actual_precision => results is empty.\n" << white;
    return std::make_pair(-1., -1.);
  }

  double worst_ratio = 1.;
  std::size_t num_correct_answers = 0;
  // For each query
  auto ground_truth_query_result = ground_truth.begin();
  for (auto const& actual_result_for_query : results)
  {
    // For each neighbor
    auto ground_truth_res = ground_truth_query_result->begin();
    for (auto const& index_and_sqdist : actual_result_for_query)
    {
      // Is this found neighbor among the true k-nearest neighbors
      for (auto r : *ground_truth_query_result) {
        if (r.first == index_and_sqdist.first) {
          ++num_correct_answers;
          break;
        }
      }

      // Skip the case where the query is a point itself => impossible to compute epsilon
      if (ground_truth_res->second < 1e-10) {
        ++ground_truth_res;
        continue;
      }
      double ratio = index_and_sqdist.second / ground_truth_res->second;
      if (ratio > worst_ratio)
        worst_ratio = ratio;
      ++ground_truth_res;
    }
    ++ground_truth_query_result;
  }
  double actual_eps = std::sqrt(worst_ratio) - 1.;
  double recall = static_cast<double>(num_correct_answers) / (results.size()*results.begin()->size());
  //std::cerr << "Actual epsilon = " << actual_eps << "\n";
  return std::make_pair(actual_eps, recall);
}

#endif // UTILITIES_264161_
