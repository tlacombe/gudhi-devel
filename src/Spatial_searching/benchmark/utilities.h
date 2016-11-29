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

#include <sstream>

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

  static void commit() {
    get().commit_current_element();
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
    subelements.push_back("Epsilon");
    subelements.push_back("Algorithm");
    subelements.push_back("Type_of_test");
    subelements.push_back("Time1_label");
    subelements.push_back("Time1");
    subelements.push_back("Time2_label");
    subelements.push_back("Time2");
    subelements.push_back("Time3_label");
    subelements.push_back("Time3");
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

  void commit_current_element() {
    m_xml.add_element(m_current_element);
    m_current_element.clear();
  }

  XML_exporter m_xml;
  XML_exporter::Element_with_map m_current_element;
};

template<
  typename Kernel, typename OutputIteratorPoints>
  bool load_points_from_file(
    const std::string &filename,
    OutputIteratorPoints points,
    std::size_t only_first_n_points = std::numeric_limits<std::size_t>::max()) {
  typedef typename Kernel::Point_d Point;

  std::ifstream in(filename);
  if (!in.is_open()) {
    std::cerr << "Could not open '" << filename << "'" << std::endl;
    return false;
  }

  Kernel k;
  Point p;
  int num_ppints;
  in >> num_ppints;

  std::size_t i = 0;
  while (i < only_first_n_points && in >> p) {
    *points++ = p;
    ++i;
  }

#ifdef DEBUG_TRACES
  std::cerr << "'" << filename << "' loaded." << std::endl;
#endif

  return true;
}

#endif // UTILITIES_264161_
