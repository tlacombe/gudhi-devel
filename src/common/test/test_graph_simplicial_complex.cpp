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

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE "graph_simplicial_complex"
#include <boost/test/unit_test.hpp>

#include <gudhi/Unitary_tests_utils.h>
#include <gudhi/graph_simplicial_complex.h>
#include <gudhi/Points_off_io.h>
#include <gudhi/distance_functions.h>

#include <boost/graph/adjacency_list.hpp>

#include <iostream>
#include <string>
#include <vector>
#include <limits>
#include <cstdio>  // for std::size_t

using Point_d = std::vector<double>;

struct Simplicial_complex_for_proximity_graph {
  using Filtration_value = double;
  using Vertex_handle = std::size_t;
};

BOOST_AUTO_TEST_CASE( alphacomplexdoc_compute_proximity_graph )
{
  // Read the OFF file (input file name given as parameter) and triangulates points
  Gudhi::Points_off_reader<Point_d> off_reader("alphacomplexdoc.off");
  // Check the read operation was correct
  BOOST_CHECK(off_reader.is_valid());
  
  // Retrieve the triangulation
  std::vector<Point_d> point_cloud = off_reader.get_point_cloud();

  using Proximity_graph = Gudhi::Proximity_graph<Simplicial_complex_for_proximity_graph>;

  std::cout << "===== PROXIMITY GRAPH =====" << std::endl;
  Proximity_graph prox_graph = Gudhi::compute_proximity_graph<Simplicial_complex_for_proximity_graph>(
      off_reader.get_point_cloud(),
      std::numeric_limits<double>::infinity(),
      Gudhi::Euclidean_distance()
  );
  std::cout << "num_edges = " << boost::num_edges(prox_graph)
            << " - num_vertices = " << boost::num_vertices(prox_graph) << std::endl;

  // Iterates on vertices of the graph
  Simplicial_complex_for_proximity_graph::Vertex_handle vh = 0;
  typename boost::graph_traits<Proximity_graph>::vertex_iterator v_it, v_it_end;
  for (std::tie(v_it, v_it_end) = boost::vertices(prox_graph); v_it != v_it_end; ++v_it) {
    Simplicial_complex_for_proximity_graph::Filtration_value filtration =
        boost::get(Gudhi::vertex_filtration_t(), prox_graph, *v_it);
    std::cout << "vh = " << vh << " - *v_it = " << *v_it << " - filtration = " << filtration << std::endl;
    BOOST_CHECK(vh == *v_it);
    BOOST_CHECK(filtration == 0.);
    ++vh;
  }

  // Iterates on edges of the graph
  typename boost::graph_traits<Proximity_graph>::edge_iterator e_it, e_it_end;
  for (std::tie(e_it, e_it_end) = boost::edges(prox_graph); e_it != e_it_end; ++e_it) {
    Simplicial_complex_for_proximity_graph::Filtration_value filtration =
        boost::get(Gudhi::edge_filtration_t(), prox_graph, *e_it);

    auto u = source(*e_it, prox_graph);
    auto v = target(*e_it, prox_graph);
    Simplicial_complex_for_proximity_graph::Filtration_value distance =
        Gudhi::Euclidean_distance()(point_cloud[u], point_cloud[v]);

    std::cout << "u = " << u << " - v = " << v << " - filtration = " << filtration
              << " - distance = " << distance << std::endl;
    GUDHI_TEST_FLOAT_EQUALITY_CHECK(distance, filtration);
  }

  using Filtered_proximity_graph = Gudhi::Filtered_proximity_graph<Simplicial_complex_for_proximity_graph>;

  std::cout << "===== FILTERED PROXIMITY GRAPH =====" << std::endl;
  Filtered_proximity_graph filtered_prox_graph =
      Gudhi::filter_proximity_graph<Simplicial_complex_for_proximity_graph>(prox_graph, 11.);
  std::cout << "num_edges = " << boost::num_edges(filtered_prox_graph)
            << " - num_vertices = " << boost::num_vertices(filtered_prox_graph) << std::endl;

  // Iterates on edges of the filtered graph
  typename boost::graph_traits<Filtered_proximity_graph>::edge_iterator fe_it, fe_it_end;
  for (std::tie(fe_it, fe_it_end) = boost::edges(filtered_prox_graph); fe_it != fe_it_end; ++fe_it) {
    Simplicial_complex_for_proximity_graph::Filtration_value filtration =
        boost::get(Gudhi::edge_filtration_t(), filtered_prox_graph, *fe_it);

    auto u = source(*fe_it, filtered_prox_graph);
    auto v = target(*fe_it, filtered_prox_graph);
    Simplicial_complex_for_proximity_graph::Filtration_value distance =
        Gudhi::Euclidean_distance()(point_cloud[u], point_cloud[v]);

    std::cout << "u = " << u << " - v = " << v << " - filtration = " << filtration
              << " - distance = " << distance << std::endl;
    GUDHI_TEST_FLOAT_EQUALITY_CHECK(distance, filtration);
    BOOST_CHECK(filtration <= 11.);
  }
  std::cout << "num_edges = " << boost::num_edges(filtered_prox_graph)
            << " - num_vertices = " << boost::num_vertices(filtered_prox_graph) << std::endl;
}
