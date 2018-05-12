/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Siargey Kachanovich
 *
 *    Copyright (C) 2018  INRIA Sophia Antipolis-Méditerranée (France)
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

#ifndef COLLAPSE_H_
#define COLLAPSE_H_

#include <boost/graph/adjacency_list.hpp>

namespace Gudhi {

/* NOTE: The elementary collapse is lazy. The coface is not removed physically, but all its adjacencies are
 * removed. Hence, the collapsed coface does not feature anywhere else in the algorithm.
 */
template <class Map,
          class Inv_map,
          class Graph,
          class CollapseTraits>
void elementary_collapse(typename Map::iterator facet_it,
                         Map* faces,
                         Inv_map* inv_faces,
                         Map* cofaces,
                         Inv_map* inv_cofaces,
                         Graph& face_coface_graph,
                         typename Map::iterator ref_it,
                         CollapseTraits& collapse_traits) {
  using Cell_comparison = typename CollapseTraits::Cell_comparison; 
  using Graph_v = typename Graph::vertex_descriptor;
  
  // check if there is one coface
  Graph_v facet_v = facet_it->second.v;
  if (boost::in_degree(facet_v, face_coface_graph) != 1)
    return;

  // check if it is valid wrt traits
  Graph_v cofacet_v = boost::source(*boost::in_edges(facet_v, face_coface_graph).first,
                                    face_coface_graph);
  typename Inv_map::iterator inv_cofacet_it = inv_cofaces->find(cofacet_v);
  typename Map::iterator cofacet_it = inv_cofacet_it->second; 
  if (!collapse_traits.is_valid_collapse(facet_it, cofacet_it))
    return;
  
  typename Graph::out_edge_iterator out_edge_it, out_edge_end;
  std::tie(out_edge_it, out_edge_end) = boost::out_edges(cofacet_v, face_coface_graph);
  std::vector<Graph_v> another_facet_vertices;
  while (out_edge_it != out_edge_end)
    another_facet_vertices.push_back(boost::target(*out_edge_it++, face_coface_graph));
  boost::clear_vertex(cofacet_v, face_coface_graph);

  // // DEBUG OUTPUT FACES
  // {
  //   std::ofstream ofs("faces.txt");
  //   ofs << "Face map:\n\n";
  //   for (auto m_pair: *faces) {
  //     ofs << "Face " << *m_pair.first << "\n";
  //     ofs << "with cofaces\n";
  //     typename Graph::in_edge_iterator in_edge_it, in_edge_end;
  //     std::tie(in_edge_it, in_edge_end) = boost::in_edges(m_pair.second.v, face_coface_graph);
  //     while (in_edge_it != in_edge_end)
  //       ofs << *(inv_cofaces->find(boost::source(*in_edge_it++, face_coface_graph))->second->first);
  //     ofs << "\n-------------------------------------------\n";
  //   } 
  //   ofs.close();
  // }
  // // DEBUG OUTPUT COFACES
  // {
  //   std::ofstream ofs("cofaces.txt");
  //   ofs << "Coface map:\n\n";
  //   for (auto m_pair: *cofaces) {
  //     ofs << "Coface " << *m_pair.first << "\n";
  //     ofs << "with faces\n";
  //     typename Graph::out_edge_iterator out_edge_it, out_edge_end;
  //     std::tie(out_edge_it, out_edge_end) = boost::out_edges(m_pair.second.v, face_coface_graph);
  //     while (out_edge_it != out_edge_end)
  //       ofs << *(inv_faces->find(boost::target(*out_edge_it++, face_coface_graph))->second->first);
  //     ofs << "\n-------------------------------------------\n";
  //   }
  //   ofs.close();
  // }
  
  for (Graph_v another_facet_v : another_facet_vertices) {
    auto inv_another_facet_it = inv_faces->find(another_facet_v);
    auto another_facet_it = inv_another_facet_it->second;
    if (inv_another_facet_it != inv_faces->end() &&
        another_facet_v != facet_v &&
        Cell_comparison()(another_facet_it->first, ref_it->first))
      elementary_collapse(another_facet_it,
                          faces,
                          inv_faces,
                          cofaces,
                          inv_cofaces,
                          face_coface_graph,
                          ref_it,
                          collapse_traits);
  }

  inv_faces->erase(inv_faces->find(facet_v));
  faces->erase(facet_it);
  boost::remove_vertex(facet_v, face_coface_graph);

  inv_cofaces->erase(inv_cofacet_it);
  cofaces->erase(cofacet_it);
  boost::remove_vertex(cofacet_v, face_coface_graph);
  
  // typedef typename Simplex_with_cofaces::List_of_facets List_of_facets;
  // Map_iterator coface_it = facet_it->second.first.coface().first;
  // // std::cout << "EColapse: facet = " << facet_it->first << ", coface = " << coface_it->first << "\n";
  // typename List_of_facets::iterator col_it = facet_it->second.first.coface().second;
  // List_of_facets& facet_list = coface_it->second.first.facets();
  // auto list_it = facet_list.begin();
  // while (list_it != facet_list.end())
  //   if (list_it != col_it) {
  //       (*list_it)->second.first.remove_coface(coface_it);
  //       if (Coface_compare<Simplex_with_cofaces>()(*list_it, *col_it))
  //         if ((*list_it)->second.first.number_of_cofaces() == 1 && (*list_it)->second.second == (*list_it)->second.first.coface().first->second.second)
  //           elementary_collapse(*list_it++);
  //         else 
  //           list_it++;
  //       else
  //         list_it++;
  //     }
  //     else
  //       list_it++;
  //   facets_->erase(facet_it);
  //   cofaces_->erase(coface_it);

  // debug note
  // break 81 if facet_it->first->position == 2359
  // break 85 if cofacet_it->first->boundary[2].first->position == 2359
}


template <class InputRange,
          class OutputIterator,
          class InputTraits,
          class CollapseTraits>
void collapse(InputRange& input_range,
              OutputIterator output_it,
              InputTraits input_traits,
              CollapseTraits collapse_traits) {
  using Graph = boost::adjacency_list< boost::listS,
                                       boost::listS,
                                       boost::bidirectionalS >;
  using Graph_v = typename Graph::vertex_descriptor;
  using Cell_type = typename CollapseTraits::Cell_type;
  using Cell_comparison = typename CollapseTraits::Cell_comparison;
  using Boundary_element = typename CollapseTraits::Boundary_element;
  struct Fields {
    Graph_v v;
    double  f;
    Fields(Graph_v v_in, double f_in) : v(v_in), f(f_in) {}
  };
  using Map = std::map<Cell_type, Fields, Cell_comparison>;
  using Inv_map = std::map<Graph_v, typename Map::iterator>;  
  Graph face_coface_graph;
  Map *faces, *cofaces;
  Inv_map *inv_faces, *inv_cofaces;
  // The common source for all cofaces. Is present during the whole execution
  Graph_v meet_v = boost::add_vertex(face_coface_graph);
  typename InputRange::iterator current_it = input_range.begin();
  if (current_it == input_range.end())
    return;
  int d = input_traits.dimension(*current_it);
  cofaces = new Map();
  inv_cofaces = new Inv_map();
  while (current_it != input_range.end() && input_traits.dimension(*current_it) == d) {
    auto cofacet_it = cofaces->find(input_traits.cell(*current_it));
    if (cofacet_it == cofaces->end()) {
      Graph_v v = boost::add_vertex(face_coface_graph);
      cofacet_it = cofaces->emplace(std::make_pair(input_traits.cell(*current_it),
                                                   Fields(v,
                                                          input_traits.filtration(*current_it)))).first;
      inv_cofaces->emplace(std::make_pair(v, cofacet_it));
    }
    else
      cofacet_it->second.f = std::min(cofacet_it->second.f, input_traits.filtration(*current_it));
    current_it++;
  }
  for (int curr_dim = d; curr_dim > 0; curr_dim--) {
    faces = new Map();
    inv_faces = new Inv_map();
    for (auto cf_pair: *cofaces)
      boost::add_edge(meet_v, cf_pair.second.v, face_coface_graph);
    while (current_it != input_range.end() && input_traits.dimension(*current_it) == curr_dim-1) {
      auto facet_it = faces->find(input_traits.cell(*current_it));
      if (facet_it == faces->end()) {
        Graph_v v = boost::add_vertex(face_coface_graph);
        if (facet_it == faces->end())
          facet_it = faces->emplace(std::make_pair(input_traits.cell(*current_it),
                                                   Fields(v,
                                                          input_traits.filtration(*current_it)))).first;
        inv_faces->emplace(std::make_pair(v, facet_it));
      }
      else
        facet_it->second.f = std::min(facet_it->second.f, input_traits.filtration(*current_it));
      current_it++;
    }
    for (auto cf_it = cofaces->begin(); cf_it != cofaces->end(); ++cf_it)
      for (Boundary_element facet: collapse_traits.boundary(cf_it)) {
        auto facet_it = faces->find(collapse_traits.facet_cell(facet));
        if (facet_it == faces->end()) {
          Graph_v v = boost::add_vertex(face_coface_graph);
          facet_it = faces->emplace(std::make_pair(collapse_traits.facet_cell(facet),
                                                   Fields(v,
                                                          collapse_traits.facet_filtration(facet)))).first;
          inv_faces->emplace(std::make_pair(v, facet_it));
        }
        else
          facet_it->second.f = std::min(facet_it->second.f, collapse_traits.facet_filtration(facet));
        boost::add_edge(cf_it->second.v, facet_it->second.v, face_coface_graph);
      }

    assert(faces->size() == inv_faces->size());
    assert(cofaces->size() == inv_cofaces->size());
    
    auto facet_it = faces->begin();
    while (facet_it != faces->end()) {
      auto prev_it = facet_it++;
      elementary_collapse(prev_it,
                          faces,
                          inv_faces,
                          cofaces,
                          inv_cofaces,
                          face_coface_graph,
                          prev_it,
                          collapse_traits);
    }
    boost::clear_vertex(meet_v, face_coface_graph);
    for (auto cf_pair: *cofaces) {
      boost::clear_vertex(cf_pair.second.v, face_coface_graph);
      boost::remove_vertex(cf_pair.second.v, face_coface_graph);
      *output_it++ = std::make_pair(cf_pair.first, cf_pair.second.f);
    }
    delete cofaces;
    delete inv_cofaces;
    cofaces = faces;
    inv_cofaces = inv_faces;
  }
  for (auto cf_pair: *cofaces)
    *output_it++ = std::make_pair(cf_pair.first, cf_pair.second.f);
  delete cofaces;
}

} // namespace Gudhi
#endif
