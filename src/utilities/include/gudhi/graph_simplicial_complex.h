


#ifndef GUDHI_GRAPH_SIMPLICIAL_COMPLEX_FILTRATION_TAG_H
#define GUDHI_GRAPH_SIMPLICIAL_COMPLEX_FILTRATION_TAG_H

/* Edge tag for Boost PropertyGraph
*/
struct edge_filtration_t {
  typedef boost::edge_property_tag kind;
};
/* Vertex tag for Boost PropertyGraph
*/
struct vertex_filtration_t {
  typedef boost::vertex_property_tag kind;
};

#endif // GUDHI_GRAPH_SIMPLICIAL_COMPLEX_FILTRATION_TAG_H
