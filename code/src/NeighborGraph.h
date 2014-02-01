/*
 *  NeighborsGeometryTraits.h
 *  Gudhi
 *
 *  Created by Clément Maria on 1/8/14.
 *  Copyright 2014 INRIA. All rights reserved.
 *
 */

/**
 * \brief Defines a metric space where one can iterate
 * over a set of neighbors.
 *
 * IsModel of boost::AdjacencyGraph
 */
 template <class MetricSpace >
 class NeighborGraph {

// boost::Graph
typedef typename MetricGraph::Vertex               vertex_descriptor;
typedef typename std::pair< vertex_descriptor , 
                            vertex_descriptor >    edge_descriptor;
/** 
* Either directed_tag or undirected_tag, depending if the
* metric in MetricSpace is symetric or not.
*/                            
typedef undirected_tag                             directed_category;
typedef disallow_parallel_edge_tag                 edge_parallel_category;
typedef adjacency_graph_tag                        traversal_category;
//boost::AdjacencyGraph
/**
* Iterator type to traverse the neighbors of a Vertex in the 1-skeleton
* of the simplicial complex.
*
* 'value_type' is vertex_descriptor.
*/
typedef unspecified                                adjacency_iterator;
typedef unspecified                                adjacency_range;

std::pair<adjacency_iterator , adjacency_iterator >
adjacent_vertices(vertex_descriptor v , NeighborGraph< MetricSpace > g);

adjacency_range adjacent_vertices_range(vertex_descriptor v ,
                                        NeighborGraph< MetricSpace > g);


/**
* Returns a Vertex that is associated to no element from
* MetricSpace.
*/
vertex_descriptor null_vertex();




  /// \name Neighbor relationship:
  /// @{
  /**
   * \brief An iterator of Vertices to traverse the geometric neighbors
   * of a Point.
   */
  //typedef unspecified Neighbor_vertex_iterator;
  
  /** \brief Range of Vertices for a set of neighbors.*/
  //typedef unspecified Neighbor_vertex_range;
  
  /** \brief Returns the range of neighbors of a Vertex.
   *
   * Arbitrary order.
   * This defines the set of edges in a flag complex and, by extension,
   * defines the flag complex.
   * There is no specified order in the traversal.
   */
  //typedef neighbor_vertex_range(Vertex v);
  /// @}
  
  /// \name Conversion methods Metric Space <-> Combinatorics:
  /// @{
  /**
   * \brief Returns the point corresponding to a Vertex.
   */
  //Point vertex_to_point(const Vertex v);
  
  /**
   * \brief Returns the Vertex assigned to a Point.
   *
   * \todo Necessary?
   */
  //Vertex point_to_vertex(const Point p);
  /// @}
  
};
