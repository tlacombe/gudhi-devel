/*
 *  Nearest_neighbors.h.h
 *  Gudhi
 *
 *  Created by Clément Maria on 4/4/14.
 *  Copyright 2013 INRIA. All rights reserved.
 *
 */

#ifndef GUDHI_NEAREST_NEIGHBORS_H
#define GUDHI_NEAREST_NEIGHBORS_H

#include <iostream>
#include "io.h"
#include "Nearest_neighbors/Euclidean_geometry.h"
#include "Nearest_neighbors/Rips_graph_naive.h"

typedef std::vector< double >             Point;
typedef std::vector< Point >              Point_range;
typedef Euclidean_geometry< Point >       Metric_space;
typedef Rips_graph_naive< Metric_space >  Neighbor_graph;

/** \brief Read a set of points from filepoints and 
  * write the list of edges of the rips graph with input 
  * threshold in filegraph
  */
void compute_rips_graph( std::string filepoints
                       , std::string filegraph 
                       , double threshold )
{
  Point_range points;                //read the points from the file
  read_points( filepoints, points ); //and turn them into a Point_range
  // Create a metric space from the points, with Euclidean metric:
  Metric_space ms;
  ms.init(points);
  // Create a NeighborGraph with the space:
  Neighbor_graph ng( ms, threshold );
  //Write the edges in filegraph
  std::ofstream outgraph;
  outgraph.open ( filegraph );
  if( !outgraph.is_open() ) {
    std::cerr << "Unable to open file " << filegraph << std::endl;
    return;}

  for(auto u : ms.space_vertex_range())
  {
    outgraph << 0 << " " << u << " " << 0. << std::endl;
    for(auto v : ng.adjacent_vertices(u)) 
    { outgraph << 1 << " " << u << " " << v << " " << ms.distance(u,v) << std::endl; }
  }
  outgraph.close();
}

#endif // GUDHI_NEAREST_NEIGHBORS_H
