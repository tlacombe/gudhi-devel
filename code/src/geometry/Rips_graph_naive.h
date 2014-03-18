/*
 *  Rips_graph_naive.h
 *  Gudhi
 *
 *  Created by Clément Maria on 1/8/14.
 *  Copyright 2014 INRIA. All rights reserved.
 *
 */

#ifndef GUDHI_RIPS_GRAPH_NAIVE_H
#define GUDHI_RIPS_GRAPH_NAIVE_H

#include "boost/iterator/filter_iterator.hpp"
#include "boost/graph/graph_concepts.hpp"

/**
 * \brief Represents Points in a metric space with
 * nearest neighbors queries within a radius of a point.
 *
 * The computation of nearest neighbors is naive, and 
 * traverses all points in the metric space to find the one
 * within a given radius.
 * In particular, it does not store any additional data 
 * structure apart from the set of points.
 *
 * \implements boost::AdjacencyGraph
 *
 * \todo Should we merge it to Euclidean_geometry?
 */
template < class MetricSpace >
class Rips_graph_naive
{
  public:
  typedef typename MetricSpace::Vertex_handle         Vertex_handle         ;     
  typedef typename MetricSpace::Space_vertex_iterator Space_vertex_iterator ; //iterate through the
  typedef typename MetricSpace::Space_vertex_range    Space_vertex_range    ; //vertices of the space.
  typedef typename MetricSpace::FT                    FT                    ; //distance type.
  //boost::Graph
  typedef typename MetricSpace::Vertex_handle         vertex_descriptor     ;
  typedef std::pair < vertex_descriptor
                    , vertex_descriptor >             edge_descriptor       ;
  typedef boost::undirected_tag                       directed_category     ;
  typedef boost::disallow_parallel_edge_tag           edge_parallel_category;
  typedef boost::adjacency_graph_tag                  traversal_category    ;
  //boost::AdjacencyGraph
/** \brief Predicate associated to a Vertex_handle v_:
  * returns true if an input Vertex_handle u is at distance
  * at most some threshold from v_.*/
  struct is_within_threshold_distance 
  {
    is_within_threshold_distance ( MetricSpace * ms
                                 , Vertex_handle v
                                 , FT            threshold )
   : ms_(ms), v_(v), threshold_(threshold) {}

    bool operator() (const Vertex_handle u) const
    { return ms_->closer_than(u,v_,threshold_); }

    MetricSpace *  ms_        ;
    Vertex_handle  v_         ;
    FT             threshold_ ;
  };

  typedef boost::filter_iterator< is_within_threshold_distance
                                , Space_vertex_iterator >  adjacency_iterator;
  typedef boost::iterator_range< adjacency_iterator >      adjacency_range;

 
  /** \brief Returns a range over all Vertices at distance at most
  * threshold from the Vertex_handle v. */
  adjacency_range adjacent_vertices(Vertex_handle v)
  { Space_vertex_range srg = ms_->space_vertex_range();
    return adjacency_range(
               adjacency_iterator(is_within_threshold_distance(ms_,v,threshold_),
                                  srg.begin(), srg.end() 
                                  )    ,
               adjacency_iterator(is_within_threshold_distance(ms_,v,threshold_),
                                  srg.end(), srg.end() 
                                  )
                          ); }


  /** \brief Constructs the class on a MetricSpace.*/
  Rips_graph_naive( MetricSpace & ms
                  , FT threshold ) 
  : ms_(&ms)
  , threshold_(threshold) {}

  /** \brief Returns the number of nodes in the graph.*/
  size_t size_graph()
  { return ms_->num_elements(); }

  private:
  MetricSpace *   ms_;
  FT              threshold_;

};

#endif // GUDHI_RIPS_GRAPH_NAIVE_H
