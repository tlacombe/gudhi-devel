/*
 *  Euclidean_geometry.h
 *  Gudhi
 *
 *  Created by Clément Maria on 1/8/14.
 *  Copyright 2014 INRIA. All rights reserved.
 *
 */

#ifndef GUDHI_EUCLIDEAN_GEOMETRY_H
#define GUDHI_EUCLIDEAN_GEOMETRY_H

#include "boost/container/container_fwd.hpp"
#include "boost/iterator/counting_iterator.hpp"
#include "boost/range/counting_range.hpp"

/**
* \brief Represents the space \f$\mathbb{R}^d\f$ with
* Euclidean distance.
*
* The Vertices are contiguous integers from 0 to n-1.
*
* \implements MetricSpace
*
* \todo template with an integer to represent any
* Lp distance.
*/
template < class Point >    //way to encode a point -> is a range of coordinates
class Euclidean_geometry {
public:
  typedef int                                            Vertex;
/** Distance value type.*/
  typedef double                                         FT;
  typedef boost::counting_iterator< Vertex >             Space_vertex_iterator;
  typedef boost::iterator_range< Space_vertex_iterator > Space_vertex_range;
 

 /** \brief Returns a range over all Vertices.*/
  Space_vertex_range space_vertex_range()
  { return Space_vertex_range( Space_vertex_iterator(0),
                               Space_vertex_iterator(point_set_.size()) ); }
  
  FT distance(Vertex u, Vertex v) 
  { return Euclidean_distance(vertex_to_point(u),
                              vertex_to_point(v)); }

  Euclidean_geometry() :
  point_set_() {}

  /** \brief A Point in a Euclidean space represented
  * by a vertor of its coordinates.*/
//  typedef std::vector< double > Point;
/**
  * \brief Initializes the space with a set of Points.
  *
  * Point_range must be a range for which the iterators
  * have 'value_type' Point.
  */
  template< class Point_range >
  void init(Point_range &points)
  { for(auto it = points.begin(); it != points.end(); ++it)
    {point_set_.push_back(*it);} }

  /** Returns the number of elements in the metric space.*/
  size_t nb_elements() { return point_set_.size(); }
  /** Returns the dimension of the Euclidean space.*/
  int dimension()
  { if(point_set_.empty()) return -1;
    return point_set_.begin().size();}
 /**
 * Returns a Vertex that is differetn from all other Vertices in the
 * space.
 */
  Vertex null_vertex() {return -1;}
  /**
  * Return true iff d(u,v) < dist.
  */
  bool closer_than(Vertex u, Vertex v, FT dist)
  { return squared_distance(u,v) <= dist*dist; }


private:
   /**
   * Returns the point corresponding to a given Vertex.
   * This allows the correspondance between combinatorial
   * vertices and geometric points.
   */
  Point & vertex_to_point(Vertex u)
  { return point_set_[u]; }
 /** 
  * \brief Euclidean distance.
  */
  FT Euclidean_distance ( Point &p1, Point &p2 )
  {  return std::sqrt(squared_distance(p1,p2)); }
  /**
  * \brief Euclidean distance squared.
  *
  * Faster to compute than Euclidean distance and
  * enough for comparison.
  */
  FT squared_distance ( Vertex u, Vertex v)
  { return squared_distance( vertex_to_point(u),
                             vertex_to_point(v));}
  FT squared_distance ( Point &p1, Point &p2 )
  { double squared_distance = 0.;
    typename Point::iterator p1_it = p1.begin();
    typename Point::iterator p2_it = p2.begin();
    for(; p1_it != p1.end(); p1_it++, p2_it++)
    {  squared_distance += (*p1_it - *p2_it) * (*p1_it - *p2_it); }
    return squared_distance; 
  }

/** 
  * \brief Point set represented as a vector of Points.
  *
  * \todo general PropertyMap Vertex -> Point.
  */
  std::vector< Point >      point_set_;   //property map Vertex -> Point

 };

#endif // GUDHI_EUCLIDEAN_GEOMETRY_H
