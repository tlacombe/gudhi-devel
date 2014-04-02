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
template < class Point >       // A Point is a range of coordinates.  
class Euclidean_geometry {
public:
/** Vertex_handle type.*/
  typedef int                                            Vertex_handle;
/** \brief Distance value type.*/
  typedef double                                         FT;
 
/** Iterator over the Vertices in the space, and corresponding range.*/
  typedef boost::counting_iterator< Vertex_handle >      Space_vertex_iterator;
  typedef boost::iterator_range< Space_vertex_iterator > Space_vertex_range;

/** \brief Construct a Euclidean Space containing no points.*/
  Euclidean_geometry() 
  : point_set_() {}

/** \brief Initializes the space with a set of Points.
  *
  * Point_range must be a range for which the iterators
  * have 'value_type' Point.*/
template < class Point_range >
void init ( Point_range & points )
{ for(auto it = points.begin(); it != points.end(); ++it)
  { point_set_.push_back(*it); } 
}

/** \brief Returns the number of elements in the metric space.*/
size_t num_elements() { return point_set_.size(); }

/** \brief Returns the dimension of the Euclidean space.*/
int dimension()
{ if(point_set_.empty()) return -1;
  return point_set_.begin().size();
}

/** \brief Returns a Vertex_handle that is different from all other 
  * Vertices in the space.*/
Vertex_handle null_vertex() { return -1; }

/** \brief Returns a range over all Vertices of the metric space.
  *
  * Iterators 'value_type' must be Vertex_handle.*/
Space_vertex_range space_vertex_range()
{ return Space_vertex_range( Space_vertex_iterator(0),
                             Space_vertex_iterator(point_set_.size()) ); }

/** \brief Returns the Euclidean distance between the points
  * associated to Vertex_handle u and Vertex_handle v.*/ 
FT distance ( Vertex_handle u
            , Vertex_handle v ) 
{ return Euclidean_distance ( point(u)
                            , point(v) ); }

/** Returns true iff d(u,v) < dist. */
  bool closer_than ( Vertex_handle u, Vertex_handle v, FT dist )
  { return squared_distance(u,v) <= dist*dist; }

private:
/** Returns the point corresponding to a given Vertex_handle.
* This allows the correspondance between combinatorial
* vertices and geometric points.*/
Point & point ( Vertex_handle u )
{ return point_set_[u]; }
/** \brief Euclidean distance.*/
FT Euclidean_distance ( Point &p1, Point &p2 )
{  return std::sqrt(squared_distance(p1,p2)); }
/** \brief Euclidean distance squared.
*
* Faster to compute than Euclidean distance and
* enough for comparison.*/
FT squared_distance ( Vertex_handle u, Vertex_handle v )
{ return squared_distance ( point(u),
                            point(v));}

FT squared_distance ( Point &p1, Point &p2 )
{ 
  double squared_distance = 0.;
  typename Point::iterator p1_it = p1.begin();
  typename Point::iterator p2_it = p2.begin();
  for(; p1_it != p1.end(); p1_it++, p2_it++)
    {  squared_distance += (*p1_it - *p2_it) * (*p1_it - *p2_it); }
  return squared_distance; }

/** \brief Point set represented as a vector of Points.*/
std::vector< Point >  point_set_;   //property map Vertex_handle -> Point

};

#endif // GUDHI_EUCLIDEAN_GEOMETRY_H
