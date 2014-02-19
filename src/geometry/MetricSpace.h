/*
 *  MetricSpace.h
 *  Gudhi
 *
 *  Created by Clément Maria on 12/10/13.
 *  Copyright 2013 INRIA. All rights reserved.
 *
 */

/** \brief Concept defining a discrete metric space. 
 *
 * The metric space contains a finite set of elements defined by
 * a type Vertex, and a distance between each pair of Vertices is
 * defined.
 *
 * Remark: the distance may be more general, like a similarity 
 * measure, and simply needs to be symmetric.
 * 
 * If the space does not have a distance, we define the discrete
 * distance satisfying: \f$\text{d}(x,y) = 0\f$ if \f$x \neq y\f$
 * and \f$0\f$ otherwise. */
 struct MetricSpace
 {
/** \brief Vertex type representing an element of the metric space.
  *
  * The set of Vertices must be totally ordered with <.
  * A Vertex must be Default Constructible, Assignable, and Equality Comparable.
  */
  typedef unspecified Vertex;
  /** \brief Iterator over Vertices.*/
  //typedef unspecified    Space_vertex_iterator ;
  /** \brief Range of Vertices.*/
  //typedef unspecified    Space_vertex_range    ;
  /** \brief Returns a range over all Vertices.
  *
  * A MetricSpace must allow to traverse all Vertices attached
  * to the Points of the MetricSpace.
  */
  //Vertex_range space_vertex_range();

/** \brief Returns a specific Vertex which is different from all
  * the Vertices represented in the MetricSpace.
  *
  * The output of null_vertex() is comparable with < with the other
  * Vertices.
  * Computing a distance involving the output of null_vertex() will
  * result in an undefined behavior.*/
  Vertex null_vertex();

/** \brief Distance value type.
  *
  * Must be totally ordered with < and allow field operations
  * with +, -, * and /. */
  typedef unspecified FT;
  
/** \brief Distance function between two Vertices.*/
  FT distance( Vertex u, Vertex v );

};
