/*
 *  MetricSpace.h
 *  GUDHI
 *
 *  Created by Clément Maria on 12/10/13.
 *  Copyright 2013 INRIA. All rights reserved.
 *
 */

/** \brief Defines a discrete metric space. 
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
 * and \f$0\f$ otherwise.
 *
 * \todo Defines the exact properties of the distance function.*/
struct MetricSpace
{
  /// \name Vertex:
  /// @{
  /**
   * \brief Vertex type.
   *
   * A Vertex is the unique combinatorial label assigned to
   * a metric Point from MetricSpace. The set of Vertices 
   * must be totally ordered with <.
   *
   * A Vertex must be Default Constructible, Assignable, and Equality Comparable.
   */
  typedef unspecified    Vertex                ;
  /** \brief Iterator over Vertices.*/
  typedef unspecified    Space_vertex_iterator ;
  /** \brief Range of Vertices.*/
  typedef unspecified    Space_vertex_range    ;
  /** \brief Returns a range over all Vertices.
  *
  * A MetricSpace must allow to traverse all Vertices attached
  * to the Points of the MetricSpace.
  */
  Vertex_range space_vertex_range();
  /// @}

  /// \name Distance:
  /// @{
  /**
  * \brief Distance value type.
  *
  * Totally ordered with <. Field operations with +, -, * and /.
  */
  typedef unspecified              FT;
  /** 
  * \brief Distance function between the Points
  * corresponding to two Vertices.
  */
  FT distance( Vertex u, Vertex v );
  /// @}






/// \name Point:
/// @{
/** \brief Element of a metric space.
*
*
* \todo Necessary?
*/
//typedef unspecified Point;
/** \brief Iterator over Points.*/
//typedef unspecified Point_iterator;
/** \brief Range of Points.*/
//typedef unspecified Point_range;
/** \brief Returns a range over all Points in the MetricSpace.*/
//Point_range point_range();
 /** 
  * \brief Distance function between two Points.
  */
//  FT distance(Point p, Point q);
/// @}

};
