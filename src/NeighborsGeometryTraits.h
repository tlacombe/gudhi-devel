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
 * \extends MetricSpace
 */
struct NeighborsGeometryTraits
{
	/// \name Vertex:
	/// @{
	/**
	 * \brief Vertex type.
	 *
	 * A Vertex is the unique combinatorial label assigned to
	 * a metric Point from MetricSpace. The set of Vertices 
	 * must be totally ordered with <.
	 */
	typedef unspecified							Vertex;
	
	/** \brief Iterator over Vertices.*/
	typedef unspecified							Vertex_iterator;
	
	/** \brief Range of Vertices.*/
	typedef unspecified							Vertex_range;
	
	/** \brief Returns a range over all Vertices.*/
	Vertex_range vertex_range();
	/// @}
	
	/// \name Neighbor relationship:
	/// @{
	/**
	 * \brief An iterator of Vertices to traverse the geometric neighbors
	 * of a Point.
	 */
	typedef unspecified Neighbor_vertex_iterator;
	
	/** \brief Range of Vertices for a set of neighbors.*/
	typedef unspecified Neighbor_vertex_range;
	
	/** \brief Returns the range of neighbors of a Vertex.
	 *
	 * Arbitrary order.
	 * This defines the set of edges in a flag complex and, by extension,
	 * defines the flag complex.
	 * There is no specified order in the traversal.
	 */
//	typedef neighbor_range(Vertex v);
	/// @}
	
	/// \name Convertion methods Metric Space <-> Combinatorics:
	/// @{
	/**
	 * \brief Returns the point corresponding to a Vertex.
	 */
	Point vertex_to_point(const Vertex v);
	
	/**
	 * \brief Returns the Vertex assigned to a Point.
	 *
	 * \todo Necessary?
	 */
	Vertex point_to_vertex(const Point p);
	/// @}
	
};
