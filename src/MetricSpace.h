/*
 *  MetricSpace.h
 *  GUDHI
 *
 *  Created by Clément Maria on 12/10/13.
 *  Copyright 2013 INRIA. All rights reserved.
 *
 */

/**
 * \brief Defines a metric space. A metric space contains
 * elements (type Point) and pairwise distances can be computed.
 *
 * Remark: the distance may be more general, like a similarity m
 * measure, and simply needs to be symmetric.
 * 
 * \todo Defines the exact properties of the distance function.
 */
struct MetricSpace
{
	/// \name Point:
	/// @{
	/** \brief Element of a metric space.*/
	typedef unspecified Point;
	
	/** \brief Iterator over Points.*/
	typedef unspecified Point_iterator;
	
	/** \brief Range of Points.*/
	typedef unspecified Point_range;
	
	/** \brief Returns a range over all points in the MetricSpace.*/
	Point_range point_range();
	/// @}
	
	/// \name Distance:
	/// @{
	/**
	 * \brief Distance value type.
	 *
	 * Totally ordered with <. Field operations with +, -, * and /.
	 */
	typedef unspecified FT;
	
	/** \brief Distance function between points.*/
	FT distance(Point p1, Point p2);
	/// @}
	
};
