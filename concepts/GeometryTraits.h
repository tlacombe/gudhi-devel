/*
 *  GeometryTraits.h
 *  GUDHI
 *
 *  Created by Cl??ment Maria on 12/10/13.
 *  Copyright 2013 __MyCompanyName__. All rights reserved.
 *
 */

/**
 Defines a metric space
 */
struct GeometryTraits
{
	/**
	 Element of a metric space, where pairwise 
	 distances can be computed
	 */
	typedef unspecified Point;
	
	/**
	 Comparable distance type, field operations
	 */
	typedef unspecified FT;
	
	/**
	 Distance function between points
	 */
	FT distance(Point p1, Point p2);
	
};
