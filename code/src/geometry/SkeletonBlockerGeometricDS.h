/*
 * SkeletonBlockerGeometricDS.h
 *
 *  Created on: Feb 20, 2014 
 *      Author: David Salinas
 *  Copyright 2013 INRIA. All rights reserved
 */

#ifndef GUDHI_SKELETONBLOCKERGEOMETRICDS_H_
#define GUDHI_SKELETONBLOCKERGEOMETRICDS_H_

/** \brief Concept that must be passed to
 * the template class Skeleton_blocker_geometric_complex
 *
 */
template<typename GT>
struct SkeletonBlockerGeometricDS : public SkeletonBlockerDS
{
	typedef GT GeometryTrait;
	typedef typename GeometryTrait::Point Point;

	class Graph_vertex : public SkeletonBlockerDS::Graph_vertex{
	private:
		Point point_;
		Point& point(){	return point_; }
		const Point& point() const {	return point_; }
	};
};



#endif /* GUDHI_SKELETONBLOCKERGEOMETRICDS_H_ */
