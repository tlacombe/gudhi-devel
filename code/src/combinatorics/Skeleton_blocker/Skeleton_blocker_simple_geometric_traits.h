/*
 * Skeleton_blocker_simple_geometric_traits.h
 *
 *  Created on: Feb 11, 2014
 *      Author: dsalinas
 */

#ifndef GUDHI_SKELETON_BLOCKERS_SIMPLE_GEOMETRIC_TRAITS_H_
#define GUDHI_SKELETON_BLOCKERS_SIMPLE_GEOMETRIC_TRAITS_H_

#include <string>
#include <sstream>
#include "Skeleton_blocker_simple_traits.h"

template<typename GT>
struct Skeleton_blocker_simple_geometric_traits : public Skeleton_blocker_simple_traits {
public:

	typedef GT GeometryTrait;
	typedef typename GeometryTrait::Point Point;
	typedef typename Skeleton_blocker_simple_traits::Root_vertex_handle Root_vertex_handle;
	typedef typename Skeleton_blocker_simple_traits::Graph_vertex Simple_vertex;

	class Simple_geometric_vertex : public Simple_vertex{
		template<class ComplexGeometricTraits> friend class Skeleton_blocker_geometric_complex;
	private:
		Point point_;
		Point& point(){	return point_; }
		const Point& point() const {	return point_; }
	};

	typedef Simple_geometric_vertex Graph_vertex;
	typedef Skeleton_blocker_simple_traits::Graph_edge Graph_edge;
};


#endif /* GUDHI_SKELETON_BLOCKERS_SIMPLE_GEOMETRIC_TRAITS_H_ */
