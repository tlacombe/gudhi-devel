/*
 * Skeleton_blocker_geometric_complex.h
 *
 *  Created on: Feb 11, 2014
 *      Author: dsalinas
 */

#ifndef SKELETON_BLOCKER_GEOMETRIC_COMPLEX_H_
#define SKELETON_BLOCKER_GEOMETRIC_COMPLEX_H_


#include "utils/Utils.h"
#include "combinatorics/Skeleton_blocker/Simplifiable_skeleton_blocker.h"

template<typename ComplexGeometricTraits>
class Skeleton_blocker_geometric_complex : public Simplifiable_Skeleton_blocker<ComplexGeometricTraits>
{
public:

	typedef typename ComplexGeometricTraits::GeometryTrait GeometryTrait;

	typedef Simplifiable_Skeleton_blocker<ComplexGeometricTraits> SimplifiableSkeletonblocker;

	typedef typename SimplifiableSkeletonblocker::Vertex_handle Vertex_handle;
	typedef typename SimplifiableSkeletonblocker::Root_vertex_handle Root_vertex_handle;

	typedef typename SimplifiableSkeletonblocker::Vertex Vertex;

	typedef typename ComplexGeometricTraits::Point Point;


	void add_vertex(const Point& point){
		Vertex_handle ad = SimplifiableSkeletonblocker::add_vertex();
		(*this)[ad].point() = point;
	}

	std::string vertices_to_string() {
		std::ostringstream stream;
		for(auto vertex : this->vertex_range()){
			stream << (*this)[vertex].get_id().vertex<<" -- ";

//			stream<<"(";
//			for (auto x : (*this)[vertex].point()){
//				stream<<x<<",";
//			}
//			stream<<")"<<std::endl;

			//<<(*this)[vertex].point()<<")\n";
		}
		stream<< std::endl;
		return stream.str();
	}

};



#endif /* SKELETON_BLOCKER_GEOMETRIC_COMPLEX_H_ */
