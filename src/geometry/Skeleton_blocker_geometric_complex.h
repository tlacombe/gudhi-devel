/*
 * Skeleton_blocker_geometric_complex.h
 *
 *  Created on: Feb 11, 2014
 *      Author: dsalinas
 */

#ifndef SKELETON_BLOCKER_GEOMETRIC_COMPLEX_H_
#define SKELETON_BLOCKER_GEOMETRIC_COMPLEX_H_


#include "utils/Utils.h"
#include "combinatorics/Skeleton_blocker/Skeleton_blocker_simplifiable_complex.h"

/**
 * @brief Class that represents a geometric complex that can be simplified.
 * The class allows access to points of vertices.
 *
 */
template<typename SkeletonBlockerGeometricDS>
class Skeleton_blocker_geometric_complex : public Skeleton_blocker_simplifiable_complex<SkeletonBlockerGeometricDS>
{
public:

	typedef typename SkeletonBlockerGeometricDS::GT GT;

	typedef Skeleton_blocker_simplifiable_complex<SkeletonBlockerGeometricDS> SimplifiableSkeletonblocker;

	typedef typename SimplifiableSkeletonblocker::Vertex_handle Vertex_handle;
	typedef typename SimplifiableSkeletonblocker::Root_vertex_handle Root_vertex_handle;

	typedef typename SimplifiableSkeletonblocker::Graph_vertex Graph_vertex;

	typedef typename SkeletonBlockerGeometricDS::Point Point;


	/**
	 * @brief Add a vertex to the complex with its associated point.
	 */
	void add_vertex(const Point& point){
		Vertex_handle ad = SimplifiableSkeletonblocker::add_vertex();
		(*this)[ad].point() = point;
	}


	/**
	 * @brief Returns the Point associated to the vertex v.
	 */
	const Point& point(Vertex_handle v) const{
		return (*this)[v].point();
	}

	/**
	 * @brief Returns the Point associated to the vertex v.
	 */
	Point& point(Vertex_handle v) {
		return (*this)[v].point();
	}

	const Point& point(Root_vertex_handle global_v) const{
		Vertex_handle local_v ( (*this)[global_v]) ;
		return (*this)[local_v].point();
	}

	Point& point(Root_vertex_handle global_v) {
		Vertex_handle local_v ( (*this)[global_v]) ;
		return (*this)[local_v].point();
	}


	//	std::string vertices_to_string() {
	//		std::ostringstream stream;
	//		for(auto vertex : this->vertex_range()){
	//			stream << (*this)[vertex].get_id().vertex<<" -- ";
	//
	////			stream<<"(";
	////			for (auto x : (*this)[vertex].point()){
	////				stream<<x<<",";
	////			}
	////			stream<<")"<<std::endl;
	//
	//			//<<(*this)[vertex].point()<<")\n";
	//		}
	//		stream<< std::endl;
	//		return stream.str();
	//	}

};



#endif /* SKELETON_BLOCKER_GEOMETRIC_COMPLEX_H_ */
