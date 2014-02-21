/*
 * Edge_profile.h
 *
 *  Created on: Feb 13, 2014
 *      Author: dsalinas
 */

#ifndef EDGE_PROFILE_H_
#define EDGE_PROFILE_H_
#include "Simplex.h"



namespace contraction {
template<typename GeometricSimplifiableComplex> class Edge_profile{

public:
	typedef typename GeometricSimplifiableComplex::Vertex_handle Vertex_handle;
	typedef typename GeometricSimplifiableComplex::Root_vertex_handle Root_vertex_handle;


	typedef typename GeometricSimplifiableComplex::Edge_handle edge_descriptor;
	typedef typename GeometricSimplifiableComplex::Vertex Vertex;
	typedef typename GeometricSimplifiableComplex::Point Point;




	Edge_profile( GeometricSimplifiableComplex& complex,edge_descriptor edge):complex_(complex),edge_(edge)
{}

	GeometricSimplifiableComplex& complex() const {
		return complex_;
	}

	edge_descriptor edge() const{
		return edge_;
	}

	Vertex& v0() const{return complex_[v0_handle()];}

	Vertex& v1() const{return complex_[v1_handle()];}

	Vertex_handle v0_handle() const{
		Root_vertex_handle root = complex_[edge_].first();
		return *complex_.get_address(root);
	}

	Vertex_handle v1_handle() const{
		Root_vertex_handle root = complex_[edge_].second();
		return *complex_.get_address(root);
	}

	const Point& p0() const {return v0().point();}

	const Point& p1() const {return v1().point();}

	friend ostream& operator << (ostream& o, const Edge_profile & v){
		o << "v0:"<<v.v0_handle() << " v1:"<<v.v1_handle();
		return o;
	}
private:

	GeometricSimplifiableComplex& complex_;

	edge_descriptor edge_;

};


}  // namespace contraction

#endif /* EDGE_PROFILE_H_ */
