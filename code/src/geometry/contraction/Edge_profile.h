/*
 * Edge_profile.h
 *
 *  Created on: Feb 13, 2014
 *      Author: dsalinas
 */

#ifndef GUDHI_EDGE_PROFILE_H_
#define GUDHI_EDGE_PROFILE_H_
//#include "combinatorics/Skeleton_blocker/Simplex.h"



namespace contraction {
template<typename GeometricSimplifiableComplex> class Edge_profile{

public:
	typedef GeometricSimplifiableComplex Complex;
	typedef typename Complex::GT GT;

	typedef typename GeometricSimplifiableComplex::Vertex_handle Vertex_handle;
	typedef typename GeometricSimplifiableComplex::Root_vertex_handle Root_vertex_handle;


	typedef typename GeometricSimplifiableComplex::Edge_handle edge_descriptor;
	typedef typename GeometricSimplifiableComplex::Graph_vertex Graph_vertex;
	typedef typename GeometricSimplifiableComplex::Graph_edge Graph_edge;
	typedef typename GeometricSimplifiableComplex::Point Point;




	Edge_profile( GeometricSimplifiableComplex& complex,edge_descriptor edge):complex_(complex),edge_handle_(edge)
{}

	GeometricSimplifiableComplex& complex() const {
		return complex_;
	}

	edge_descriptor edge_handle() const{
		return edge_handle_;
	}

	Graph_edge& edge() const{
		return complex_[edge_handle_];
	}


	Graph_vertex& v0() const{return complex_[v0_handle()];}
	Graph_vertex& v1() const{return complex_[v1_handle()];}


	Vertex_handle v0_handle() const{
		Root_vertex_handle root = complex_[edge_handle_].first();
		return *complex_.get_address(root);
	}

	Vertex_handle v1_handle() const{
		Root_vertex_handle root = complex_[edge_handle_].second();
		return *complex_.get_address(root);
	}

	const Point& p0() const {return complex_.point(v0_handle());}

	const Point& p1() const {return complex_.point(v1_handle());}

	friend std::ostream& operator << (std::ostream& o, const Edge_profile & v){
		o << "v0:"<<v.v0_handle() << " v1:"<<v.v1_handle();
		return o;
	}
private:

	GeometricSimplifiableComplex& complex_;

	edge_descriptor edge_handle_;

};


}  // namespace contraction

#endif /* GUDHI_EDGE_PROFILE_H_ */
