/*
 * Simple_skeleton_blockers_traits.h
 *
 *  Created on: Feb 11, 2014
 *      Author: dsalinas
 */

#ifndef GUDHI_SIMPLE_SKELETON_BLOCKERS_TRAITS_H_
#define GUDHI_SIMPLE_SKELETON_BLOCKERS_TRAITS_H_

#include <string>
#include <sstream>
#include "Simplex.h"





struct Simple_complex_DS_traits{
	/**
	 * @brief global and local handle similar to <a href="http://www.boost.org/doc/libs/1_38_0/libs/graph/doc/subgraph.html">boost subgraphs</a>.
	 * In gross, vertices are stored in a vector.
	 * For the root simplicial complex, the local and global descriptors are the same.
	 * For a subcomplex and one of its vertices 'v', its local descriptor is its position is the vector of these subcomplex
	 * whereas its global address is the position of 'v' in the root simplicial complex.
	 */
	struct Root_vertex_handle{
		typedef int boost_vertex_handle;
		explicit Root_vertex_handle(boost_vertex_handle val = -1 ):vertex(val){}
		boost_vertex_handle vertex;

		bool operator==( const Root_vertex_handle& other) const{
			return this->vertex == other.vertex;
		}

		bool operator<( const Root_vertex_handle& other) const{
			return this->vertex < other.vertex;
		}

		friend ostream& operator << (ostream& o, const Root_vertex_handle & v)
		{
			o << v.vertex;
			return o;
		}
	};

	struct Vertex_handle{
		typedef int boost_vertex_handle;
		Vertex_handle(boost_vertex_handle val=-1):vertex(val){}

		boost_vertex_handle vertex;

		bool operator==( const Vertex_handle& other) const{
			return this->vertex == other.vertex;
		}

		bool operator<(const Vertex_handle& other) const{
			return this->vertex < other.vertex;
		}

		friend ostream& operator << (ostream& o, const Vertex_handle & v)
		{
			o << v.vertex;
			return o;
		}
	};

	//typedef Simple_vertex<Root_vertex_handle> Vertex;


	class Simple_vertex {
		bool is_active_;
		Root_vertex_handle id_;
	public:
		virtual ~Simple_vertex(){}

		void activate(){is_active_=true;}
		void deactivate(){is_active_=false;}
		bool is_active() const{return is_active_;}
		void set_id(Root_vertex_handle i){id_=i;}
		Root_vertex_handle get_id() const{return id_;}

		virtual string to_string() const {
			ostringstream res;
			res<< id_;
			return res.str();
		}

		friend ostream& operator << (ostream& o, const Simple_vertex & v){
			o << v.to_string();
			return o;
		}
	};
	typedef Simple_vertex Vertex;


	// todo pkoi on stocke les id plutot que des addresses??
	// urgent a voir si on aurait pas mieux fait de stocker des addresses
	class Simple_edge {
		Root_vertex_handle a_;
		Root_vertex_handle b_;
		int id_;
	public:

		Simple_edge():a_(-1),b_(-1),id_(-1){}

		int& id(){return id_;}
		int id() const {return id_;}

		void setId(Root_vertex_handle a,Root_vertex_handle b){
			a_ = a;
			b_ = b;
		}

		Root_vertex_handle first() const {
			return a_;
		}

		Root_vertex_handle second() const {
			return b_;
		}

		friend ostream& operator << (ostream& o, const Simple_edge & v){
			o << "("<<v.a_<<","<<v.b_<<" - id = "<<v.id();
			return o;
		}
	};

	typedef Simple_edge Edge;


	typedef Simplex<Vertex_handle> Simplex_handle;
	typedef Simplex<Root_vertex_handle> Root_simplex_handle;

};



template<typename GT>
struct Simple_complex_geometry_traits : public Simple_complex_DS_traits {
public:

	typedef GT GeometryTrait;
	typedef typename GeometryTrait::Point Point;
	typedef typename Simple_complex_DS_traits::Root_vertex_handle Root_vertex_handle;
	typedef typename Simple_complex_DS_traits::Vertex Simple_vertex;

	class Simple_geometric_vertex : public Simple_vertex{
	private:
		Point point_;
	public:
		Point& point(){	return point_; }

		const Point& point() const {	return point_; }
	};

	typedef Simple_geometric_vertex Vertex;
	typedef Simple_complex_DS_traits::Edge Edge;
};



#endif /* GUDHI_SIMPLE_SKELETON_BLOCKERS_TRAITS_H_ */
