/*
 * Skeleton_blocker_complex.h
 *
 *  Created on: Jan, 2014
 *      Author: dsalinas
 */

#ifndef GUDHI_SKELETON_BLOCKER_COMPLEX_H
#define GUDHI_SKELETON_BLOCKER_COMPLEX_H

#include <map>
#include <vector>
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <memory>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>

#include "Skeleton_blocker_link_complex.h"
#include "Skeleton_blocker_link_superior.h"
#include "Skeleton_blocker_sub_complex.h"
#include "Simplex.h"

#include "Skeleton_blocker_complex_visitor.h"
#include "utils/Utils.h"

/**
 *@class Skeleton_blocker_complex
 *@brief Abstract Simplicial Complex represented with a skeleton/blockers pair.
 *
 * A simplicial complex is completely defined by :
 * - the graph of its 1-skeleton;
 * - its set of blockers.
 *
 * The graph is a boost graph templated with ComplexDS::Vertex and ComplexDS::Edge.
 *
 * One can access vertices through ComplexDS::Vertex_handle, edges through ComplexDS::Edge_handle and
 * simplices through Simplex_handle.
 * @todo TODO Simplex_handle are not classic handle
 *
 * The ComplexDS::Root_vertex_handle serves in the case of a subcomplex (see class Skeleton_blocker_sub_complex)
 * to access to the address of one vertex in the parent complex.
 *
 * @todo TODO constants iterator
 */
template<class ComplexDS>
class Skeleton_blocker_complex
{
	template<class ComplexType> friend class Skeleton_blocker_link_complex;
	template<class ComplexType> friend class Skeleton_blocker_link_superior;
	template<class ComplexType> friend class Skeleton_blocker_sub_complex;

public:


	typedef typename ComplexDS::Vertex Vertex;
	typedef typename ComplexDS::Edge Edge;

	typedef typename ComplexDS::Vertex_handle Vertex_handle;
	typedef typename ComplexDS::Root_vertex_handle Root_vertex_handle;
	typedef typename Root_vertex_handle::boost_vertex_handle boost_vertex_handle;
	typedef typename ComplexDS::Simplex_handle Simplex_handle;
	typedef typename ComplexDS::Root_simplex_handle Root_simplex_handle;
	typedef typename Root_simplex_handle::Simplex_vertex_const_iterator Root_simplex_iterator;
	typedef typename Simplex_handle::Simplex_vertex_const_iterator Simplex_handle_iterator;


protected:

	typedef typename boost::adjacency_list
			< boost::listS,
			boost::vecS,
			boost::undirectedS,
			Vertex,
			Edge
			> Graph;
	typedef typename boost::graph_traits<Graph>::vertex_iterator boost_vertex_iterator;
	typedef typename boost::graph_traits<Graph>::edge_iterator boost_edge_iterator;

protected:
	typedef typename boost::graph_traits<Graph>::adjacency_iterator boost_adjacency_iterator;

public:
	typedef typename boost::graph_traits<Graph>::edge_descriptor Edge_handle;



	typedef std::multimap<Vertex_handle,Simplex_handle *> BlockerMap;
	typedef typename std::multimap<Vertex_handle,Simplex_handle *>::value_type BlockerPair;
	typedef typename std::multimap<Vertex_handle,Simplex_handle *>::iterator BlockerMapIterator;
	typedef typename std::multimap<Vertex_handle,Simplex_handle *>::const_iterator BlockerMapConstIterator;

protected:
	int num_vertices_;
	int num_blockers_;

	typedef Skeleton_blocker_complex_visitor<Vertex_handle> Visitor;
	//	typedef Visitor* Visitor_ptr;
	Visitor* visitor;

	/**
	 * @brief If 'x' is a Vertex_handle of a vertex in the complex then degree[x] = d is its degree
	 * This quantity is updated when adding/removing edge (useful because the operation
	 * list.size() is done in linear time).
	 */
	std::vector<boost_vertex_handle> degree;
	Graph skeleton; /** 1-skeleton of the simplicial complex. */

	/** Each vertex can access to the blockers passing through it. */
	BlockerMap blocker_map;




public:

	/** @name Constructors / Destructors / Initialization
	 */
	//@{
	Skeleton_blocker_complex(int num_vertices_ = 0,Visitor* visitor_=NULL):visitor(visitor_){
		clear();
		for (int i=0; i<num_vertices_; ++i){
			add_vertex();
		}
	}

	Skeleton_blocker_complex(const Skeleton_blocker_complex& copy){
		clear();
		// xxx need an abstract factory or something to copy the visitor
		visitor = NULL;
		degree = copy.degree;
		skeleton = Graph(copy.skeleton);
		// we copy the blockers
		if (!copy.blocker_map.empty()){
			for (BlockerMapConstIterator dp=copy.blocker_map.begin(); dp!=copy.blocker_map.end(); ++dp){
				if ( (*dp).first == ( (*dp).second )->first_vertex() ){
					Simplex_handle* sigma = new Simplex_handle(*(*dp).second);
					add_blocker(sigma);
				}
			}
		}
	}

	Skeleton_blocker_complex& operator=(const Skeleton_blocker_complex& copy){
		clear();
		// xxx need an abstract factory or something to copy the visitor
		visitor = NULL;
		degree = copy.degree;
		skeleton = Graph(copy.skeleton);
		// we copy the blockers
		if (!copy.blocker_map.empty()){
			for (BlockerMapConstIterator dp=copy.blocker_map.begin(); dp!=copy.blocker_map.end(); ++dp){
				if ( (*dp).first == ( (*dp).second )->first_vertex() ){
					Simplex_handle* sigma = new Simplex_handle(*((*dp).second));
					add_blocker(sigma);
				}
			}
		}
		return *this;
	}


	/**
	 * The destructor delete all blockers allocated.
	 */
	virtual ~Skeleton_blocker_complex(){
		clear();
	}

	/**
	 * Clears the simplicial complex. After a call to this function,
	 * blockers are destroyed. The 1-skeleton and the set of blockers
	 * are both empty.
	 */
	virtual void clear(){
		// xxx for now the responsabilty of freeing the visitor is for
		// the user... not great do a shared_ptr instead
		visitor = NULL;

		degree.clear();
		num_vertices_ =0;

		// Desallocate the blockers

		while (!blocker_map.empty()){
			delete_blocker(blocker_map.begin()->second);
		}
		num_blockers_ = 0;


		blocker_map.clear();
		skeleton.clear();
	}

	void set_visitor(Visitor* other_visitor){
		visitor = other_visitor;
	}

	//@}


	/** @name Vertices operations
	 */
	//@{


public:

	/**
	 * @brief Return a local Vertex_handle of a vertex given a global one.
	 * @remark Assume that the vertex is present in the complex.
	 */
	Vertex_handle operator[](Root_vertex_handle global) const{
		auto local(get_address(global));
		assert(local);
		return *local;
	}

	Vertex& operator[](Vertex_handle address){
		return skeleton[address.vertex];
	}

	const Vertex& operator[](Vertex_handle address) const{
		return skeleton[address.vertex];
	}

	/**
	 * @brief Adds a vertex to the simplicial complex and returns its Vertex_handle.
	 */
	virtual Vertex_handle add_vertex(){
		Vertex_handle address(boost::add_vertex(skeleton));
		(*this)[address].activate();
		// safe since we now that we are in the root complex and the field 'address' and 'id'
		// are identical for every vertices
		(*this)[address].set_id(Root_vertex_handle(address.vertex));
		num_vertices_++;
		degree.push_back(0);
		if (visitor) visitor->on_add_vertex(address);
		return address;
	}


	/**
	 * @brief Remove a vertex from the simplicial complex
	 * @remark In fact, it just deactivates the vertex.
	 */
	void remove_vertex(Vertex_handle address){
		// We remove b
		boost::clear_vertex(address.vertex,skeleton);
		(*this)[address].deactivate();
		num_vertices_--;
		degree[address.vertex]=-1;
		if (visitor) visitor->on_remove_vertex(address);
	}

	/**
	 * @return true iff the simplicial complex contains the vertex u
	 */
	bool contains_vertex(Vertex_handle u) const{
		if (u.vertex<0 || u.vertex>=boost::num_vertices(skeleton)) return false;
		return (*this)[u].is_active();
	}

	/**
	 * @return true iff the simplicial complex contains the vertex u
	 */
	bool contains_vertex(Root_vertex_handle u) const{
		boost::optional<Vertex_handle> address = get_address(u);
		if (!address)
			return false;
		else
			return (*this)[address].is_active();
	}

	/**
	 * @return true iff the simplicial complex contains all vertices of simplex sigma
	 */
	bool contains_vertices(Simplex_handle* sigma) const{
		for (auto i = sigma->begin() ; i != sigma->end() ; ++i)
			if (!contains_vertex(*i))
				return false;
		return true;
	}


	/**
	 * Given an Id return the address of the vertex having this Id in the complex.
	 * For a simplicial complex, the address is the id but it may not be the case for a SubComplex.
	 *
	 */
	virtual boost::optional<Vertex_handle> get_address(Root_vertex_handle id) const{
		boost::optional<Vertex_handle> res;
		if ( id.vertex<= boost::num_vertices(skeleton) ) res = Vertex_handle(id.vertex);
		return res;
	}



	/**
	 * return the id of a vertex of adress local present in the graph
	 */
	Root_vertex_handle get_id(Vertex_handle local) const{
		return (*this)[local].get_id();
	}

	//@}

	/** @name Edges operations
	 */
	//@{

public:

	/**
	 * @brief returns a pair that consists in :
	 * - an edge handle that is valid and describes ab iff the edge is present.
	 * - a boolean which is true is the edge was founded.
	 */
	std::pair<Edge_handle,bool> operator[](const std::pair<Vertex_handle,Vertex_handle>& ab){
		return boost::edge(ab.first.vertex,ab.second.vertex,skeleton);
	}

	Edge& operator[](Edge_handle edge_descriptor){
		return skeleton[edge_descriptor];
	}

	const Edge& operator[](Edge_handle edge_descriptor) const{
		return skeleton[edge_descriptor];
	}

	/**
	 * @brief returns the simplex made with the two vertices of the edge
	 */
	Simplex_handle get_vertices(Edge_handle edge_descriptor) const{
		auto edge((*this)[edge_descriptor]);
		return Simplex_handle((*this)[edge.first()],(*this)[edge.second()]);
	}

	/**
	 * @brief Adds an edge between vertices a and b
	 */
	Edge_handle add_edge(Vertex_handle a, Vertex_handle b){
		assert(contains_vertex(a) && contains_vertex(b));
		std::pair<Edge_handle,bool> pair_descr_bool = (*this)[std::make_pair(a,b)];
		Edge_handle edge_descr;
		bool edge_present = pair_descr_bool.second;
		if (edge_present== false)
		{
			edge_descr = boost::add_edge(a.vertex,b.vertex,skeleton).first;
			(*this)[edge_descr].setId(get_id(a),get_id(b));
			degree[a.vertex]++;
			degree[b.vertex]++;
			if (visitor) visitor->on_add_edge(a,b);
		}
		return edge_descr;
	}

	/**
	 * @brief Adds all edges of simplex sigma to the simplicial complex.
	 */
	void add_edges(Simplex_handle & sigma){
		Simplex_handle_iterator i, j;
		for (i = sigma.begin() ; i != sigma.end() ; ++i)
			for (j = i, j++ ; j != sigma.end() ; ++j)
				add_edge(*i,*j);
	}

	/**
	 * @brief Removes edge ab from the simplicial complex.
	 */
	virtual Edge_handle remove_edge(Vertex_handle a, Vertex_handle b){
		bool found;
		Edge_handle edge;
		tie(edge,found) = boost::edge(a.vertex,b.vertex,skeleton);
		if (found)
		{
			if (visitor) visitor->on_remove_edge(a,b);
			//		if (heapCollapse.Contains(edge)) heapCollapse.Delete(edge);
			boost::remove_edge(a.vertex,b.vertex,skeleton);
			degree[a.vertex]--;
			degree[b.vertex]--;
		}
		return edge;
	}

	/**
	 * @return true iff the simplicial complex contains an edge between
	 * vertices a and b
	 */
	bool contains_edge(Vertex_handle a, Vertex_handle b) const{
		//if (a.vertex<0 || b.vertex <0) return false;
		return boost::edge(a.vertex,b.vertex,skeleton).second;
	}


	/**
	 * @return true iff the simplicial complex contains all vertices
	 * of simplex sigma
	 */
	bool contains_vertices(const Simplex_handle & sigma) const{
		for (auto vertex : sigma)
			if(!contains_vertex(vertex)) return false;
		return true;
	}

	/**
	 * @return true iff the simplicial complex contains all vertices
	 * and all edges of simplex sigma
	 */
	bool contains_edges(const Simplex_handle & sigma) const{
		for (auto i = sigma.begin() ; i != sigma.end() ; ++i){
			if(!contains_vertex(*i)) return false;
			auto j=i;
			for (++j ; j != sigma.end() ; ++j){
				if (!contains_edge(*i,*j))
					return false;
			}
		}
		return true;
	}
	//@}

	/** @name Blockers operations
	 */
	//@{
	/**
	 * Adds the 2-blocker abc
	 */
	void add_blocker(Vertex_handle a, Vertex_handle b, Vertex_handle c){
		Simplex_handle * sigma = new Simplex_handle(a,b,c);
		add_blocker(sigma);
		if (visitor) visitor->on_add_blocker(*sigma);
	}

	/**
	 * Adds the 3-blocker abcd
	 */
	void add_blocker(Vertex_handle a, Vertex_handle b, Vertex_handle c, Vertex_handle d){
		Simplex_handle * sigma = new Simplex_handle(a,b,c,d);
		add_blocker(sigma);
	}

	/**
	 * Adds the simplex s to the set of blockers
	 */
	void add_blocker(Simplex_handle * sigma){
		if (contains_blocker(*sigma))
		{
			//		cout << "ATTEMPT TO ADD A BLOCKER ALREADY THERE ---> BLOCKER IGNORED" << endl;
			return;
		}
		else{
			num_blockers_++;
			auto vertex = sigma->begin();
			while(vertex != sigma->end())
			{
				blocker_map.insert(BlockerPair(*vertex,sigma));
				++vertex;
			}
		}
	}

	/**
	 * Removes sigma from the blocker map of vertex v
	 */
	void remove_blocker(const Simplex_handle * sigma, Vertex_handle v){
		Complex_const_blocker_iterator blocker;
		for (
				blocker = const_blocker_range(v).begin();
				blocker != const_blocker_range(v).end();
				++blocker
		){
			if (*blocker == sigma) break;
		}
		if (*blocker != sigma){
			std::cerr << "bug ((*blocker).second == sigma) ie try to remove a blocker not present\n";
			assert(false);
		}
		else{
			blocker_map.erase(blocker.currentPosition);
		}
	}

	/**
	 * Removes the simplex s from the set of blockers.
	 * s has to belongs to the set of blockers
	 */
	void remove_blocker(const Simplex_handle * sigma){
		if (visitor) visitor->on_remove_blocker(*sigma);
		for (auto vertex : *sigma){
			remove_blocker(sigma,vertex);
		}
		num_blockers_--;
	}



	/**
	 * Removes the simplex s from the set of blockers
	 * and desallocate s.
	 */
	void delete_blocker(Simplex_handle * sigma){
		remove_blocker(sigma);
		delete sigma;
	}




	/**
	 * @return true iff s is a blocker of the simplicial complex
	 */
	bool contains_blocker(const Simplex_handle & s) const{
		if (s.dimension()<2)
			return false;

		Vertex_handle a = s.first_vertex();

		for ( Complex_const_blocker_iterator	blocker = const_blocker_range(a).begin();
				blocker != const_blocker_range(a).end();
				++blocker
		){
			if ( s == *((*blocker)) )
				return true;
		}


		return false;
	}


	/**
	 * @returns the list of blockers of the simplex
	 */
	std::list<Simplex_handle*> get_blockers_list(){
		std::list <Simplex_handle*> res;
		for (BlockerMapConstIterator dp=blocker_map.begin(); dp!=blocker_map.end(); ++dp){
			// check if it is the first time we encounter the blocker
			if ( (*dp).first == ( (*dp).second )->first_vertex() ){
				res.push_back((*dp).second);
			}
		}
		return res;
	}


private:
	/**
	 * @return true iff a blocker of the simplicial complex
	 * is a face of sigma.
	 */
	bool blocks(const Simplex_handle & sigma) const{
		BlockerMapConstIterator lb, ub, tau;
		for (auto vi : sigma)
		{
			// We start by listing all blockers that pass through
			// one vertex of sigma.
			lb=blocker_map.lower_bound(vi);
			ub=blocker_map.upper_bound(vi);
			for (tau = lb; tau != ub; ++tau)
				if ( (*tau).first == ( (*tau).second )->first_vertex() )
					if ( sigma.contains( *(*tau).second ) )
						return true;
		}

		return false;
	}

	//@}


	/** @name Neighbourhood access
	 */
	//@{


protected:
	/**
	 * @brief Adds to simplex n the neighbours of v:
	 * \f$ n \leftarrow n \cup N(v) \f$
	 * If 'keep_only_superior' is true then only vertices that are greater than v are added.
	 */
	virtual void add_neighbours(Vertex_handle v, Simplex_handle & n,bool keep_only_superior=false) const{
		boost_adjacency_iterator ai, ai_end;
		for (tie(ai, ai_end) = adjacent_vertices(v.vertex, skeleton); ai != ai_end; ++ai){
			if (keep_only_superior){
				if (*ai>v.vertex)
					n.add_vertex(Vertex_handle(*ai));
			}
			else
				n.add_vertex(Vertex_handle(*ai));
		}
	}

	/**
	 * @brief Add to simplex n all vertices which are
	 * neighbours of alpha: \f$ n \leftarrow n \cup N(alpha) \f$
	 * If 'keep_only_superior' is true then only vertices that are greater than alpha are added.
	 *
	 */
	virtual void add_neighbours(const Simplex_handle &alpha, Simplex_handle & n,bool keep_only_superior=false) const{
		n.clear();
		// ----------------------------
		// Compute vertices in the link
		// we compute the intersection of N(alpha_i) and store it in n
		// ----------------------------
		auto alpha_vertex = alpha.begin();
		add_neighbours(*alpha_vertex,n,keep_only_superior);
		for (alpha_vertex = (alpha.begin())++ ; alpha_vertex != alpha.end() ; ++alpha_vertex)
		{
			keep_neighbours(*alpha_vertex,n,keep_only_superior);
		}
	}

	/**
	 * @brief Eliminates from simplex n all vertices which are
	 * not neighbours of v: \f$ n \leftarrow n \cap N(v) \f$
	 */
	virtual void keep_neighbours(Vertex_handle v, Simplex_handle& n,bool keep_only_superior=false) const{
		Simplex_handle nv;
		add_neighbours(v,nv,keep_only_superior);
		n.intersection(nv);
	}

	/**
	 * @brief Eliminates from simplex n all vertices which are
	 * neighbours of v: \f$ n \leftarrow n \setminus N(v) \f$
	 */
	virtual void remove_neighbours(Vertex_handle v, Simplex_handle & n,bool keep_only_superior=false) const{
		Simplex_handle nv;
		add_neighbours(v,nv,keep_only_superior);
		n.difference(nv);
	}



	//@}


	/** @name Operations on the simplicial complex
	 */
	//@{
public:

	/**
	 * @brief Compute the local vertices of 's' in the current complex
	 * If one of them is not present in the complex then the return value is uninitialized.
	 */
	boost::optional<Simplex_handle>	get_simplex_address(const Root_simplex_handle& s) const
	{
		boost::optional<Simplex_handle> res;

		Simplex_handle s_address;
		//Root_simplex_const_iterator i;
		for (auto i = s.begin() ; i != s.end() ; ++i)
		{
			boost::optional<Vertex_handle> address = get_address(*i);
			if (!address)
				return res;
			else
				s_address.add_vertex(*address);
		}
		res = s_address;
		return res;
	}

	/**
	 * @return a simplex with vertices which are the id of vertices of the
	 * argument
	 */
	Root_simplex_handle get_id(const Simplex_handle& local_simplex) const{
		Root_simplex_handle global_simplex;
		for (auto x = local_simplex.begin(); x!= local_simplex.end();++x){
			global_simplex.add_vertex(get_id(*x));

		}
		return global_simplex;

	}


	/**
	 * @return true iff the simplex s belongs to the simplicial
	 * complex.
	 */
	virtual bool contains(const Simplex_handle & s) const{
		if (s.dimension() == -1 ) return false;
		else
			if (s.dimension() ==0 ){
				return contains_vertex(s.first_vertex());
			}
			else
				return ( contains_edges(s) && !blocks(s) );
	}

	/*
	 * @return true iff the complex is empty
	 */
	bool empty() const{
		return num_vertices()==0;
	}

	/*
	 * @return the number of vertices in the complex
	 */
	int num_vertices() const{
		return num_vertices_;
	}

	/*
	 * @return the number of edges in the complex
	 */
	int num_edges() const{
		return boost::num_edges(skeleton);
	}

	/*
	 * @return the number of blockers in the complex
	 */
	int num_blockers() const{
		return num_blockers_;
	}

	/*
	 * @return true iff the graph of the 1-skeleton of the complex is complete
	 */
	bool complete() const{
		return (num_vertices()*(num_vertices()-1))/2 == num_edges();
	}

	/**
	 * @return the number of connected components in the graph of the 1-skeleton
	 */
	int num_connected_components(){
		int num_vert_collapsed = skeleton.vertex_set().size() - num_vertices();
		std::vector<int> component(skeleton.vertex_set().size());
		return boost::connected_components(this->skeleton,&component[0]) - num_vert_collapsed;
	}

	/**
	 * Test if the complex is a cone
	 * runs in O(n) if n is the number of vertices
	 */
	bool is_cone() const{
		if (num_vertices()==0) return false;
		if (num_vertices()==1) return true;
		std::pair<boost_vertex_iterator, boost_vertex_iterator> vp;
		for (vp = vertices(skeleton); vp.first != vp.second; ++vp.first)
			if ((skeleton[*vp.first]).is_active()){
				Vertex_handle v(*vp.first);
				// we check if the current vertex belongs to a blocker
				if (blocker_map.find(v)==blocker_map.end()){
					//				if (out_degree(v, skeleton) == num_vertices() -1)
					if (degree[v.vertex] == num_vertices() -1)

						// we check if the current vertex is linked to all others vertices of the complex
						return true;
				}
			}
		return false;
	}

	/** @name Vertex, Edge, simplex and blockers iterators
	 */
	//@{

	//friend class Complex_vertex_iterator;
	/**
	 * \brief Range over the vertices of the simplicial complex.
	 * Methods .begin() and .end() return a Complex_vertex_iterator.
	 */
	class Complex_vertex_range;

	/**
	 * \brief Iterator over the simplices of a
	 * simplicial complex.
	 *
	 * 'value_type' must be Simplex_handle.
	 */
	class Complex_vertex_iterator;
	/**
	 * Returns a Complex_vertex_range over all vertices of the complex
	 */
	Complex_vertex_range vertex_range() const
	{return Complex_vertex_range(this);}

	/**
	 * \brief Range over the edges of the simplicial complex.
	 * Methods .begin() and .end() return a Complex_edge_iterator.
	 */
	class Complex_edge_range;
	/**
	 * \brief Iterator over the edges of the
	 * simplicial complex.
	 *
	 */
	class Complex_edge_iterator;

	/**
	 * Returns a Complex_edge_range over all edges of the
	 * simplicial complex
	 */
	Complex_edge_range edge_range() const
	{return Complex_edge_range(this);}



	/**
	 * \brief Range over triangle around a vertex of the simplicial complex
	 * Methods .begin() and .end() return a Triangle_around_vertex_iterator.
	 *
	 */
	template<typename LinkType>
	class Triangle_around_vertex_range;

	/**
	 * \brief Iterator over the triangle around a vertex 'v' of the simplicial complex.
	 * The template LinkType has to either Skeleton_blocker_link_complex
	 * or Skeleton_blocker_link_superior if one just wants the triangles whose
	 * vertices are greater than 'v'.
	 */
	template<typename LinkType>
	class Triangle_around_vertex_iterator;

	/**
	 * Returns a Triangle_around_vertex_range over all triangles around
	 * a vertex of the simplicial complex.
	 * The template LinkType has to either Skeleton_blocker_link_complex
	 * or Skeleton_blocker_link_superior if one just wants the triangles whose
	 * vertices are greater than 'v'.
	 */
	/*template<typename LinkType>
	Triangle_around_vertex_range<LinkType> triangle_range(Vertex_handle v) const
	{return Triangle_around_vertex_range<LinkType>(this,v);}
	 */

	typedef Skeleton_blocker_link_complex<Skeleton_blocker_complex<ComplexDS> > Link;
	typedef Skeleton_blocker_link_superior<Skeleton_blocker_complex<ComplexDS> > Superior_link;

	typedef Triangle_around_vertex_iterator<Superior_link> Superior_triangle_around_vertex_iterator;

	Triangle_around_vertex_range<Link> triangle_range(Vertex_handle v)
									{return Triangle_around_vertex_range<Link>(this,v);}

	Triangle_around_vertex_range<Superior_link> superior_triangle_range(Vertex_handle v)
									{return Triangle_around_vertex_range<Superior_link>(this,v);}




	//////////////
	/**
	 * \brief Range over triangles of the simplicial complex
	 * Methods .begin() and .end() return a Triangle_vertex_iterator.
	 *
	 */
	class Triangle_range;

	/**
	 * \brief Iterator over the triangles of the simplicial complex.
	 */
	class Triangle_iterator;

	/**
	 * Returns a Triangle_range over all triangles of the simplicial complex.
	 */
	Triangle_range triangle_range()
	{return Triangle_range(this);}

	/**
	 * \brief Range over the blockers of the simplicial complex adjacent to a vertex.
	 * Methods .cbegin() and .cend() return a Complex_const_blocker_iterator.
	 */
	friend class Complex_const_blocker_range;
	class Complex_const_blocker_range;

	/**
	 * \brief Iterator over the blockers adjacent to a vertex
	 *
	 */
	class Complex_const_blocker_iterator;
	/**
	 * Returns a Complex_simplex_range over all blockers of the complex adjacent to the vertex 'v'
	 */
	Complex_const_blocker_range const_blocker_range(Vertex_handle v) const
	{return Complex_const_blocker_range(this,v);}



public:
	//	/**
	//	 * \brief Range over the blockers of the simplicial complex adjacent to a vertex.
	//	 * Methods .begin() and .end() return a Complex_blocker_iterator.
	//	 */
	//template<typename ComplexType> friend class Complex_blocker_range;


	/**
	 * \brief Iterator over the blockers adjacent to a vertex
	 *
	 */
	template<typename MapType, typename ReturnType>
	class Complex_blocker_iterator;

	typedef Complex_blocker_iterator<BlockerMapIterator,Simplex_handle*> MyBlockerIterator;
	typedef Complex_blocker_iterator<BlockerMapConstIterator, Simplex_handle*> MyConstBlockerIterator;


	// ComplexType should be Skeleton_blocker_complex& or const Skeleton_blocker_complex&
	template <typename ComplexType>
	class Complex_blocker_range{
	private:
		Vertex_handle v_;
		ComplexType complex_;
	public:
		Complex_blocker_range(ComplexType complex,Vertex_handle v):v_(v),complex_(complex){}

		MyBlockerIterator begin(){
			return MyBlockerIterator(complex_.blocker_map.lower_bound(v_));

		}
		MyBlockerIterator end()	{
			return MyBlockerIterator(complex_.blocker_map.upper_bound(v_));
		}

		MyConstBlockerIterator begin() const{
			return MyConstBlockerIterator(complex_.blocker_map.lower_bound(v_));

		}
		MyConstBlockerIterator end() const{
			return MyConstBlockerIterator(complex_.blocker_map.upper_bound(v_));
		}
	};

	typedef Complex_blocker_range<Skeleton_blocker_complex&> ComplexBlockerRange;
	typedef Complex_blocker_range<const Skeleton_blocker_complex&> ConstComplexBlockerRange;

	/**
	 * Returns a Complex_simplex_range over all blockers of the complex adjacent to the vertex 'v'
	 */
	ConstComplexBlockerRange blocker_range(Vertex_handle v) const
	{return ConstComplexBlockerRange(*this,v);}

	ComplexBlockerRange blocker_range(Vertex_handle v)
	{return ComplexBlockerRange(*this,v);}

	//@}

	/** @name Print and IO methods
	 */
	//@{
public:
	std::string to_string(){
		std::ostringstream stream;
		stream<<num_vertices()<<" vertices:\n"<<vertices_to_string()<<std::endl;
		stream<<num_edges()<<" edges:\n"<<edges_to_string()<<std::endl;
		stream<<num_blockers()<<" blockers:\n"<<blockers_to_string()<<std::endl;
		return stream.str();
	}

	virtual std::string vertices_to_string() {
		std::ostringstream stream;
		for(auto vertex : vertex_range())
			stream << "("<<(*this)[vertex].get_id()<<"),";
		stream<< std::endl;
		return stream.str();
	}

	std::string edges_to_string() {
		std::ostringstream stream;
		for(auto edge : edge_range())
			stream <<edge<<" id = "<< (*this)[edge].id()<< std::endl;
		stream<< std::endl;
		return stream.str();
	}


	std::string blockers_to_string() const{
		std::ostringstream stream;
		for (auto bl:blocker_map){
			stream << bl.first << " => " << bl.second << ":"<<*bl.second <<"\n";
		}
		return stream.str();
	}

	//@}

};

//#include "src/AbstractComplex/SimplicialComplex.cxx"

#include "Skeleton_blocker_complex_iterators.h"



#endif

