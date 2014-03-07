/*
 * Simplifiable_skeleton_blockers.h
 *
 *  Created on: Feb 4, 2014
 *      Author: dsalinas
 */

#ifndef GUDHI_SIMPLIFIABLE_SKELETON_BLOCKERS_H_
#define GUDHI_SIMPLIFIABLE_SKELETON_BLOCKERS_H_

#include "Skeleton_blocker_sub_complex.h"


/**
 *  \brief Class that allows simplification operation on a simplicial complex represented
 *  by a skeleton/blockers pair.
 */
template<typename Traits>
class Simplifiable_Skeleton_blocker : public Skeleton_blocker_complex<Traits>
{
	template<class ComplexType> friend class Skeleton_blocker_sub_complex;

public:


	typedef Skeleton_blocker_complex<Traits> SkeletonBlockerComplex;

	typedef typename SkeletonBlockerComplex::Edge Edge;

	typedef typename SkeletonBlockerComplex::boost_adjacency_iterator boost_adjacency_iterator;
	typedef typename SkeletonBlockerComplex::Edge_handle Edge_handle;
	typedef typename SkeletonBlockerComplex::boost_vertex_handle boost_vertex_handle;
	typedef typename SkeletonBlockerComplex::Vertex_handle Vertex_handle;
	typedef typename SkeletonBlockerComplex::Root_vertex_handle Root_vertex_handle;
	typedef typename SkeletonBlockerComplex::Simplex_handle Simplex_handle;
	typedef typename SkeletonBlockerComplex::Root_simplex_handle Root_simplex_handle;
	typedef typename SkeletonBlockerComplex::Blocker_handle Blocker_handle;


	typedef typename SkeletonBlockerComplex::Root_simplex_iterator Root_simplex_iterator;
	typedef typename SkeletonBlockerComplex::Simplex_handle_iterator Simplex_handle_iterator;
	typedef typename SkeletonBlockerComplex::BlockerMap BlockerMap;
	typedef typename SkeletonBlockerComplex::BlockerPair BlockerPair;
	typedef typename SkeletonBlockerComplex::BlockerMapIterator BlockerMapIterator;
	typedef typename SkeletonBlockerComplex::BlockerMapConstIterator BlockerMapConstIterator;

	typedef typename SkeletonBlockerComplex::Visitor Visitor;


	/** @name Constructors / Destructors / Initialization
	 */
	//@{
	Simplifiable_Skeleton_blocker(int num_vertices_ = 0,Visitor* visitor_=NULL):
		Skeleton_blocker_complex<Traits>(num_vertices_,visitor_){	}
	//@}



	virtual ~Simplifiable_Skeleton_blocker(){
	}

	/**
	 * Returns true iff the blocker 'sigma' is popable.
	 * To define popable, let us call 'L' the complex that
	 * consists in the current complex without the blocker 'sigma'.
	 * A blocker 'sigma' is then "popable" if the link of 'sigma'
	 * in L is reducible.
	 */
	virtual bool is_popable_blocker(Simplex_handle *sigma){
		this->remove_blocker(sigma);
		Skeleton_blocker_link_complex<Simplifiable_Skeleton_blocker<Traits> > link_blocker(*this,*sigma);
		this->add_blocker(sigma);
		return link_blocker.is_reducible()==CONTRACTIBLE;
	}



	/**
	 * Removes all the popable blockers of the complex and delete them.
	 * Doing so, it push edge that may have become contractible because of the suppression
	 * of some blockers.
	 * It also update the maximum dimension of blockers since it passes through every blockers
	 * @returns the number of popable blockers deleted
	 */
	int remove_popable_blockers(){
		int nb_blockers_deleted=0;
		bool blocker_popable_found=true;
		while (blocker_popable_found){
			blocker_popable_found=false;
			std::list <Blocker_handle> blockers = this->get_blockers_list();
			for (auto block : blockers)
			{
				if (is_popable_blocker(block)) {
					++nb_blockers_deleted;
					this->delete_blocker(block);
					blocker_popable_found = true;
					// some blockers may have become popable
					// xxx do something more efficient with a queue
					// that receives the neigborhood of an edge and
					// then the neigborhood of blockers that are poped.
					break;
				}
			}
		}
		return nb_blockers_deleted;
	}


	/**
	 * Remove the star of the vertex 'v'
	 */
	void remove_star(Vertex_handle v){
		while (this->degree(v) > 0)
		{
			Vertex_handle w( * (adjacent_vertices(v.vertex, this->skeleton).first));
			this->remove_edge(v,w);
		}
		this->remove_vertex(v);
	}

	/**
	 * Remove the star of the edge connecting vertices a and b.
	 * @returns the number of blocker that have been removed
	 */
	int remove_star(Vertex_handle a, Vertex_handle b){
		//	cout << "CollapseEdge("<<a<<","<<b<<")\n";
		std::list <Blocker_handle> blockers_to_delete;
		Edge_handle edge_descriptor;
		Simplex_handle edge(a,b);

		// we remove the blockers that are not consistent anymore
		for (auto blocker : this->blocker_range(a))
		{
			if (blocker->contains(b)){
				blockers_to_delete.push_back(blocker);
			}
		}
		int nb_blockers_removed = blockers_to_delete.size();
		while (!blockers_to_delete.empty())
		{
			this->delete_blocker(blockers_to_delete.back());
			blockers_to_delete.pop_back();
		}
		// we remove the edge
		this->remove_edge(a,b);

		return nb_blockers_removed;
	}




	/**
	 * Remove the star of the edge 'e'.
	 * @return true if the set of blockers has changed during contraction.
	 */
	bool remove_star(Edge_handle & e);

	/**
	 * Remove the star of the simplex 'sigma' which needs to belong to the complex
	 */
	void remove_star(Simplex_handle& sigma){
		if (sigma.dimension()==0){
			remove_star(sigma.first_vertex());
		}
		else if (sigma.dimension()==1){
			remove_star(sigma.first_vertex(),sigma.last_vertex());
		}
		else {
			std::cerr << "not implemented yet"<<std::endl;
			assert(false);
		}
	}

	/**
	 * Test if the complex is reducible using a strategy defined in the class
	 * (by default it tests if the complex is a cone)
	 * Note that NO could be returned if some invariant ensures that the complex
	 * is not a point (for instance if the euler characteristic is different from 1).
	 * This function will surely have to return MAYBE in some case because the
	 * associated problem is undecidable but it in practice, it can often
	 * be solved with the help of geometry.
	 */
	enum contractible_status{ NOT_CONTRACTIBLE,MAYBE_CONTRACTIBLE,CONTRACTIBLE, };
	virtual contractible_status is_reducible() const{
		if (this->is_cone()) return CONTRACTIBLE;
		else return MAYBE_CONTRACTIBLE;
//		return this->is_cone();
	}


	/** @Edge contraction operations
	 */
	//@{



private:
	/**
	 * @return true iff the link condition at edge ab is satisfied
	 * or equivalently iff no blocker contains ab.
	 */
	bool link_condition(Vertex_handle a, Vertex_handle b) const{
		for (auto blocker : this->const_blocker_range(a))
			if ( blocker->contains(b) ){
				return false;
			}
		return true;
	}

public:
	/**
	 * @return true iff the link condition at edge e is satisfied
	 * or equivalently iff no blocker contains e.
	 */
	bool link_condition(Edge_handle & e) const{
		const Edge& edge = (*this)[e];
		Vertex_handle a(*this->get_address(edge.first()));
		Vertex_handle b(*this->get_address(edge.second()));
		return link_condition(a,b);
	}

protected:
	/**
	 * Compute simplices beta such that a.beta is an order 0 blocker
	 * that may be used to construct a new blocker after contracting ab.
	 * Suppose that the link condition Link(ab) = Link(a) inter Link(b)
	 * is satisfied.
	 */
	void tip_blockers(Vertex_handle a, Vertex_handle b, std::vector<Simplex_handle> & buffer) const{
		for (auto const & blocker : this->const_blocker_range(a))
		{
			Simplex_handle beta = (*blocker);
			beta.remove_vertex(a);
			buffer.push_back(beta);
		}

		Simplex_handle n;
		this->add_neighbours(b,n);
		this->remove_neighbours(a,n);
		n.remove_vertex(a);


		for (Vertex_handle y : n)
		{
			Simplex_handle beta;
			beta.add_vertex( y );
			buffer.push_back(beta);
		}
	}


private:

	/**
	 * @brief "Replace" the edge 'bx' by the edge 'ax'.
	 * Assume that the edge 'bx' was present whereas 'ax' was not.
	 * Precisely, it does not replace edges, but remove 'bx' and then add 'ax'.
	 * The visitor 'on_swaped_edge' is called just after edge 'ax' had been added
	 * and just before edge 'bx' had been removed. That way, it can
	 * eventually access to information of 'ax'.
	 */
	void swap_edge(Vertex_handle a,Vertex_handle b,Vertex_handle x){
		this->add_edge(a,x);
		if (this->visitor) this->visitor->on_swaped_edge(a,b,x);
		this->remove_edge(b,x);
	}


private:

	/**
	 * @brief removes all blockers passing through the edge 'ab'
	 */
	void remove_blockers(Vertex_handle a, Vertex_handle b){
		std::vector<Blocker_handle> blocker_to_delete;
		for (auto blocker : this->blocker_range(a))
			if (blocker->contains(b)) blocker_to_delete.push_back(blocker);
		while (!blocker_to_delete.empty())
		{
			this->delete_blocker(blocker_to_delete.back());
			blocker_to_delete.pop_back();
		}
	}

public:
	/**
	 * Contracts the edge connecting vertices a and b.
	 * @return true if the set of blockers has changed during contraction.
	 * @remark If the link condition Link(ab) = Link(a) inter Link(b) is not satisfied,
	 * it removes first all blockers passing through 'ab'
	 */
	void contract_edge(Vertex_handle a, Vertex_handle b){
		if (a.vertex>b.vertex) std::swap(a,b);

		// if some blockers passes through 'ab', we remove them.
		if (!link_condition(a,b)){
			remove_blockers(a,b);
		}

		// Suppose ab is contracted to c
		// we first compute blockers of c of the type c.alpha.beta with:
		// 1) a.beta an order 0 blocker
		// 2) b.alpha an order 0 blocker
		// 3) every proper face of alpha.beta in Link(a) \cup Link(b)
		// 4) alpha.beta in K

		typedef Skeleton_blocker_link_complex<Skeleton_blocker_complex<Traits> > LinkComplexType;

		LinkComplexType link_a(*this,a);
		LinkComplexType link_b(*this,b);

		std::vector<Simplex_handle> vector_alpha, vector_beta;

		tip_blockers(a,b,vector_alpha);
		tip_blockers(b,a,vector_beta);

		std::vector<Simplex_handle *> blocker_to_add;
		for (auto alpha = vector_alpha.begin(); alpha != vector_alpha.end(); ++alpha){
			for (auto beta = vector_beta.begin(); beta != vector_beta.end(); ++beta)
			{
				Simplex_handle sigma = *alpha; sigma.union_vertices(*beta);
				Root_simplex_handle sigma_id = this->get_id(sigma);
				if ( this->contains(sigma) &&
						proper_faces_in_union<SkeletonBlockerComplex>(sigma_id,link_a,link_b))
				{
					Simplex_handle * blocker = new Simplex_handle(sigma);
					blocker->add_vertex(a);
					bool found=false;
					// we check that the blocker is not already there
					/**
					 * @todo TODO faire un set
					 */
					for(auto block : blocker_to_add)
					{
						if(*block==*blocker){
							found = true;
							break;
						}
					}
					if (!found)
						blocker_to_add.push_back(blocker);
				}
			}
		}

		// We delete all blockers a.tau and b.tau of K
		// after the edge contraction ab ----> a
		std::vector<Simplex_handle *> blocker_to_delete;
		for(auto bl : this->blocker_range(a))
			blocker_to_delete.push_back(bl);
		for(auto bl : this->blocker_range(b))
			blocker_to_delete.push_back(bl);
		while (!blocker_to_delete.empty())
		{
			this->delete_blocker(blocker_to_delete.back());
			blocker_to_delete.pop_back();
		}

		// We now proceed to the edge contraction ab ----> a,
		// modifying the simplicial complex

		// We update the set of edges
		this->remove_edge(a,b);


		// For all edges {b,x} incident to b,
		// we remove {b,x} and add {a,x} if not already there.
		while (this->degree(b)> 0)
		{
			Vertex_handle x(*(adjacent_vertices(b.vertex, this->skeleton).first));
			if(!this->contains_edge(a,x))
				// we 'replace' the edge 'bx' by the edge 'ax'
				this->swap_edge(a,b,x);
			else
				this->remove_edge(b,x);
		}

		// We remove b
		this->remove_vertex(b);
		// We notify the visitor that all edges incident to 'a' had changed
		boost_adjacency_iterator v, v_end;

		for (tie(v, v_end) = adjacent_vertices(a.vertex, this->skeleton); v != v_end; ++v)
			if (this->visitor) this->visitor->on_changed_edge(a,Vertex_handle(*v));

		// We add the new blockers through c
		for(auto block : blocker_to_add)
			this->add_blocker(block);
	}

	//@}


};



#endif /* GUDHI_SIMPLIFIABLE_SKELETON_BLOCKERS_H_ */
