/*
 * Skeleton_blocker_contractor.h
 *
 *  Created on: Feb 11, 2014
 *      Author: dsalinas
 */

#ifndef GUDHI_SKELETON_BLOCKER_CONTRACTOR_H_
#define GUDHI_SKELETON_BLOCKER_CONTRACTOR_H_

#include <cassert>

// todo remove the queue to be independent from cgal
#include <CGAL/Modifiable_priority_queue.h>
#include <boost/scoped_array.hpp>
#include <boost/scoped_ptr.hpp>

#include "Edge_profile.h"
#include "policies/Cost_policy.h"
#include "policies/Edge_length_cost.h" //xxx remove
#include "policies/Placement_policy.h"
#include "policies/Middle_placement.h" //xxx remove

#include "policies/Valid_contraction_policy.h"
#include "policies/Dummy_valid_contraction.h" //xxx remove
#include "policies/Link_condition_valid_contraction.h" //xxx remove

#include "policies/Contraction_visitor.h"

#include "combinatorics/Skeleton_blocker/Skeleton_blocker_complex_visitor.h"

#include "utils/Utils.h"


namespace contraction {

/**
 *@class Skeleton_blocker_contractor
 *@brief Class that allows to contract iteratively edges of a simplicial complex.
 *
 * @details Basically, the simplification algorithm consists in iteratively picking the
 * edge with lowest cost and performing an edge contraction if the contraction is valid.
 * This class is policy based (and much inspired from the edge collapse package of CGAL).
 *
 * Policies that can be changed are :
 *  - the cost policy : how much cost an edge contraction
 *  - the placement policy : where will be placed the contraction point
 *  - the valid contraction policy : is the contraction valid. For instance, it can be
 *  a topological condition (link condition) or a geometrical condition (normals messed up).
 *
 * TODO expliquer la pile
 */
//@TODO constants iterator
template<class GeometricSimplifiableComplex
//   ,class ShouldStop
//        ,class Visitor
>
class Skeleton_blocker_contractor : public Dummy_complex_visitor<typename GeometricSimplifiableComplex::Vertex_handle>{

	GeometricSimplifiableComplex& complex_;

public:
	typedef typename GeometricSimplifiableComplex::Vertex Vertex;
	typedef typename GeometricSimplifiableComplex::Vertex_handle Vertex_handle;
	typedef typename GeometricSimplifiableComplex::Simplex_handle Simplex_handle;
	typedef typename GeometricSimplifiableComplex::Simplex_handle_iterator Simplex_handle_iterator;



	typedef typename GeometricSimplifiableComplex::Root_vertex_handle Root_vertex_handle;

	typedef typename GeometricSimplifiableComplex::Edge Edge;
	typedef typename GeometricSimplifiableComplex::Edge_handle edge_descriptor;
	typedef typename GeometricSimplifiableComplex::Complex_edge_iterator EdgeIterator;
	typedef typename GeometricSimplifiableComplex::Point Point;

	typedef Edge_profile<GeometricSimplifiableComplex> Profile;


	typedef Cost_policy<Profile> Cost_policy_;
	typedef Placement_policy<Profile> Placement_policy_;
	typedef Valid_contraction_policy<Profile> Valid_contraction_policy_;




	typedef Contraction_visitor<GeometricSimplifiableComplex> Contraction_visitor_;

	typedef boost::optional<double> Cost_type;
	typedef boost::optional<Point> Placement_type ;

	typedef size_t size_type;

	typedef Skeleton_blocker_contractor Self ;



private:

	struct Compare_id
	{
		Compare_id() : algorithm_(0) {}

		Compare_id( Self const* aAlgorithm ) : algorithm_(aAlgorithm) {}

		bool operator() ( edge_descriptor const& a, edge_descriptor const& b ) const
		{
			return algorithm_->get_undirected_edge_id(a) < algorithm_->get_undirected_edge_id(b);
		}

		Self const* algorithm_ ;
	} ;

	struct Compare_cost
	{
		Compare_cost() : algorithm_(0) {}

		Compare_cost( Self const* aAlgorithm ) : algorithm_(aAlgorithm) {}

		bool operator() ( edge_descriptor const& a, edge_descriptor const& b ) const
		{
			// NOTE: A cost is an optional<> value.
			// Absent optionals are ordered first; that is, "none < T" and "T > none" for any defined T != none.
			// In consequence, edges with undefined costs will be promoted to the top of the priority queue and poped out first.
			return algorithm_->get_data(a).cost() < algorithm_->get_data(b).cost();
		}

		Self const* algorithm_ ;
	} ;

	struct Undirected_edge_id : boost::put_get_helper<size_type, Undirected_edge_id>
	{
		typedef boost::readable_property_map_tag category;
		typedef size_type                        value_type;
		typedef size_type                        reference;
		typedef edge_descriptor                  key_type;

		Undirected_edge_id() : algorithm_(0) {}

		Undirected_edge_id( Self const* aAlgorithm ) : algorithm_(aAlgorithm) {}

		size_type operator[] ( edge_descriptor const& e ) const { return algorithm_->get_undirected_edge_id(e); }

		Self const* algorithm_ ;
	} ;

	typedef CGAL::Modifiable_priority_queue<edge_descriptor,Compare_cost,Undirected_edge_id> PQ ;
	typedef typename PQ::handle pq_handle ;


	// An Edge_data is associated with EVERY edge in the complex (collapsable or not).
	// It relates the edge with the PQ-handle needed to update the priority queue
	// It also relates the edge with a policy-based cache
	class Edge_data
	{
	public :

		Edge_data() : PQHandle_(),cost_() {}

		Cost_type const& cost() const { return cost_ ; }
		Cost_type      & cost()       { return cost_ ; }

		pq_handle PQ_handle() const { return PQHandle_ ;}

		bool is_in_PQ() const { return PQHandle_ != PQ::null_handle() ; }

		void set_PQ_handle( pq_handle h ) { PQHandle_ = h ; }

		void reset_PQ_handle() { PQHandle_ = PQ::null_handle() ; }

	private:
		pq_handle PQHandle_ ;
		Cost_type cost_ ;

	} ;
	typedef Edge_data* Edge_data_ptr ;
	typedef boost::scoped_array<Edge_data> Edge_data_array ;


	int get_undirected_edge_id ( edge_descriptor edge ) const {
		return complex_[edge].index() ;
	}

	Edge_data& get_data ( edge_descriptor const& edge ) const
	{
		return edge_data_array_[get_undirected_edge_id(edge)];
	}

	Cost_type get_cost(const Profile & profile){
		return (*cost_policy_)(profile,get_placement(profile));
	}

	Profile create_profile(edge_descriptor edge){
		return Profile(complex_,edge);
	}


	void insert_in_PQ( edge_descriptor const& edge, Edge_data& data )
	{
		data.set_PQ_handle(heap_PQ_->push(edge));
		++current_num_edges_heap_;
	}

	void update_in_PQ( edge_descriptor const& edge, Edge_data& data )
	{
		data.set_PQ_handle(heap_PQ_->update(edge,data.PQ_handle())) ;
	}

	void remove_from_PQ( edge_descriptor const& edge, Edge_data& data )
	{
		data.set_PQ_handle(heap_PQ_->erase(edge,data.PQ_handle()));
		--current_num_edges_heap_;
	}

	boost::optional<edge_descriptor> pop_from_PQ()	{
		boost::optional<edge_descriptor> rEdge = heap_PQ_->extract_top();
		if ( rEdge )
		{
			get_data(*rEdge).reset_PQ_handle();
		}
		return rEdge ;
	}





private:

	/**
	 * @brief Collect edges.
	 *
	 * Iterates over all edges of the simplicial complex and
	 * 1) inserts them in the priority queue sorted according to the Cost policy.
	 * 2) set the id() field of each edge
	 */
	void collect_edges(){
		//
		// Loop over all the edges in the complex in the heap
		//
		size_type lSize = complex_.num_edges();
		std::cerr  << "Collecting edges ..."<<std::endl;
		std::cerr  << lSize<<" edges "<<std::endl;

		edge_data_array_.reset( new Edge_data[lSize] ) ;

		heap_PQ_.reset( new PQ (lSize, Compare_cost(this), Undirected_edge_id(this) ) ) ;


		std::size_t id = 0 ;

		EdgeIterator edge_it;
		for(edge_it = complex_.edge_range().begin();
				edge_it != complex_.edge_range().end();
				++edge_it
		){
			edge_descriptor edge = *edge_it;
			complex_[edge].index() = id++;
			Profile const& lProfile = create_profile(edge);
			Edge_data& data = get_data(edge);
			data.cost() = get_cost(lProfile) ;
			++initial_num_edges_heap_;
			insert_in_PQ(edge,data);
			contraction_visitor_->on_collected(lProfile,data.cost());

		}
	}

	bool should_stop(double lCost,const Profile &profile){
		return false;
	}

	boost::optional<Point> get_placement(const Profile& profile){
		return (*placement_policy_)(profile);
	}

	bool is_contraction_valid( Profile const& profile, Placement_type placement ){
		return (*valid_contraction_policy_)(profile,placement);
	}


public:
	/**
	 * \brief Contract edges.
	 *
	 * While the heap is not empty, it extracts the edge with the minimum
	 * cost in the heap then try to contract it.
	 * It stops when the Stop policy says so or when the number of contractions
	 * given by 'num_max_contractions' is reached (if this number is positive).
	 */
	void contract_edges(int num_max_contractions=-1){

		DBG("Contract edges");
		DBGVALUE(complex_.num_vertices());
		int num_contraction = 0 ;
		//
		// Pops and processes each edge from the PQ
		//
		boost::optional<edge_descriptor> edge ;
		while ( (edge = pop_from_PQ())&& ((num_contraction<num_max_contractions)||(num_max_contractions<0)))
		{
			Profile const& profile = create_profile(*edge);
			Cost_type cost = get_data(*edge).cost();
			contraction_visitor_->on_selected(profile,cost,0,0);

			DBGMSG("---- Pop edge - num vertices :",complex_.num_vertices());

			if (cost)
			{
				DBGMSG("lCost",*cost);
				if (should_stop(*cost,profile) )
				{
					contraction_visitor_->on_stop_condition_reached(profile);
					DBG("should_stop");
					break ;
				}
				Placement_type placement = get_placement(profile);
				if ( is_contraction_valid(profile,placement) && placement )
				{
					DBG("contraction_valid");
					// The external function Get_new_vertex_point() is allowed to return
					// an absent point if there is no way to place the vertex
					// satisfying its constrians.
					// In that case the remaining vertex is simply left unmoved.
					contract_edge(profile,placement);
					++ num_contraction;
				}
				else
				{
					DBG("contraction not valid");
					contraction_visitor_->on_non_valid(profile);
				}
			}
			else
			{
				DBG("uncomputable cost");
			}
		}
	}

	/**
	 * @brief Returns an edge_descriptor and a Placement_type. This pair consists in
	 * the edge with the lowest cost in the heap together with its placement.
	 * The returned value is initialized iff the heap is non-empty.
	 */
	boost::optional<std::pair<edge_descriptor,Placement_type > > top_edge(){
		boost::optional<std::pair<edge_descriptor,Placement_type > > res;

		if(!heap_PQ_->empty()) {
			auto edge = heap_PQ_->top();
			Profile const& profile = create_profile(edge);
			Placement_type placement = get_placement(profile);
			res = make_pair(edge,placement);
			DBGMSG("top edge:",complex_[edge]);

		}
		return res;
	}

	Skeleton_blocker_contractor(GeometricSimplifiableComplex& complex)
	:complex_(complex),
	 cost_policy_(new Edge_length_cost<Profile>),
	 placement_policy_(new Middle_placement<Profile>),
	 valid_contraction_policy_(new Link_condition_valid_contraction<Profile>),
	 contraction_visitor_(new Contraction_visitor_()),
	 initial_num_edges_heap_(0),
	 current_num_edges_heap_(0)
	{
		contraction_visitor_->on_started(complex);
		complex_.set_visitor(this);
		collect_edges();
	}

	Skeleton_blocker_contractor(GeometricSimplifiableComplex& complex,
			Cost_policy_ *cost_policy_,
			Placement_policy_ * placement_policy_,
			Valid_contraction_policy_ * valid_contraction_policy_,
			Contraction_visitor_* contraction_visitor_ = new Contraction_visitor_()
	):
		complex_(complex),
		cost_policy_(cost_policy_),
		placement_policy_(placement_policy_),
		valid_contraction_policy_(valid_contraction_policy_),
		contraction_visitor_(contraction_visitor_),
		initial_num_edges_heap_(0),
		current_num_edges_heap_(0)
	{
		complex_.set_visitor(this);
		collect_edges();
	}





private:
	void contract_edge(const Profile& profile, Placement_type placement ) {
		contraction_visitor_->on_contracting(profile,placement);

		profile.v0().point() = *placement;
		profile.v1().point() = *placement; // remark optional since v1 would deactivated

		complex_.contract_edge(profile.v0_handle(),profile.v1_handle());

		// the visitor could do something as complex_.remove_popable_blockers();
		contraction_visitor_->on_contracted(profile,placement);
	}

	/**
	 * @brief we update the cost that has changed and the position in the heap
	 */
	void on_changed_edge(Vertex_handle a,Vertex_handle b) override{
		//1-get the edge_descriptor corresponding to ab
		//2-change the data in mEdgeArray[ab.id()]
		//3-update the heap
		edge_descriptor edge = (complex_[std::make_pair(a,b)]).first;
		Edge_data& data = get_data(edge);
		Profile const& profile = create_profile(edge);
		data.cost() = get_cost(profile) ;
		if ( data.is_in_PQ()){
			update_in_PQ(edge,data);
		}
		else{
			insert_in_PQ(edge,data);
		}
	}


	void on_remove_edge(Vertex_handle a,Vertex_handle b) override{

		edge_descriptor lEdge = (complex_[std::make_pair(a,b)]).first;
		Edge_data& lData = get_data(lEdge) ;
		if ( lData.is_in_PQ() )
		{
			remove_from_PQ(lEdge,lData) ;
		}
	}

	/**
	 * @brief Called when the edge 'ax' has been added while the edge 'bx'
	 * is still there but will be removed on next instruction.
	 * We assign the index of 'bx' to the edge index of 'ax'
	 */
	void on_swaped_edge(Vertex_handle a,Vertex_handle b,Vertex_handle x) override{
		std::pair<edge_descriptor,bool> ax_pair = complex_[std::make_pair(a,x)];
		std::pair<edge_descriptor,bool> bx_pair = complex_[std::make_pair(b,x)];
		assert(ax_pair.second && bx_pair.second);
		complex_[ax_pair.first].index() =complex_[bx_pair.first].index();
	}

	/**
	 * @brief Called when a blocker is removed.
	 * All the edges that passes through the blocker may be edge-contractible
	 * again and are thus reinserted in the heap.
	 */
	void on_delete_blocker(const Simplex_handle * blocker) override{
		// we go for all pairs xy that belongs to the blocker
		// note that such pairs xy are necessarily edges of the complex
		// by definition of a blocker
//		DBGVALUE(*blocker);

//		// boucle 2
//		//      gros bug -> pourquoi ce code loop???
//		// blocker est constant et ne devrais pas etre modifié pourtant
//		for ( Simplex_handle_iterator x = blocker->begin(); x!= blocker->end(); ++x){
//			DBGMSG("\n\nloopx, bl:",*blocker);
//			DBGMSG("loopx, x:",*x);
//
//			for(Simplex_handle_iterator y = x ; ++y != blocker->end(); ){
//				auto edge_descr = complex_[std::make_pair(*x,*y)].first;
//				Edge_data& data = get_data(edge_descr);
//				Profile const& profile = create_profile(edge_descr);
//
//				// cette ligne fait looper
//				data.cost() = get_cost(profile) ;
//				if ( !data.is_in_PQ() ){
//					insert_in_PQ(edge_descr,data);
//				}
//
//				DBGMSG("  loopy, x:",*x);
//				DBGMSG("  loopy, y:",*y);
//			}
//
//			DBGMSG("loopx end, x:",*x);
//		}

////		//boucle 3
////		// todo bug ce code ne loop pas mais moins efficace
		Simplex_handle blocker_copy(*blocker);
		for (auto x = blocker_copy.begin(); x!= blocker_copy.end(); ++x){
			for(auto y=x ; ++y != blocker_copy.end(); ){
				auto edge_descr = complex_[std::make_pair(*x,*y)].first;
				Edge_data& data = get_data(edge_descr);
				Profile const& profile = create_profile(edge_descr);
				data.cost() = get_cost(profile) ;

				// If the edge is already in the heap
				// its priority has not changed.
				// If the edge is not present, we reinsert it
				// remark : we could also reinsert the edge
				// only if it is valid
				if ( !data.is_in_PQ() ){
					insert_in_PQ(edge_descr,data);
				}
			}
		}
	}


private:

	Cost_policy_ * cost_policy_;
	Placement_policy_ * placement_policy_;
	Valid_contraction_policy_* valid_contraction_policy_;

	Contraction_visitor_* contraction_visitor_;

	Edge_data_array edge_data_array_ ;

	boost::scoped_ptr<PQ> heap_PQ_ ;
	int initial_num_edges_heap_;
	int current_num_edges_heap_;

};

}  // namespace contraction


#endif /* GUDHI_SKELETON_BLOCKER_CONTRACTOR_H_ */
