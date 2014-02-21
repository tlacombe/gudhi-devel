/*
 * Skeleton_blocker_contractor.h
 *
 *  Created on: Feb 11, 2014
 *      Author: dsalinas
 */

#ifndef SKELETON_BLOCKER_CONTRACTOR_H_
#define SKELETON_BLOCKER_CONTRACTOR_H_

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
#include "Skeleton_blocker_complex_visitor.h"

#include "Utils.h"


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
	typedef typename GeometricSimplifiableComplex::Root_vertex_handle Root_vertex_handle;
	typedef typename GeometricSimplifiableComplex::Edge Edge;
	typedef typename GeometricSimplifiableComplex::Edge_handle edge_descriptor;
	typedef typename GeometricSimplifiableComplex::Complex_edge_iterator EdgeIterator;
	typedef typename GeometricSimplifiableComplex::Point Point;

	typedef Edge_profile<GeometricSimplifiableComplex> Profile;

	typedef boost::optional<double> Cost_type;
	typedef boost::optional<Point> Placement_type ;

	typedef size_t size_type;

	typedef Skeleton_blocker_contractor Self ;

private:

	struct Compare_id
	{
		Compare_id() : mAlgorithm(0) {}

		Compare_id( Self const* aAlgorithm ) : mAlgorithm(aAlgorithm) {}

		bool operator() ( edge_descriptor const& a, edge_descriptor const& b ) const
		{
			return mAlgorithm->get_undirected_edge_id(a) < mAlgorithm->get_undirected_edge_id(b);
		}

		Self const* mAlgorithm ;
	} ;

	struct Compare_cost
	{
		Compare_cost() : mAlgorithm(0) {}

		Compare_cost( Self const* aAlgorithm ) : mAlgorithm(aAlgorithm) {}

		bool operator() ( edge_descriptor const& a, edge_descriptor const& b ) const
		{
			// NOTE: A cost is an optional<> value.
			// Absent optionals are ordered first; that is, "none < T" and "T > none" for any defined T != none.
			// In consequence, edges with undefined costs will be promoted to the top of the priority queue and poped out first.
			return mAlgorithm->get_data(a).cost() < mAlgorithm->get_data(b).cost();
		}

		Self const* mAlgorithm ;
	} ;

	struct Undirected_edge_id : boost::put_get_helper<size_type, Undirected_edge_id>
	{
		typedef boost::readable_property_map_tag category;
		typedef size_type                        value_type;
		typedef size_type                        reference;
		typedef edge_descriptor                  key_type;

		Undirected_edge_id() : mAlgorithm(0) {}

		Undirected_edge_id( Self const* aAlgorithm ) : mAlgorithm(aAlgorithm) {}

		size_type operator[] ( edge_descriptor const& e ) const { return mAlgorithm->get_undirected_edge_id(e); }

		Self const* mAlgorithm ;
	} ;

	typedef CGAL::Modifiable_priority_queue<edge_descriptor,Compare_cost,Undirected_edge_id> PQ ;
	typedef typename PQ::handle pq_handle ;


	// An Edge_data is associated with EVERY edge in the complex (collapsable or not).
	// It relates the edge with the PQ-handle needed to update the priority queue
	// It also relates the edge with a policy-based cache
	class Edge_data
	{
	public :

		Edge_data() : mPQHandle(),mCost() {}

		Cost_type const& cost() const { return mCost ; }
		Cost_type      & cost()       { return mCost ; }

		pq_handle PQ_handle() const { return mPQHandle ;}

		bool is_in_PQ() const { return mPQHandle != PQ::null_handle() ; }

		void set_PQ_handle( pq_handle h ) { mPQHandle = h ; }

		void reset_PQ_handle() { mPQHandle = PQ::null_handle() ; }

	private:
		pq_handle mPQHandle ;
		Cost_type mCost ;

	} ;
	typedef Edge_data* Edge_data_ptr ;
	typedef boost::scoped_array<Edge_data> Edge_data_array ;


	int get_undirected_edge_id ( edge_descriptor aEdge ) const {
		return complex_[aEdge].id() ;
	}

	Edge_data& get_data ( edge_descriptor const& aEdge ) const
	{
		return mEdgeDataArray[get_undirected_edge_id(aEdge)];
	}

	Cost_type get_cost(const Profile & profile){
		return (*cost_policy)(profile,get_placement(profile));
	}

	Profile create_profile(edge_descriptor lEdge){
		return Profile(complex_,lEdge);
	}


	void insert_in_PQ( edge_descriptor const& aEdge, Edge_data& aData )
	{
		aData.set_PQ_handle(mPQ->push(aEdge));
	}

	void update_in_PQ( edge_descriptor const& aEdge, Edge_data& aData )
	{
		aData.set_PQ_handle(mPQ->update(aEdge,aData.PQ_handle())) ;
	}

	void remove_from_PQ( edge_descriptor const& aEdge, Edge_data& aData )
	{
		aData.set_PQ_handle(mPQ->erase(aEdge,aData.PQ_handle()));
	}

	boost::optional<edge_descriptor> pop_from_PQ()	{
		boost::optional<edge_descriptor> rEdge = mPQ->extract_top();
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
		size_type lSize = complex_.get_num_edges();
		std::cerr  << "Collecting edges ..."<<std::endl;
		std::cerr  << lSize<<" edges "<<std::endl;

		//mInitialEdgeCount = mCurrentEdgeCount = lSize;

		//mEdgeDataArray.reset( new Edge_data[lSize] ) ;

		mEdgeDataArray.reset( new Edge_data[lSize] ) ;


		mPQ.reset( new PQ (lSize, Compare_cost(this), Undirected_edge_id(this) ) ) ;


		std::size_t id = 0 ;

		EdgeIterator edge_it;
		for(edge_it = complex_.edge_range().begin();
				edge_it != complex_.edge_range().end();
				++edge_it
		){
			edge_descriptor lEdge = *edge_it;
			complex_[lEdge].id() = id++;
			Profile const& lProfile = create_profile(lEdge);
			Edge_data& lData = get_data(lEdge);
			lData.cost() = get_cost(lProfile) ;
			insert_in_PQ(lEdge,lData);
		}
	}

	bool should_stop(double lCost,const Profile &lProfile,int mInitialEdgeCount,int mCurrentEdgeCount){
		return false;
	}

	boost::optional<Point> get_placement(const Profile& lProfile){
		return (*placement_policy)(lProfile);
	}

	bool is_collapse_valid( Profile const& aProfile, Placement_type aPlacement ){
		return (*valid_contraction_policy)(aProfile);
	}


public:
	/**
	 * \brief Contract edges.
	 *
	 * While the heap is not empty, it extracts the edge with the minimum
	 * cost in the heap then try to collapse it.
	 * It stops when the Stop policy says so or when the number of collapses
	 * given by 'num_collapses' is reached.
	 */
	void contract_edges(int num_max_loop){

		DBG("collapse_edges");
		DBGVALUE(complex_.get_num_vertices());
		int num_loop = 0 ;
		//
		// Pops and processes each edge from the PQ
		//
		boost::optional<edge_descriptor> lEdge ;
		while ( (lEdge = pop_from_PQ())&& ((num_loop<num_max_loop)||(num_max_loop<0)))
		{
			++ num_loop;

			DBG("\n\n--------Pop edge");

			Profile const& lProfile = create_profile(*lEdge);
			Cost_type lCost = get_data(*lEdge).cost();
			//Visitor.OnSelected(lProfile,lCost,mInitialEdgeCount,mCurrentEdgeCount);
			if (lCost)
			{
				DBGMSG("lCost",*lCost);
				if (should_stop(*lCost,lProfile,mInitialEdgeCount,mCurrentEdgeCount) )
				{
					//Visitor.OnStopConditionReached(lProfile);
					DBG("should_stop");
					break ;
				}
				Placement_type lPlacement = get_placement(lProfile);
				assert(lPlacement);
				if ( is_collapse_valid(lProfile,lPlacement) )
				{
					DBG("collapse_valid");
					// The external function Get_new_vertex_point() is allowed to return an absent point if there is no way to place the vertex
					// satisfying its constrians. In that case the remaining vertex is simply left unmoved.
					contract_edge(lProfile,lPlacement);
				}
				else
				{
					DBG("collapse not valid");
					//Visitor.OnNonCollapsable(lProfile);
				}
			}
			else
			{
				DBG("uncomputable cost");
			}
		}
	}

	/**
	 * @brief Returns an edge_descriptor which is the edge with the lowest
	 * cost in the heap.
	 * The value is initialized iff the heap is non-empty.
	 */
	boost::optional<edge_descriptor> top_edge(){
		boost::optional<edge_descriptor> res;

		if(!mPQ->empty()) {
			res = mPQ->top();
			DBGMSG("top edge:",complex_[*res]);
		}
		return res;
	}

	Skeleton_blocker_contractor(GeometricSimplifiableComplex& complex)
	:complex_(complex),cost_policy(new Edge_length_cost<Profile>),
	 placement_policy(new Middle_placement<Profile>),
	 valid_contraction_policy(new Link_condition_valid_contraction<Profile>)

	//placement_policy(std::move(placement_policy_)),
	//valid_contraction_policy(std::move(valid_contraction_policy_))
	{
		complex_.set_visitor(this);
		collect_edges();
	}

	Skeleton_blocker_contractor(GeometricSimplifiableComplex& complex,
			std::shared_ptr<Cost_policy<Profile> > &cost_policy_,
			std::shared_ptr<Placement_policy<Profile> >& placement_policy_,
			std::shared_ptr<Valid_contraction_policy<Profile> > &valid_contraction_policy_
	):
		complex_(complex),
		cost_policy(std::move(cost_policy_)),
		placement_policy(std::move(placement_policy_)),
		valid_contraction_policy(std::move(valid_contraction_policy_))

	{
		complex_.set_visitor(this);
		collect_edges();
	}





private:
	void contract_edge( Profile const& aProfile, Placement_type aPlacement ) {
		DBGMSG("collapse edge:",aProfile);
		if (aPlacement)	aProfile.v0().point() = *aPlacement;
		complex_.contract_edge(aProfile.v0_handle(),aProfile.v1_handle());
	}

	/**
	 * @brief we update the cost that has changed and the position in the heap
	 */
	void on_changed_edge(Vertex_handle a,Vertex_handle b){
		//1-get the edge_descriptor corresponding to ab
		//2-change the data in mEdgeArray[ab.id()]
		//3-update the heap
		edge_descriptor lEdge = (complex_[make_pair(a,b)]).first;
		Edge_data& lData = get_data(lEdge);
		Profile const& lProfile = create_profile(lEdge);
		lData.cost() = get_cost(lProfile) ;
		if ( lData.is_in_PQ()){
			update_in_PQ(lEdge,lData);
		}
		else{
			insert_in_PQ(lEdge,lData);
		}
	}


	void on_add_edge(Vertex_handle a,Vertex_handle b){

	}

	void on_remove_edge(Vertex_handle a,Vertex_handle b){

		edge_descriptor lEdge = (complex_[make_pair(a,b)]).first;
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
	void on_swaped_edge(Vertex_handle a,Vertex_handle b,Vertex_handle x){
		std::pair<edge_descriptor,bool> ax_pair = complex_[make_pair(a,x)];
		std::pair<edge_descriptor,bool> bx_pair = complex_[make_pair(b,x)];
		assert(ax_pair.second && bx_pair.second);
		complex_[ax_pair.first].id() =complex_[bx_pair.first].id();
	}




private:

	std::shared_ptr<Cost_policy<Profile> > cost_policy;

	std::shared_ptr<Placement_policy<Profile> > placement_policy;

	std::shared_ptr<Valid_contraction_policy<Profile> > valid_contraction_policy;

	Edge_data_array mEdgeDataArray ;

	boost::scoped_ptr<PQ> mPQ ;

	std::size_t mInitialEdgeCount ;
	std::size_t mCurrentEdgeCount ;

};

}  // namespace contraction


#endif /* SKELETON_BLOCKER_CONTRACTOR_H_ */
