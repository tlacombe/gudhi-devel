/*
 * Contraction_visitor.h
 *
 *  Created on: Feb 13, 2014
 *      Author: dsalinas
 */

#ifndef CONTRACTION_VISITOR_H_
#define CONTRACTION_VISITOR_H_


/*
 * ComplexVisitor.h
 *
 *  Created on: Dec 11, 2013
 *      Author: dsalinas
 */

#ifndef COMPLEXVISITOR_H_
#define COMPLEXVISITOR_H_


namespace contraction {

/**
 *@class Contraction_visitor
 *@brief Interface for a visitor of the edge contraction process.
 */
template <typename ComplexType>
class Contraction_visitor {
public:
	virtual ~Contraction_visitor(){};

	virtual void on_add_vertex(Address){};

	/**
	 * @brief Called before the edge contraction process starts.
	 */
	virtual void on_started (ComplexType & complex){}

	/**
	 * @brief Called when the edge contraction process finishes.
	 */
	virtual void on_finished (ComplexType &complex){}


	/**
	 * @brief Called when the StopPredicate returned true (but not if the algorithm terminates because the surface could not be simplified any further).
	 */
	virtual void on_stop_condition_reached (ComplexType &complex){}


	/**
	 * @brief Called during the collecting phase (when a cost is assigned to the edges), for each edge collected.
	 */
 	virtual void on_collected (Profile const &profile, boost::optional< FT > cost){}

 	//Called during the processing phase (when edges are collapsed), for each edge that is selected. More...
 	virtual void on_selected (Profile const &profile, boost::optional< FT > cost, size_type initial_count, size_type current_count){}


 	//Called when an edge is about to be collapsed and replaced by a vertex whose position is *placement. More...
 	virtual void 	on_collapsing(Profile const &profile, boost::optional< Point > placement){}

 	// 	Called for each selected edge which cannot be contracted because the ValidContractionPredicate is false
 	virtual void 	on_non_valid(Profile const &profile){}

};

}  // namespace contraction

#endif /* CONTRACTION_VISITOR_H_ */
