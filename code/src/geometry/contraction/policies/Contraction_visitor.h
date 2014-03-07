/*
 * Contraction_visitor.h
 *
 *  Created on: Feb 13, 2014
 *      Author: dsalinas
 */

#ifndef GUDHI_CONTRACTION_VISITOR_H_
#define GUDHI_CONTRACTION_VISITOR_H_

#include "geometry/contraction/Edge_profile.h"
#include "boost/optional.hpp"

namespace contraction {

/**
 *@class Contraction_visitor
 *@brief Interface for a visitor of the edge contraction process.
 */
template <typename ComplexType>
class Contraction_visitor {
public:
	//typedef typename ComplexType::GeometryTrait GT;
	typedef Edge_profile<ComplexType> Profile;
	typedef double FT;
	typedef typename ComplexType::Point Point;


	virtual ~Contraction_visitor(){};

	/**
	 * @brief Called before the edge contraction process starts.
	 */
	virtual void on_started (ComplexType & complex){}

	/**
	 * @brief Called when the StopPredicate returned true (but not if the algorithm terminates because the surface could not be simplified any further).
	 */
	virtual void on_stop_condition_reached (const Profile &profile){}


	/**
	 * @brief Called during the collecting phase (when a cost is assigned to the edges), for each edge collected.
	 */
	virtual void on_collected (const Profile &profile, boost::optional< FT > cost){}

	/**
	 * @brief Called during the processing phase (when edges are collapsed), for each edge that is selected.
	 */
	virtual void on_selected (const Profile &profile, boost::optional< FT > cost, int initial_count, int current_count){}


	/**
	 * @brief Called when an edge is about to be contracted and replaced by a vertex whose position is *placement.
	 */
	virtual void on_contracting(const Profile &profile, boost::optional< Point > placement){}

	/**
	 * @brief Called when after an edge has been contracted onto a new point placement.
	 * A possibility would to remove popable blockers at this point for instance.
	 */
	virtual void on_contracted(const Profile  &profile, boost::optional< Point > placement){}


	/**
	 * @brief Called for each selected edge which cannot be contracted because the ValidContractionPredicate is false
	 */
	virtual void on_non_valid(const Profile  &profile){}

};

}  // namespace contraction

#endif /* GUDHI_CONTRACTION_VISITOR_H_ */
