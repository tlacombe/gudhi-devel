/*
 * Cost_policy.h
 *
 *  Created on: Feb 13, 2014
 *      Author: dsalinas
 */

#ifndef COST_POLICY_H_
#define COST_POLICY_H_

#include <boost/optional.hpp>


namespace contraction {

template< typename EdgeProfile> class Cost_policy{
public:
	typedef typename EdgeProfile::Point Point;
	typedef typename EdgeProfile::Vertex Vertex;

	typedef boost::optional<double> Cost_type;

	virtual Cost_type operator()(const EdgeProfile& profile, const boost::optional<Point>& placement)=0;
	virtual ~Cost_policy(){};
};

}  // namespace contraction
#endif /* COST_POLICY_H_ */
