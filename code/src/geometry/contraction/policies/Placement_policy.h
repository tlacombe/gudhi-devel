/*
 * Placement_policy.h
 *
 *  Created on: Feb 13, 2014
 *      Author: dsalinas
 */

#ifndef PLACEMENT_POLICY_H_
#define PLACEMENT_POLICY_H_

#include <boost/optional.hpp>

namespace contraction {

template< typename EdgeProfile> class Placement_policy{
public:
	typedef typename EdgeProfile::Point Point;
	typedef boost::optional<Point> Placement_type;

	virtual Placement_type operator()(const EdgeProfile& profile)=0;
	virtual ~Placement_policy(){};
};


}  // namespace contraction


#endif /* PLACEMENT_POLICY_H_ */
