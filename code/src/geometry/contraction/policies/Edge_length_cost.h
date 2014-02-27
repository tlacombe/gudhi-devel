/*
 * Edge_length_cost.h
 *
 *  Created on: Feb 13, 2014
 *      Author: dsalinas
 */

#ifndef GUDHI_EDGE_LENGTH_COST_H_
#define GUDHI_EDGE_LENGTH_COST_H_

#include "Cost_policy.h"



namespace contraction {

template< typename EdgeProfile> class Edge_length_cost : public Cost_policy<EdgeProfile>{
public:
	typedef typename Cost_policy<EdgeProfile>::Cost_type Cost_type;
	typedef typename EdgeProfile::Point Point;
	Cost_type operator()(const EdgeProfile& profile, const boost::optional<Point>& placement){
		Cost_type res;
		//const Point& a = profile.p0();
		//const Point& b = profile.p1();
		//res = CGAL::squared_distance(a,b);
		// todo length
		return 1.0;
	}

};

}  // namespace contraction

#endif /* GUDHI_EDGE_LENGTH_COST_H_ */
