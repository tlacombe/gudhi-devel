/*
 * Valid_contraction_policy.h
 *
 *  Created on: Feb 13, 2014
 *      Author: dsalinas
 */

#ifndef VALID_CONTRACTION_POLICY_H_
#define VALID_CONTRACTION_POLICY_H_

namespace contraction {
template< typename EdgeProfile> class Valid_contraction_policy{
public:
	typedef typename EdgeProfile::Point Point;
	typedef typename EdgeProfile::edge_descriptor edge_descriptor;
	typedef typename EdgeProfile::Vertex Vertex;

	virtual bool operator()(const EdgeProfile& profile)=0;
	virtual ~Valid_contraction_policy(){};
};

}  // namespace contraction



#endif /* VALID_CONTRACTION_POLICY_H_ */
