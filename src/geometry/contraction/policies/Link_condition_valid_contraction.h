/*
 * Link_condition_valid_contraction.h
 *
 *  Created on: Feb 13, 2014
 *      Author: dsalinas
 */

#ifndef GUDHI_LINK_CONDITION_VALID_CONTRACTION_H_
#define GUDHI_LINK_CONDITION_VALID_CONTRACTION_H_

#include "utils/Utils.h"


namespace contraction {



template< typename EdgeProfile> class Link_condition_valid_contraction : public Valid_contraction_policy<EdgeProfile>{
public:
	typedef typename EdgeProfile::edge_descriptor edge_descriptor;
	//typedef typename EdgeProfile::edge_descriptor edge_descriptor;
	bool operator()(const EdgeProfile& profile){
		edge_descriptor edge(profile.edge());
		DBGMSG("Link_condition_valid_contraction:",profile.complex().link_condition(edge));
		return profile.complex().link_condition(edge);
	}
};
}  // namespace contraction

#endif /* GUDHI_LINK_CONDITION_VALID_CONTRACTION_H_ */
