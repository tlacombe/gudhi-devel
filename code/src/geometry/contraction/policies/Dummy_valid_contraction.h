/*
 * Dummy_valid_contraction.h
 *
 *  Created on: Feb 13, 2014
 *      Author: dsalinas
 */

#ifndef GUDHI_DUMMY_VALID_CONTRACTION_H_
#define GUDHI_DUMMY_VALID_CONTRACTION_H_

namespace contraction {




template< typename EdgeProfile> class Dummy_valid_contraction : public Valid_contraction_policy<EdgeProfile>{
public:
	bool operator()(const EdgeProfile& profile){
		return true;
	}
};

}  // namespace contraction



#endif /* GUDHI_DUMMY_VALID_CONTRACTION_H_ */
