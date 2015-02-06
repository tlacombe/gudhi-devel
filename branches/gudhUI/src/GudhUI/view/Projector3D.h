/*
 * Projector.h
 *
 *  Created on: Aug 27, 2014
 *      Author: dsalinas
 */

#ifndef PROJECTOR3D_H_
#define PROJECTOR3D_H_

#include "model/Complex_typedefs.h"


class Projector3D{
public:
	typedef Geometry_trait::Point Point;
	typedef Geometry_trait::Point_3 Point_3;
	typedef typename Complex::Vertex_handle Vertex_handle;
	virtual Point_3 operator()(Vertex_handle) const  = 0;

	virtual ~Projector3D(){
	}

};



#endif /* PROJECTOR3D_H_ */
