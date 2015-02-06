/*
 * FirstCoordProjector.h
 *
 *  Created on: Aug 27, 2014
 *      Author: dsalinas
 */

#ifndef FIRSTCOORDPROJECTOR_H_
#define FIRSTCOORDPROJECTOR_H_

#include "utils/UI_utils.h"
#include "Projector3D.h"

template<typename Complex>
class FirstCoordProjector3D : public Projector3D{
	typedef Projector3D::Point Point;
	typedef Projector3D::Point_3 Point_3;
	const Complex& complex_;

public:

	FirstCoordProjector3D(const Complex& complex):complex_(complex){}

	Point_3 operator()(Vertex_handle v) const{
//		assert(p.dimension()>=3);
		const auto & p = complex_.point(v);
		return Point_3(p.x(),p.y(),p.z());
	}
};

#endif /* FIRSTCOORDPROJECTOR_H_ */
