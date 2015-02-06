/*
 * Projection_from_points.h
 *  Created on: Feb 6, 2015
 * This file is part of the Gudhi Library. The Gudhi library 
 *    (Geometric Understanding in Higher Dimensions) is a generic C++ 
 *    library for computational topology.
 *
 *    Author(s):       David Salinas
 *
 *    Copyright (C) 2014  INRIA Sophia Antipolis-Méditerranée (France)
 *
 *    This program is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *
 *    This program is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 * 
 */


#ifndef PROJECTION_FROM_POINTS_H_
#define PROJECTION_FROM_POINTS_H_

#include <iostream>
#include "utils/UI_utils.h"
#include "Projector3D.h"


class Projection_from_points : public Projector3D{
	typedef Projector3D::Point Point;
	typedef Projector3D::Point_3 Point_3;

	std::vector<Point_3> projected_points_;

public:

	Projection_from_points(const std::vector<Point_3>& projected_points):
		projected_points_(projected_points){
	}

	Point_3 operator()(Vertex_handle v) const{
		assert(p.dimension()>=3);
		return projected_points_[v];
	}

	void print(){
		for(auto p : projected_points_)
			std::cout << p << std::endl;
	}
};





#endif /* PROJECTION_FROM_POINTS_H_ */
