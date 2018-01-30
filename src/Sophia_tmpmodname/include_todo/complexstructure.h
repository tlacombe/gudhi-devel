/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2018  INRIA Sophia Antipolis-Méditerranée (France)
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
 */

#ifndef COMPLEXSTRUCTURE
#define COMPLEXSTRUCTURE

#include <vector>
#include <string>

namespace Gudhi {
namespace tmp_package_name {

class ComplexStructure
{
public:
    //enum streamType : int {VERTICES, FACES, ACTIVITY};

	virtual ~ComplexStructure() {}
	using vertex = double;
	using index = double;
	using simplex_base = std::vector<vertex>;

    virtual bool insert_simplex(simplex_base *numVertices) = 0;
    virtual vertex get_smallest_closed_star(vertex s1, vertex s2, std::vector<simplex_base*> *closedStar) = 0;  //returns vertex with smallest closed star; closedStar is ordered
                                                                                                                //and simplices in closedStar are independent of the ones in Complex
    virtual bool remove_simplex(simplex_base *numVertices) = 0;
    virtual index get_boundary(simplex_base *simplex, std::vector<index> *boundary) = 0;

    virtual index contract_vertices(vertex v, vertex u, double timestamp, std::vector<std::vector<index>*> *boundaries, std::vector<index> *&inactiveInsertionNumbers) = 0;
	virtual bool exists(simplex_base *simplex) = 0;
	virtual void print(std::string outputFileName) = 0;
    virtual double get_size() const = 0;
    //virtual double get_max_size() const = 0;
    //virtual double get_tower_width() const = 0;
    //virtual double get_filtration_size() const = 0;
    virtual vertex get_vertex_number(double num) const = 0;
};

}
}

#endif // COMPLEXSTRUCTURE

