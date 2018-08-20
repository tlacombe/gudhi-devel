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

#ifndef CONCEPT_COMPLEX_STRUCTURE_H_
#define CONCEPT_COMPLEX_STRUCTURE_H_

/** @file ComplexStructure.h
 * @brief Contains @ref Gudhi::tower_to_filtration::ComplexStructure concept.
 */

namespace Gudhi {
namespace tower_to_filtration {

/**
 * @brief Data structure for a complex. This concept describes the methodes required by the module.
 */
class ComplexStructure{
public:
    /**
     * @brief Constructor with no parameters
     */
    ComplexStructure();

    /**
     * @brief Insert a simplex in the complex, whose facets were already inserted in the complex.
     * @param numVertices simplex to be inserted, described as the vector of its vertex identifiers in increasing order.
     * @return true if the simplex was not already inserted in the complex, false otherwise.
     */
    bool insert_simplex(std::vector<double> *numVertices);
    /**
     * @brief Remove a simplex and all its cofaces from the complex.
     * @param vertex simplex to be removed. This modules only needs to call the function on vertices.
     * @return true if removal was successful, false otherwise.
     */
    bool remove_simplex(std::vector<double> *vertex);
    /**
     * @brief Remove a simplex and all its cofaces from the complex.
     * @param vertex simplex to be removed. This modules only needs to call the function on vertices.
     * @param removedIndices pointer to an (empty) vector of doubles ; if this parameter is given, the identifiers of all the removed simplices are stored in the vector.
     * @return true if removal was successful, false otherwise.
     */
    bool remove_simplex(std::vector<double> *vertex, std::vector<double> *removedIndices);
    /**
     * @brief Compute the smallest closed star of the vertices v and u and returns the corresponding vertex.
     *
     * The star of a simplex is the set of simplices containing this simplex. The closed star is the smallest subcomplex containing the star.
     *
     * @param v identifier of the first vertex
     * @param u identifier of the second vertex
     * @param closedStar pointer to an empty vector of simplices ; the method stores here the simplices contained in the smallest closed star, ordered such that possible faces
     * of a simplex have smaller indices than the simplex itself.
     * The simplices are represented by pointers to vectors containing the vertex indentifiers of the simplex in increasing order.
     * Simplices will be modified and deleted later on and therefore should not be the same than the one stored in the complex.
     * @return The identifier of the vertex which had the smallest closed star.
     */
    double get_smallest_closed_star(double v, double u, std::vector<std::vector<double>*> *closedStar);
    /**
     * @brief Computes the boundary of a simplex and returns the identifier of the simplex.
     * @param simplex simplex which boundary will be computed, described as the vector of its vertex identifiers in increasing order.
     * @param boundary pointer to an (empty) vector of doubles ; the method stores here the simplices in the boundary in increasing order of insertion.
     * @return The identifier of the simplex.
     */
    double get_boundary(std::vector<double> *simplex, std::vector<double> *boundary);
    /**
     * @brief Returns the current size if the complex.
     * @return Current size of the complex.
     */
    double get_size() const;
    /**
     * @brief Returns the biggest identifier of a simplex currently in the complex.
     * @return Current biggest identifier of a simplex.
     */
    double get_max_index() const;
    /**
     * @brief Returns the maximal size the complex had at some point until now.
     * @return The maximal size the complex had at some point until now.
     */
    double get_max_size() const;
    /**
     * @brief Returns the maximal dimension of the simplices currently in the complex.
     * @return Maximal dimension of the simplices currently in the complex.
     */
    int get_max_dimension() const;
};

}
}

#endif  // CONCEPT_COMPLEX_STRUCTURE_H_
