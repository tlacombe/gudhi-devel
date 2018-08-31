/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2018  TU Graz (Austria)
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
    typename vertex;	    /**< Type for vertex identifiers. */
    typename index;	    /**< Type for simplex identifiers. */
    typename size_type;	    /**< Type for size mesure. */

    /**
     * @brief Constructor with no parameters
     */
    ComplexStructure();

    /**
     * @brief Insert a simplex in the complex, whose facets were already inserted in the complex.
     * @param simplex simplex to be inserted, described as the vector of its vertex identifiers. The module will give them in increasing order.
     * @return true if the simplex was not already inserted in the complex, false otherwise.
     *	    This module uses the result in an `if` statement, the output does not have to be `bool` if the test result is the same.
     */
    bool insert_simplex(std::vector<vertex> &simplex);
    /**
     * @brief Inserts the edge @p uv (and its vertices if not inserted) and all its possible cofaces.
     * @param u vertex identifier of the first vertex of the edge ; the module calls the method such that @p u < @p v.
     * @param v vertex identifier of the second vertex of the edge ; the module calls the method such that @p u < @p v.
     * @param maxDim the maximal dimension of the cofaces to be inserted ; if -1, then there is no limit.
     * @param addedSimplices pointer to an (empty) vector of simplices ; the method stores in insertion order all newly inserted simplices here.
     * @param addedIndices pointer to an (empty) vector of identifiers ;
     *	    the method stores there the identifiers of the simplices in @p addedSimplices in the corresponding order.
     * @param boundaries pointer to an (empty) vector of boundary identifiers ;
     *	    the method stores there the identifiers of the boundary simplices of the simplices in @p addedSimplices in the corresponding order.
     * @return true if the edge was not already inserted in the complex, false otherwise.
     */
    bool insert_edge_and_expand(vertex u, vertex v, int maxDim, std::vector<std::vector<vertex> > *addedSimplices, std::vector<index> *addedIndices, std::vector<std::vector<index>*> *boundaries);
    /**
     * @brief Remove a simplex and all its cofaces from the complex.
     * @param vertex simplex to be removed. This modules only needs to call the function on vertices.
     * @param removedIndices pointer to an (empty) vector of simplex identifiers ; if this parameter is given, the identifiers of all the removed simplices are stored in the vector.
     */
    void remove_simplex(std::vector<vertex> &vertex, std::vector<index> *removedIndices);
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
    vertex get_smallest_closed_star(vertex v, vertex u, std::vector<std::vector<vertex>*> *closedStar);
    /**
     * @brief Computes the boundary of a simplex and returns the identifier of the simplex.
     * @param simplex simplex which boundary will be computed, described as the vector of its vertex identifiers in increasing order.
     * @param boundary pointer to an (empty) vector of simplex identifiers ; the method stores here the simplices in the boundary in increasing order of insertion.
     * @return The identifier of the simplex.
     */
    index get_boundary(std::vector<vertex> &simplex, std::vector<index> *boundary);
    /**
     * @brief Returns the current size if the complex.
     * @return Current size of the complex.
     */
    size_type get_size() const;
    /**
     * @brief Returns the biggest identifier of a simplex currently in the complex.
     * @return Current biggest identifier of a simplex.
     */
    index get_max_index() const;
};

}
}

#endif  // CONCEPT_COMPLEX_STRUCTURE_H_
