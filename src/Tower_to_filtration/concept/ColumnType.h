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

#ifndef CONCEPT_COLUMN_TYPE_H_
#define CONCEPT_COLUMN_TYPE_H_

/** @file ColumnType.h
 * @brief Contains @ref Gudhi::tower_to_filtration::ColumnType concept.
 */

namespace Gudhi {
namespace tower_to_filtration {

/**
 * @brief Data structure for the columns of the boundary matrix used by @ref Gudhi::tower_to_filtration::Persistence<ComplexStructure,ColumnType>.
 * This concept describes the methodes required by the module.
 */
class ColumnType{
public:
    /**
     * @brief Constructor
     * @param dim dimension of the simplex whose boundary is represented by this column.
     */
    ColumnType(int dim);

    /**
     * @brief Insert a cell at the end of the column.
     * @param cell value of the cell to be inserted.
     */
    void push_back(double cell);
    /**
     * @brief Return the pivot of the column, which is the last nonzero value.
     * @return The last nonzero value of the column.
     */
    double get_pivot();
    /**
     * @brief Sum this column and another column. The result replaces this column.
     * @param columnToAdd column to add up.
     */
    void add(ColumnType *columnToAdd);
    /**
     * @brief Erase useless cells from the column. See @cite KerberS17 for more details.
     * @param latest private member of @ref Gudhi::tower_to_filtration::Persistence<ComplexStructure,ColumnType>::Boundary_matrix.
     * @param isActivePositive private member of @ref Gudhi::tower_to_filtration::Persistence<ComplexStructure,ColumnType>::Boundary_matrix.
     * @param columns private member of @ref Gudhi::tower_to_filtration::Persistence<ComplexStructure,ColumnType>::Boundary_matrix.
     */
    void clean(std::unordered_map<double, double> *latest, std::unordered_map<double, std::pair<bool, bool> *> *isActivePositive,
               std::unordered_map<double, List_column *> *columns);
};

}
}

#endif  // CONCEPT_COLUMN_TYPE_H_
