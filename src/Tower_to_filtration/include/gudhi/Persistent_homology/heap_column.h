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

#ifndef HEAPCOLUMN_H
#define HEAPCOLUMN_H

/** @file heap_column.h
 * @brief Contains @ref Gudhi::tower_to_filtration::Heap_column class.
 */

#include <iostream>
#include <algorithm>
#include <cmath>
#include <vector>
#include <unordered_map>

namespace Gudhi {
namespace tower_to_filtration {

/**
 * @brief Column type which is based on `std::vector<double>` and `std::make_heap`. Fulfills the requirements of the @ref Gudhi::tower_to_filtration::ColumnType concept.
 * Its coefficient are in @f$\mathbb{Z}_2@f$ only. Therefore the values of the vector are the indices of non-zero cells in increasing order.
 */
class Heap_column
{
public:
    Heap_column(int dim);
    ~Heap_column();

public:
    void add(Heap_column *columnToAdd);
    void set(std::vector<double> *newColumn);
    double get_pivot();
    double pop_pivot();
    double get_size();
    /**
     * @brief Returns the stored dimension.
     * @return The dimension.
     */
    int get_dim() const { return dim_; }
    void clean(std::unordered_map<double, double> *latest, std::unordered_map<double, std::pair<bool, bool> *> *isActivePositive,
               std::unordered_map<double, Heap_column*> *columns);
    void push_back(double value);
    double at(double index);

private:
    int dim_;                           /**< Dimension of the column. */
    std::vector<double> *column_;       /**< Data container of the column. */
    std::vector<double> *temp_column_;  /**< Temporary data container of the column. */
    double insertsSinceLastPrune_;      /**< Number of insertion since the last time the heap was pruned. */

    void prune();
};

/**
 * @brief Constructor
 * @param dim dimension to be stored as the column dimension.
 */
inline Heap_column::Heap_column(int dim) : dim_(dim), insertsSinceLastPrune_(0)
{
    column_ = new std::vector<double>();
    std::make_heap(this->column_->begin(), this->column_->end());
    temp_column_ = new std::vector<double>();
}

/**
 * @brief Destructor
 */
inline Heap_column::~Heap_column()
{
    delete column_;
    delete temp_column_;
}

/**
 * @brief Replaces the column values by the sum of this column and @p columnToAdd.
 *
 * Adding two columns can create duplication of values in the column, which makes the computation faster but more space consumming.
 * This function also calls a pruning function to get rid of duplicates whenever the size of the column exceets two times the number of insertion since last pruning.
 *
 * @param columnToAdd column to sum with.
 */
inline void Heap_column::add(Heap_column *columnToAdd)
{
    double size = columnToAdd->get_size();
    for (double i = 0; i < size; i++) {
        column_->push_back(columnToAdd->at(i));
        std::push_heap(column_->begin(), column_->end());
    }
    insertsSinceLastPrune_ += size;

    if (2 * insertsSinceLastPrune_ > (double)column_->size()) prune();
}

/**
 * @brief Replaces the content of the column with the content of @p newColumn.
 * @param newColumn new column content. Its values represent the index of non-zero cells in increasing order (@f$\mathbb{Z}_2@f$ coefficients).
 */
inline void Heap_column::set(std::vector<double> *newColumn)
{
    column_->clear();
    column_->insert(column_->end(), newColumn->begin(), newColumn->end());
    std::make_heap(column_->begin(), column_->end());
}

/**
 * @brief Insert a cell at the end of the column.
 * @param value value of the cell to be inserted.
 */
inline void Heap_column::push_back(double value)
{
    column_->push_back(value);
    std::push_heap(column_->begin(), column_->end());
}

/**
 * @brief Returns element at index @p index in the underlying vector.
 *
 * The vector being used at a heap, the indices do NOT correspond to the order of the element in the column!
 * Usefull for traversal of the column, when the order does not matter.
 *
 * @param index desired index.
 * @return Element at index @p index in the underlying vector.
 */
inline double Heap_column::at(double index)
{
    return column_->at(index);
}

/**
 * @brief Returns the pivot of the column, i.e. the index of the last nonzero value.
 * @return The pivot of the column.
 */
inline double Heap_column::get_pivot()
{
    double pivot = pop_pivot();
    if (pivot != -1){
        column_->push_back(pivot);
        std::push_heap(column_->begin(), column_->end());
    }
    return pivot;
}

/**
 * @brief Returns and removes the pivot of the column (i.e. the index of the last nonzero value) from it.
 *
 * It will also remove every double representation of the pivot in the heap.
 *
 * @return The pivot of the column.
 */
inline double Heap_column::pop_pivot()
{
    if (column_->empty()) {
        return -1;
    } else {
        double pivot = column_->front();
        std::pop_heap(column_->begin(), column_->end());
        column_->pop_back();
        while (!column_->empty() && column_->front() == pivot) {
            std::pop_heap(column_->begin(), column_->end());
            column_->pop_back();
            if (column_->empty()) {
                return -1;
            } else {
                pivot = column_->front();
                std::pop_heap(column_->begin(), column_->end());
                column_->pop_back();
            }
        }
        return pivot;
    }
}

/**
 * @brief Returns the number of nonzero values in the column.
 * @return The number of nonzero values in the column.
 */
inline double Heap_column::get_size()
{
    prune();
    return column_->size();
}

/**
 * @brief Erase useless cells from the column. See @cite KerberS17 for more details.
 * @param latest private member of @ref Gudhi::tower_to_filtration::Persistence::Boundary_matrix.
 * @param isActivePositive private member of @ref Gudhi::tower_to_filtration::Persistence::Boundary_matrix.
 * @param columns private member of @ref Gudhi::tower_to_filtration::Persistence::Boundary_matrix.
 */
inline void Heap_column::clean(std::unordered_map<double, double> *latest, std::unordered_map<double, std::pair<bool, bool> *> *isActivePositive,
                        std::unordered_map<double, Heap_column*> *columns)
{
    std::vector<double> *tmp = new std::vector<double>();
    tmp->push_back(pop_pivot());
    double max = pop_pivot();
    while (max != -1){
        if (latest->find(max) != latest->end() && !isActivePositive->at(max)->first){
            push_back(max);
            add(columns->at(latest->at(max)));
        } else if (isActivePositive->at(max)->second || isActivePositive->at(max)->first) {
            tmp->push_back(max);
        }
        max = pop_pivot();
    }
    std::reverse(tmp->begin(), tmp->end());
    delete column_;
    column_ = tmp;
    std::make_heap(column_->begin(), column_->end());
}

/**
 * @brief Prunes the heap representation of the column, i.e. removes duplications which appears when adding two columns.
 */
inline void Heap_column::prune()
{
    if (insertsSinceLastPrune_ == 0) return;

    std::vector<double> *tempCol = temp_column_;
    tempCol->clear();
    double pivot = pop_pivot();
    while (pivot != -1) {
        tempCol->push_back(pivot);
        pivot = pop_pivot();
    }
    temp_column_ = column_;
    column_ = tempCol;
    std::reverse(column_->begin(), column_->end());
    std::make_heap(column_->begin(), column_->end());

    insertsSinceLastPrune_ = 0;
}

}
}

#endif // HEAPCOLUMN_H
