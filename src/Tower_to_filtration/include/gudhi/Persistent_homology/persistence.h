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

#ifndef PERSISTENCE_H
#define PERSISTENCE_H

#include <vector>
#include <list>
#include <unordered_map>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <cmath>

#include <gudhi/tower_converter.h>

namespace Gudhi {
namespace tower_to_filtration {

template<class ComplexStructure, class ColumnType>
class Persistence
{
public:
    Persistence(double reductionInterval, std::string persistencePairsFileName);
	~Persistence();

    class Boundary_matrix
	{
	public:
        Boundary_matrix(std::string persistencePairsFileName);
        ~Boundary_matrix();

        void insert_column(double insertionNumber, std::vector<double> *boundary, double timestamp);
        void insert_vertex(double insertionNumber, double timestamp);
        void reduce(double start);
        void clear_out();
        void mark_inactive(std::vector<double> *insertionNumbers);
        void mark_inactive(double insertionNumber);

        double get_last_insert_number() const;
        int get_max_dim() const;

	private:
        std::unordered_map<double, ColumnType*> *columns_;
        std::unordered_map<double, double> *latest_;
        std::unordered_map<double, std::pair<bool, bool>*> *isActivePositive_;
        std::unordered_map<double, double> *timestamps_;
        double lastInsertNumber_;
        int maxDim_;
        std::ofstream *persistencePairsFile_;

        void clear_column(double columnIndex);
        void print_persistence_pair(int dim, double birth, double death);
	};

    bool add_insertion(std::vector<double> *simplex, double timestamp);
    void add_contraction(double v, double u, double timestamp);
    void finalize_reduction();

    void print_filtration_data();

private:
    Tower_converter<ComplexStructure> *converter_;
    Boundary_matrix *matrix_;
    double reductionInterval_;
    double lastReduction_;

    void compute_partial_persistence();
};

template<class ComplexStructure, class ColumnType>
Persistence<ComplexStructure, ColumnType>::Persistence(double reductionInterval, std::string persistencePairsFileName) : reductionInterval_(reductionInterval), lastReduction_(-1)
{
    converter_ = new Tower_converter<ComplexStructure>();
    matrix_ = new Boundary_matrix(persistencePairsFileName);
    if (reductionInterval_ < 1) reductionInterval_ = 1;
}

template<class ComplexStructure, class ColumnType>
Persistence<ComplexStructure, ColumnType>::~Persistence()
{
    delete converter_;
    delete matrix_;
}

template<class ComplexStructure, class ColumnType>
bool Persistence<ComplexStructure, ColumnType>::add_insertion(std::vector<double> *simplex, double timestamp)
{
    std::vector<double> boundary;
    double insertionNum;

    if (!converter_->add_insertion(simplex, timestamp, &boundary, &insertionNum)) return false;
    if (simplex->size() == 1) {
        matrix_->insert_vertex(insertionNum, timestamp);
        return true;
    }
    matrix_->insert_column(insertionNum, &boundary, timestamp);

    if (fmod(insertionNum, reductionInterval_) == 0) {
        compute_partial_persistence();
    }

    return true;
}

template<class ComplexStructure, class ColumnType>
void Persistence<ComplexStructure, ColumnType>::add_contraction(double v, double u, double timestamp)
{
    std::vector<std::vector<double>*> boundaries;
    std::vector<double> insertionNumbers;
    bool reduce = false;

    double first = converter_->add_contraction(v, u, timestamp, &boundaries, &insertionNumbers);

    for (std::vector<std::vector<double>*>::size_type i = 0; i < boundaries.size(); i++){
        matrix_->insert_column(first + i, boundaries.at(i), timestamp);
        delete boundaries.at(i);
        if (fmod((first + i), reductionInterval_) == 0) reduce = true;
    }
    matrix_->mark_inactive(&insertionNumbers);

    if (reduce) compute_partial_persistence();
}

template<class ComplexStructure, class ColumnType>
inline void Persistence<ComplexStructure, ColumnType>::finalize_reduction()
{
    if (lastReduction_ != matrix_->get_last_insert_number()) matrix_->reduce(lastReduction_ + 1);
    lastReduction_ = matrix_->get_last_insert_number();
}

template<class ComplexStructure, class ColumnType>
inline void Persistence<ComplexStructure, ColumnType>::print_filtration_data()
{
    converter_->print_filtration_data();
}

template<class ComplexStructure, class ColumnType>
inline void Persistence<ComplexStructure, ColumnType>::compute_partial_persistence()
{
    matrix_->reduce(lastReduction_ + 1);
    matrix_->clear_out();
    lastReduction_ = matrix_->get_last_insert_number();
}

template<class ComplexStructure, class ColumnType>
Persistence<ComplexStructure, ColumnType>::Boundary_matrix::Boundary_matrix(std::string persistencePairsFileName) : lastInsertNumber_(-1), maxDim_(-1)
{
    columns_ = new std::unordered_map<double, ColumnType*>();
    latest_ = new std::unordered_map<double, double>();
    isActivePositive_ = new std::unordered_map<double, std::pair<bool, bool>*>();
    timestamps_ = new std::unordered_map<double, double>();
    persistencePairsFile_ = new std::ofstream(persistencePairsFileName);
    if (!persistencePairsFile_->is_open()){
        std::cout << "Persistence Pairs File could not be open\n";
        exit(0);
    }
}

template<class ComplexStructure, class ColumnType>
Persistence<ComplexStructure, ColumnType>::Boundary_matrix::~Boundary_matrix()
{
    for (auto it = columns_->begin(); it != columns_->end(); it++){
        delete it->second;
    }
    delete columns_;
    delete latest_;
    for (auto it = isActivePositive_->begin(); it != isActivePositive_->end(); it++){
        delete it->second;
    }
    delete isActivePositive_;
    persistencePairsFile_->close();
    delete persistencePairsFile_;
    delete timestamps_;
}

template<class ComplexStructure, class ColumnType>
void Persistence<ComplexStructure, ColumnType>::Boundary_matrix::insert_column(double insertionNumber, std::vector<double> *boundary, double timestamp)
{
    ColumnType *col = new ColumnType(boundary->size() - 1);
    isActivePositive_->emplace(insertionNumber, new std::pair<bool, bool>(true, true));

    for (int i = 0; i < (int)boundary->size(); i++){
        col->push_back(boundary->at(i));
    }

    columns_->emplace(insertionNumber, col);

    lastInsertNumber_ = insertionNumber;
    if (maxDim_ < (int)boundary->size() - 1) maxDim_ = boundary->size() - 1;
    timestamps_->emplace(insertionNumber, timestamp);
}

template<class ComplexStructure, class ColumnType>
inline void Persistence<ComplexStructure, ColumnType>::Boundary_matrix::insert_vertex(double insertionNumber, double timestamp)
{
    isActivePositive_->emplace(insertionNumber, new std::pair<bool, bool>(true, true));
    timestamps_->emplace(insertionNumber, timestamp);
}

template<class ComplexStructure, class ColumnType>
void Persistence<ComplexStructure, ColumnType>::Boundary_matrix::reduce(double start)
{
    for (int d = maxDim_; d > 0; d--){
        for (double i = start; i <= lastInsertNumber_; i++){
            if (columns_->find(i) != columns_->end() && columns_->at(i)->get_dim() == d){
                ColumnType *curr = columns_->at(i);
                double pivot = curr->get_pivot();

                while (pivot != -1 && latest_->find(pivot) != latest_->end()){
                    curr->add(columns_->at(latest_->at(pivot)));
                    pivot = curr->get_pivot();
                }

                if (pivot != -1){
                    isActivePositive_->at(i)->second = false;
                    latest_->emplace(pivot, i);
                    clear_column(pivot);
                    print_persistence_pair(d - 1, pivot, i);
                } else {
                    clear_column(i);
                }
            }
        }
    }
}

template<class ComplexStructure, class ColumnType>
void Persistence<ComplexStructure, ColumnType>::Boundary_matrix::clear_out()
{
    double r;
    double c;

    for (auto it = columns_->begin(); it != columns_->end(); it++){
        c = it->first;
        ColumnType *column = columns_->at(c);
        r = column->get_pivot();
        if (isActivePositive_->at(r)->first){
            column->clean(latest_, isActivePositive_, columns_);
        }
    }

    for (auto it = columns_->begin(), next_it = columns_->begin(); it != columns_->end(); it = next_it)
    {
        next_it = it; ++next_it;
        c = it->first;
        r = columns_->at(c)->get_pivot();
        if (!isActivePositive_->at(r)->first)
        {
            latest_->erase(r);
            delete isActivePositive_->at(r);
            isActivePositive_->erase(r);
            clear_column(c);
            if (!isActivePositive_->at(c)->first) {
                delete isActivePositive_->at(c);
                isActivePositive_->erase(c);
            }
        }
    }

    for (auto it = isActivePositive_->begin(), next_it = isActivePositive_->begin(); it != isActivePositive_->end(); it = next_it)
    {
        next_it = it; ++next_it;
        c = it->first;
        if (columns_->find(c) == columns_->end() && !it->second->second && !it->second->first)
        {
            delete isActivePositive_->at(c);
            isActivePositive_->erase(c);
        }
    }
}

template<class ComplexStructure, class ColumnType>
inline void Persistence<ComplexStructure, ColumnType>::Boundary_matrix::mark_inactive(std::vector<double> *insertionNumbers)
{
    for (std::vector<double>::size_type i = 0; i < insertionNumbers->size(); i++){
        isActivePositive_->at(insertionNumbers->at(i))->first = false;
    }
}

template<class ComplexStructure, class ColumnType>
inline void Persistence<ComplexStructure, ColumnType>::Boundary_matrix::mark_inactive(double insertionNumber)
{
    isActivePositive_->at(insertionNumber)->first = false;
}

template<class ComplexStructure, class ColumnType>
inline double Persistence<ComplexStructure, ColumnType>::Boundary_matrix::get_last_insert_number() const
{
    return lastInsertNumber_;
}

template<class ComplexStructure, class ColumnType>
inline int Persistence<ComplexStructure, ColumnType>::Boundary_matrix::get_max_dim() const
{
    return maxDim_;
}

template<class ComplexStructure, class ColumnType>
inline void Persistence<ComplexStructure, ColumnType>::Boundary_matrix::clear_column(double columnIndex)
{
    if (columns_->find(columnIndex) == columns_->end()) return;
    delete columns_->at(columnIndex);
    columns_->erase(columnIndex);
}

template<class ComplexStructure, class ColumnType>
inline void Persistence<ComplexStructure, ColumnType>::Boundary_matrix::print_persistence_pair(int dim, double birth, double death)
{
    if (timestamps_->at(birth) != timestamps_->at(death)) *persistencePairsFile_ << std::setprecision(std::numeric_limits<double>::digits10 + 1)
                                                                              << dim << " " << timestamps_->at(birth) << " " << timestamps_->at(death) << std::endl;
    timestamps_->erase(birth);
    timestamps_->erase(death);
}

}
}

#endif // PERSISTENCE_H
