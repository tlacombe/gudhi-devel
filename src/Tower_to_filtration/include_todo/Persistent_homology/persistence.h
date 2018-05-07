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

#include <tower_converter.h>
#include <Persistent_homology/list_column.h>
#include <Persistent_homology/heap_column.h>

#define GUDHI_LIST

namespace Gudhi {
namespace tower_to_filtration {

template<class ComplexStructure>
class Persistence
{
public:
    Persistence(double reductionInterval, std::string persistencePairsFileName);
	~Persistence();

    template<class ColumnType, typename ContainerType>
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
#ifdef GUDHI_LIST
    Boundary_matrix<List_column, std::list<double> > *matrix_;
#else
    Boundary_matrix<Heap_column, std::vector<double> > *matrix_;
#endif
    double reductionInterval_;
    double lastReduction_;

    void compute_partial_persistence();
};

template<class ComplexStructure>
Persistence<ComplexStructure>::Persistence(double reductionInterval, std::string persistencePairsFileName) : reductionInterval_(reductionInterval), lastReduction_(-1)
{
    converter_ = new Tower_converter<ComplexStructure>();
#ifdef GUDHI_LIST
    matrix_ = new Boundary_matrix<List_column, std::list<double> >(persistencePairsFileName);
#else
    matrix_ = new Boundary_matrix<Heap_column, std::vector<double> >(persistencePairsFileName);
#endif
    if (reductionInterval_ < 1) reductionInterval_ = 1;
}

template<class ComplexStructure>
Persistence<ComplexStructure>::~Persistence()
{
    delete converter_;
    delete matrix_;
}

template<class ComplexStructure>
bool Persistence<ComplexStructure>::add_insertion(std::vector<double> *simplex, double timestamp)
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

template<class ComplexStructure>
void Persistence<ComplexStructure>::add_contraction(double v, double u, double timestamp)
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

template<class ComplexStructure>
inline void Persistence<ComplexStructure>::finalize_reduction()
{
    if (lastReduction_ != matrix_->get_last_insert_number()) matrix_->reduce(lastReduction_ + 1);
    lastReduction_ = matrix_->get_last_insert_number();
}

template<class ComplexStructure>
inline void Persistence<ComplexStructure>::print_filtration_data()
{
    converter_->print_filtration_data();
}

template<class ComplexStructure>
inline void Persistence<ComplexStructure>::compute_partial_persistence()
{
    matrix_->reduce(lastReduction_ + 1);
    matrix_->clear_out();
    lastReduction_ = matrix_->get_last_insert_number();
}

template<class ComplexStructure>
template<class ColumnType, typename ContainerType>
Persistence<ComplexStructure>::Boundary_matrix<ColumnType, ContainerType>::Boundary_matrix(std::string persistencePairsFileName) : lastInsertNumber_(-1), maxDim_(-1)
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

template<class ComplexStructure>
template<class ColumnType, typename ContainerType>
Persistence<ComplexStructure>::Boundary_matrix<ColumnType, ContainerType>::~Boundary_matrix()
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

template<class ComplexStructure>
template<class ColumnType, typename ContainerType>
void Persistence<ComplexStructure>::Boundary_matrix<ColumnType, ContainerType>::insert_column(double insertionNumber, std::vector<double> *boundary, double timestamp)
{
    ContainerType *boundaryCells = new ContainerType();
    isActivePositive_->emplace(insertionNumber, new std::pair<bool, bool>(true, true));

    for (int i = 0; i < (int)boundary->size(); i++){
        boundaryCells->push_back(boundary->at(i));
    }

    columns_->emplace(insertionNumber, new ColumnType(boundaryCells, boundary->size() - 1));

    lastInsertNumber_ = insertionNumber;
    if (maxDim_ < (int)boundary->size() - 1) maxDim_ = boundary->size() - 1;
    timestamps_->emplace(insertionNumber, timestamp);
}

template<class ComplexStructure>
template<class ColumnType, typename ContainerType>
inline void Persistence<ComplexStructure>::Boundary_matrix<ColumnType, ContainerType>::insert_vertex(double insertionNumber, double timestamp)
{
    isActivePositive_->emplace(insertionNumber, new std::pair<bool, bool>(true, true));
    timestamps_->emplace(insertionNumber, timestamp);
}

template<class ComplexStructure>
template<class ColumnType, typename ContainerType>
void Persistence<ComplexStructure>::Boundary_matrix<ColumnType, ContainerType>::reduce(double start)
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

template<class ComplexStructure>
template<class ColumnType, typename ContainerType>
void Persistence<ComplexStructure>::Boundary_matrix<ColumnType, ContainerType>::clear_out()
{
    double r;
#ifndef GUDHI_LIST
    ContainerType tmp;
#endif

    double c;
    for (auto it = columns_->begin(); it != columns_->end(); it++){
        c = it->first;
        ColumnType *column = columns_->at(c);
        r = column->get_pivot();
        if (isActivePositive_->at(r)->first){
#ifdef GUDHI_LIST
            typename ContainerType::reverse_iterator it;
            typename ContainerType::iterator it2;
            it = column->get_reverse_begin_iterator();
            it++;
            while (it != column->get_reverse_end_iterator()){
                if (latest_->find(*it) != latest_->end() && !isActivePositive_->at(*it)->first){
                    column->add(columns_->at(latest_->at(*(it--))));
                } else if (!isActivePositive_->at(*it)->second && !isActivePositive_->at(*it)->first) {
                    it2 = (++it).base();
                    it--; it--;
                    column->erase(&it2);
                }
                it++;
            }
#else
            tmp.clear();
            tmp.push_back(column->pop_pivot());
            double max = column->pop_pivot();
            while (max != -1){
                if (latest_->find(max) != latest_->end() && !isActivePositive_->at(max)->first){
                    column->insert(max);
                    column->add(columns_->at(latest->at(max)));
                } else if (isActivePositive_->at(max)->second || isActivePositive_->at(max)->first) {
                    tmp.push_back(max);
                }
                max = column->pop_pivot();
            }
            std::reverse(tmp.begin(), tmp.end());
            column->set(&tmp);
#endif
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

template<class ComplexStructure>
template<class ColumnType, typename ContainerType>
inline void Persistence<ComplexStructure>::Boundary_matrix<ColumnType, ContainerType>::mark_inactive(std::vector<double> *insertionNumbers)
{
    for (std::vector<double>::size_type i = 0; i < insertionNumbers->size(); i++){
        isActivePositive_->at(insertionNumbers->at(i))->first = false;
    }
}

template<class ComplexStructure>
template<class ColumnType, typename ContainerType>
inline void Persistence<ComplexStructure>::Boundary_matrix<ColumnType, ContainerType>::mark_inactive(double insertionNumber)
{
    isActivePositive_->at(insertionNumber)->first = false;
}

template<class ComplexStructure>
template<class ColumnType, typename ContainerType>
inline double Persistence<ComplexStructure>::Boundary_matrix<ColumnType, ContainerType>::get_last_insert_number() const
{
    return lastInsertNumber_;
}

template<class ComplexStructure>
template<class ColumnType, typename ContainerType>
inline int Persistence<ComplexStructure>::Boundary_matrix<ColumnType, ContainerType>::get_max_dim() const
{
    return maxDim_;
}

template<class ComplexStructure>
template<class ColumnType, typename ContainerType>
inline void Persistence<ComplexStructure>::Boundary_matrix<ColumnType, ContainerType>::clear_column(double columnIndex)
{
    if (columns_->find(columnIndex) == columns_->end()) return;
    delete columns_->at(columnIndex);
    columns_->erase(columnIndex);
}

template<class ComplexStructure>
template<class ColumnType, typename ContainerType>
inline void Persistence<ComplexStructure>::Boundary_matrix<ColumnType, ContainerType>::print_persistence_pair(int dim, double birth, double death)
{
    if (timestamps_->at(birth) != timestamps_->at(death)) *persistencePairsFile_ << std::setprecision(std::numeric_limits<double>::digits10 + 1)
                                                                              << dim << " " << timestamps_->at(birth) << " " << timestamps_->at(death) << std::endl;
    timestamps_->erase(birth);
    timestamps_->erase(death);
}

}
}

#endif // PERSISTENCE_H
