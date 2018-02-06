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

#ifndef LISTCOLUMN_H
#define LISTCOLUMN_H

#include <iostream>
#include <list>

namespace Gudhi {
namespace tmp_package_name {

class List_column
{
public:
    List_column(std::list<double> *column, int dim);
    ~List_column();

    void add(List_column *columnToAdd);
    void erase(std::list<double>::iterator *pos){ column_->erase(*pos); }
    std::list<double>::iterator get_begin_iterator(){ return column_->begin(); }
    std::list<double>::reverse_iterator get_reverse_begin_iterator(){ return column_->rbegin(); }
    std::list<double>::iterator get_end_iterator(){ return column_->end(); }
    std::list<double>::reverse_iterator get_reverse_end_iterator(){ return column_->rend(); }
    double get_size(){ return column_->size(); }
    int get_dim() const{ return dim_; }
    double get_pivot();

private:
    int dim_;
    std::list<double> *column_;
};

List_column::List_column(std::list<double> *column, int dim) : dim_(dim), column_(column)
{}

List_column::~List_column()
{
    delete column_;
}

void List_column::add(List_column *columnToAdd)
{
    std::list<double>::iterator itToAdd = columnToAdd->get_begin_iterator(), itTarget = column_->begin();
    while (itToAdd != columnToAdd->get_end_iterator() && itTarget != column_->end()){
        if (*itToAdd == *itTarget){
            column_->erase(itTarget++);
            itToAdd++;
        } else if (*itToAdd < *itTarget){
            column_->insert(itTarget, *itToAdd);
            itToAdd++;
        } else {
            itTarget++;
        }
    }
    while (itToAdd != columnToAdd->get_end_iterator()){
        column_->push_back(*itToAdd);
        itToAdd++;
    }
}

inline double List_column::get_pivot(){
    if (column_->empty()) return -1;
    return column_->back();
}

}
}

#endif // LISTCOLUMN_H
