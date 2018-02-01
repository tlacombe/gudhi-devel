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

public:
    void add(List_column *columnToAdd);
    int get_dim() const;
    std::list<double>::iterator get_begin_iterator();
    std::list<double>::reverse_iterator get_reverse_begin_iterator();
    std::list<double>::iterator get_end_iterator();
    std::list<double>::reverse_iterator get_reverse_end_iterator();
	void erase(std::list<double>::iterator *pos);
    double get_pivot();
    double get_size();

private:
    int dim_;
    std::list<double> *column_;
};

}
}

#endif // LISTCOLUMN_H
