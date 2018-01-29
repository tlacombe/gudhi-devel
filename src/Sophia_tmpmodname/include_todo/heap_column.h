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

#include <iostream>
#include <algorithm>
#include <cmath>
#include <vector>

namespace Gudhi {
namespace tmp_package_name {

class Heap_column
{
public:
    Heap_column(std::vector<double> *column, int dim);
    ~Heap_column();

public:
    void add(Heap_column *columnToAdd);
    int get_dim() const;
    double get_pivot();
    double pop_pivot();
    double get_size();
	double at(double index);
	void set(std::vector<double> *newColumn);
	void insert(double value);

	void print();

private:
    int dim_;
    std::vector<double> *column_;
    std::vector<double> *temp_column_;
    double insertsSinceLastPrune_;

	void prune();
};

}
}

#endif // HEAPCOLUMN_H
