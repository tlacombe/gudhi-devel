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

//#include "hash_complex.h"
#include <list_column.h>
#include <heap_column.h>

#define GUDHI_LIST

namespace Gudhi {
namespace tmp_package_name {

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
		void reduce(double start);
        void clear_out();
        void mark_inactive(std::vector<double> *insertionNumbers);
        void mark_inactive(double insertionNumber);
        void insert_vertex(double insertionNumber, double timestamp);

        double get_last_insert_number() const;
        int get_max_dim() const;

	private:
#ifdef GUDHI_LIST
		typedef std::list<double> ContainerType;
		typedef List_column ColumnType;
#else
		typedef std::vector<double> ContainerType;
		typedef HeapColumn ColumnType;
#endif

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

    bool conv_insert_simplex(std::vector<double> *simplex, double timestamp);
    void conv_contract_vertices(double u, double v, double timestamp);
    void insert_simplex(double insertionNum, std::vector<double> *boundary, double timestamp);
    void finalize_reduction();
    void mark_inactive(double insertionNumber);

    void print_complex_data();

private:
    //ComplexStructure *complex_;
    Boundary_matrix *matrix_;
    double reductionInterval_;
    double lastReduction_;

    void compute_partial_persistence();
};

}
}

#endif // PERSISTENCE_H
