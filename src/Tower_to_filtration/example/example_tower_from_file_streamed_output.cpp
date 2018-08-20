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

#include <iostream>
#include "gudhi/tower_converter.h"
#include "gudhi/hash_complex.h"
#include "gudhi/tc_reading_utilities.h"

using namespace Gudhi::tower_to_filtration;	//module namespace

/**
 * @brief Print usage of example file.
 */
void print_usage(){
    std::cout << "Usage:\n";
    std::cout << "  ./example_tower_from_file_streamed_output input_file_name\n";
}

bool get_next_filtration_step(std::stringstream *ss, int *dim, int *filtrationValue, std::vector<double> *vertices);	/**< Reads the filtration output stream. */
void do_something_with_next_filtration_steps(std::stringstream *ss);							/**< Example of how to use the output stream. */

int main(int argc, char *argv[])
{
    if (argc != 2){
	print_usage();
	return 0;
    }

    std::ifstream file(argv[1]);
    std::stringstream ss;
    Tower_converter<Hash_complex> tc(&ss);  // by default: output with vertices of simplices, for faces use 'Tower_converter::FACES' as second argument.
					    // see documentation for stream output format.
    std::string line;

    if (file.is_open()){
	std::vector<double> vertices;
	double timestamp = -1;
	double defaultTimestamp = 0;
	while (getline(file, line, '\n')){
		operationType type = read_operation<Hash_complex>(line, &vertices, &timestamp); // read_operation() function in gudhi/tc_reading_utilities.h, see documentation for file format.
	    if (timestamp != -1) defaultTimestamp = timestamp;

		if (type == INCLUSION){							    // see enumeration in gudhi/tc_reading_utilities.h
		if (tc.add_insertion(vertices, defaultTimestamp)) defaultTimestamp++;	    // add insertion.
		} else if (type == CONTRACTION) {						    // see enumeration in gudhi/tc_reading_utilities.h
		tc.add_contraction(vertices.at(0), vertices.at(1), defaultTimestamp);	    // add contraction.
		defaultTimestamp++;
	    }

	    do_something_with_next_filtration_steps(&ss);				    // use output stream.

	    timestamp = -1;
	}

	file.close();
    } else {
	std::cout << "Unable to open input file.\n";
	file.setstate(std::ios::failbit);
	return 0;
    }
}

void do_something_with_next_filtration_steps(std::stringstream *ss){
    int dim;
    int filtrationValue;
    std::vector<double> vertices;

    while (get_next_filtration_step(ss, &dim, &filtrationValue, &vertices)){	// reads the next filtration operation from output stream
	// do something with result:
	// dim: dimension of new inserted simplex
	// filtrationValue: filtration value of new inserted simplex
	// vertices: vertices of new inserted simplex
	//	in the case where 'Tower_converter::FACES' was used, it contains the faces of the new inserted simplex

	// for example, print information:

	std::cout << dim << " ";
	for (int i = 0; i < vertices.size(); i++){
	    std::cout << vertices.at(i) << " ";
	}
	std::cout << filtrationValue << "\n";
    }

    ss->clear();     //to enable rewriting in ss.
}

bool get_next_filtration_step(std::stringstream *ss, int *dim, int *filtrationValue, std::vector<double> *vertices){
    std::string line;
    if (getline(*ss, line, '\n')){
	std::stringstream nss(line);
	vertices->clear();
	int buf;
	nss >> *dim;
	while (nss >> buf) vertices->push_back(buf);
	*filtrationValue = vertices->back();
	vertices->pop_back();
	return true;
    } else {
	return false;
    }
}

