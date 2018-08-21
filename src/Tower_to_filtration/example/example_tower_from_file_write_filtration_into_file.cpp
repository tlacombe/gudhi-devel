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

#include <iostream>
#include "gudhi/tower_converter.h"
#include "gudhi/hash_complex.h"
#include "gudhi/tc_reading_utilities.h"

using namespace Gudhi::tower_to_filtration;	//module namespace

/**
 * @brief Print usage of example file
 */
void print_usage(){
    std::cout << "Usage:\n";
    std::cout << "  ./example_tower_from_file_write_filtration_into_file input_file_name output_file_name\n";
}

int main(int argc, char *argv[])
{
    if (argc != 3){
	print_usage();
	return 0;
    }

    std::ifstream file(argv[1]);
    Tower_converter<Hash_complex> tc(argv[2]);	    // by default: output with vertices of simplices, for faces use 'Tower_converter::FACES' as second argument.
						    // see documentation for output file format.

    if (file.is_open()){
	file >> tc;				    // >> function in gudhi/tc_reading_utilities.h, see documentation for input file format.
	file.close();
    } else {
	std::cout << "Unable to open input file\n";
	file.setstate(std::ios::failbit);
	return 0;
    }
}

