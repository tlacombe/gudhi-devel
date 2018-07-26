#include <iostream>
#include "gudhi/tower_converter.h"
#include "gudhi/hash_complex.h"
#include "gudhi/utilities.h"

using namespace Gudhi::tower_to_filtration;	//module namespace

/**
 * @brief Print usage of example file
 */
void print_usage(){
    std::cout << "Usage:\n";
    std::cout << "  ./example_elementary_input_write_filtration_into_file output_file_name\n";
}

int main(int argc, char *argv[])
{
    if (argc != 2){
	print_usage();
	return 0;
    }

    Tower_converter<Hash_complex> tc(argv[1]);	// by default: output with vertices of simplices, for faces use 'Tower_converter::FACES' as second argument.
						// see documentation for output file format.
    std::vector<double> vertices;

    vertices.push_back(0);
    tc.add_insertion(&vertices, 0);		// add vertex 0 at time 0

    vertices.at(0) = 1;
    tc.add_insertion(&vertices, 1);		// add vertex 1 at time 1

    vertices.at(0) = 2;
    tc.add_insertion(&vertices, 2);		// add vertex 2 at time 2

    vertices.at(0) = 3;
    tc.add_insertion(&vertices, 3);		// add vertex 3 at time 3

    vertices.at(0) = 4;
    tc.add_insertion(&vertices, 4);		// add vertex 4 at time 4

    vertices.at(0) = 5;
    tc.add_insertion(&vertices, 5);		// add vertex 5 at time 5

    vertices.at(0) = 6;
    tc.add_insertion(&vertices, 6);		// add vertex 6 at time 6

    vertices.at(0) = 7;
    tc.add_insertion(&vertices, 7);		// add vertex 7 at time 7

    vertices.at(0) = 8;
    tc.add_insertion(&vertices, 8);		// add vertex 8 at time 8

    vertices.at(0) = 9;
    tc.add_insertion(&vertices, 9);		// add vertex 9 at time 9

    tc.add_contraction(5, 0, 10);		// add contraction of 5 and 0 at time 10, 5 is not valid from now on.

    vertices.at(0) = 10;
    tc.add_insertion(&vertices, 11);		// add vertex 10 at time 11

    vertices.at(0) = 0;
    vertices.push_back(7);
    tc.add_insertion(&vertices, 12);		// add edge (0,7) at time 12

    vertices.at(0) = 0;
    vertices.at(1) = 8;
    tc.add_insertion(&vertices, 13);		// add edge (0,8) at time 13

    vertices.at(0) = 1;
    vertices.at(1) = 3;
    tc.add_insertion(&vertices, 14);		// add edge (1,3) at time 14

    tc.add_contraction(6, 1, 15);		// add contraction of 6 and 1 at time 15, 6 is not valid from now on.

    vertices.pop_back();
    vertices.at(0) = 11;
    tc.add_insertion(&vertices, 16);		// add vertex 11 at time 16

    vertices.at(0) = 1;
    vertices.push_back(7);
    tc.add_insertion(&vertices, 17);		// add edge (1,7) at time 17

    tc.add_contraction(8, 1, 18);		// add contraction of 8 and 1 at time 18, 8 is not valid from now on.

    // ...
}

