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
    std::cout << "  ./example_elementary_input_streamed_output\n";
}

bool get_next_filtration_step(std::stringstream *ss, int *dim, int *filtrationValue, std::vector<double> *vertices);	/**< Reads the filtration output stream. */
void do_something_with_next_filtration_steps(std::stringstream *ss);							/**< Example of how to use the output stream. */

int main(int argc, char *argv[])
{
    if (argc != 1){
	print_usage();
	return 0;
    }

    std::stringstream ss;
    Tower_converter<Hash_complex> tc(&ss);	    // by default: output with vertices of simplices, for faces use 'Tower_converter::FACES' as second argument.
						    // see documentation for stream output format.
    std::vector<double> vertices;

    vertices.push_back(0);
    tc.add_insertion(&vertices, 0);		    // add vertex 0 at time 0
    do_something_with_next_filtration_steps(&ss);   // use output stream

    vertices.at(0) = 1;
    tc.add_insertion(&vertices, 1);		    // add vertex 1 at time 1
    do_something_with_next_filtration_steps(&ss);   // use output stream

    vertices.at(0) = 2;
    tc.add_insertion(&vertices, 2);		    // add vertex 2 at time 2
    do_something_with_next_filtration_steps(&ss);   // use output stream

    vertices.at(0) = 3;
    tc.add_insertion(&vertices, 3);		    // add vertex 3 at time 3
    do_something_with_next_filtration_steps(&ss);   // use output stream

    vertices.at(0) = 4;
    tc.add_insertion(&vertices, 4);		    // add vertex 4 at time 4
    do_something_with_next_filtration_steps(&ss);   // use output stream

    vertices.at(0) = 5;
    tc.add_insertion(&vertices, 5);		    // add vertex 5 at time 5
    do_something_with_next_filtration_steps(&ss);   // use output stream

    vertices.at(0) = 6;
    tc.add_insertion(&vertices, 6);		    // add vertex 6 at time 6
    do_something_with_next_filtration_steps(&ss);   // use output stream

    vertices.at(0) = 7;
    tc.add_insertion(&vertices, 7);		    // add vertex 7 at time 7
    do_something_with_next_filtration_steps(&ss);   // use output stream

    vertices.at(0) = 8;
    tc.add_insertion(&vertices, 8);		    // add vertex 8 at time 8
    do_something_with_next_filtration_steps(&ss);   // use output stream

    vertices.at(0) = 9;
    tc.add_insertion(&vertices, 9);		    // add vertex 9 at time 9
    do_something_with_next_filtration_steps(&ss);   // use output stream

    tc.add_contraction(5, 0, 10);		    // add contraction of 5 and 0 at time 10, 5 is not valid from now on.
    do_something_with_next_filtration_steps(&ss);   // use output stream

    vertices.at(0) = 10;
    tc.add_insertion(&vertices, 11);		    // add vertex 10 at time 11
    do_something_with_next_filtration_steps(&ss);   // use output stream

    vertices.at(0) = 0;
    vertices.push_back(7);
    tc.add_insertion(&vertices, 12);		    // add edge (0,7) at time 12
    do_something_with_next_filtration_steps(&ss);   // use output stream

    vertices.at(0) = 0;
    vertices.at(1) = 8;
    tc.add_insertion(&vertices, 13);		    // add edge (0,8) at time 13
    do_something_with_next_filtration_steps(&ss);   // use output stream

    vertices.at(0) = 1;
    vertices.at(1) = 3;
    tc.add_insertion(&vertices, 14);		    // add edge (1,3) at time 14
    do_something_with_next_filtration_steps(&ss);   // use output stream

    tc.add_contraction(6, 1, 15);		    // add contraction of 6 and 1 at time 15, 6 is not valid from now on.
    do_something_with_next_filtration_steps(&ss);   // use output stream

    vertices.pop_back();
    vertices.at(0) = 11;
    tc.add_insertion(&vertices, 16);		    // add vertex 11 at time 16
    do_something_with_next_filtration_steps(&ss);   // use output stream

    vertices.at(0) = 1;
    vertices.push_back(7);
    tc.add_insertion(&vertices, 17);		    // add edge (1,7) at time 17
    do_something_with_next_filtration_steps(&ss);   // use output stream

    tc.add_contraction(8, 1, 18);		    // add contraction of 8 and 1 at time 18, 8 is not valid from now on.
    do_something_with_next_filtration_steps(&ss);   // use output stream

    // ...
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

