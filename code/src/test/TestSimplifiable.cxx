/*
 * TestSimplifiable.cxx
 *
 *  Created on: Feb 4, 2014
 *      Author: dsalinas
 */




#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <fstream>
#include <sstream>
#include "Test.h"
//#include "Skeleton_blocker/Simplex.h"
#include "Skeleton_blocker/Skeleton_blocker_complex.h"
#include "Skeleton_blocker/Skeleton_blocker_complex_iterators.h"
#include "Skeleton_blocker/Skeleton_blocker_simplifiable_complex.h"
#include "Skeleton_blocker/Skeleton_blocker_simple_traits.h"


using namespace std;


template<typename ComplexType> class Skeleton_blocker_sub_complex;
typedef Skeleton_blocker_simplifiable_complex<Skeleton_blocker_simple_traits> Complex;
typedef typename Complex::Vertex_handle Vertex_handle;
typedef typename Complex::Root_vertex_handle Root_vertex_handle;
typedef Skeleton_blocker_simplex<Vertex_handle> Simplex_handle;
// true iff v \in complex
bool assert_vertex(Complex &complex,Vertex_handle v){
	Simplex_handle simplex(v);
	bool test = complex.contains(simplex);
	assert(test);
	return test;
}

// true iff the blocker (a,b,c) is in complex
bool assert_blocker(Complex &complex,Root_vertex_handle a,Root_vertex_handle b,Root_vertex_handle c){

	return complex.contains_blocker(Simplex_handle(*complex.get_address(a),*complex.get_address(b),*complex.get_address(c)));
	//return complex.contains_blocker((a),(b),(c));
}

void build_complete(int n,Complex& complex){
	complex.clear();
	for(int i=0;i<n;i++)
		complex.add_vertex();
	for(int i=0;i<n;i++)
		for(int j=0;j<i;j++)
			complex.add_edge(Vertex_handle(i),Vertex_handle(j));
}

bool test_contraction1(){
	enum { a, b, x, y, z, n };
	Complex complex(n);
	build_complete(n,complex);
	complex.remove_edge(Vertex_handle(b),Vertex_handle(z));
	complex.add_blocker(Vertex_handle(a),Vertex_handle(x),Vertex_handle(y));
	complex.add_blocker(Vertex_handle(b),Vertex_handle(x),Vertex_handle(y));

	// Print result
	cerr << "complex before complex"<< complex.to_string()<<endl;

	cerr <<endl<<endl;
	complex.contract_edge(Vertex_handle(a),Vertex_handle(b));
	// Print result
	cerr << "ContractEdge(0,1)\n";
	PRINT(complex.to_string());

	// verification
	for (int i=0;i<5;i++)
		if (i!=1) assert_vertex(complex,Vertex_handle(i));
	bool test1 = !complex.contains_edge(Vertex_handle(a),Vertex_handle(b));
	bool test2 = assert_blocker(complex,Root_vertex_handle(a),Root_vertex_handle(x),Root_vertex_handle(y));
	bool test3 = complex.num_edges()==6;
	bool test4 = complex.num_blockers()==1;
	Simplex_handle sigma;
	sigma.add_vertex(Vertex_handle(a));
	sigma.add_vertex(Vertex_handle(x));
	sigma.add_vertex(Vertex_handle(y));
	sigma.add_vertex(Vertex_handle(z));
	bool test5 = !(complex.contains(sigma));
	return test1&&test2&&test3&&test4&&test5;
}


bool test_contraction2(){
	// Make convenient labels for the vertices
	enum { a, b, x, y, z, n };
	Complex complex(n);
	build_complete(n,complex);
	complex.remove_edge(Vertex_handle(b),Vertex_handle(x));
	Simplex_handle blocker;
	blocker.add_vertex(Vertex_handle(a));
	blocker.add_vertex(Vertex_handle(y));
	blocker.add_vertex(Vertex_handle(z));

	complex.add_blocker(blocker);

	// Print result
	cerr << "complex complex"<< complex.to_string();
	cerr <<endl<<endl;
	complex.contract_edge(Vertex_handle(a),Vertex_handle(b));

	cerr << "complex.ContractEdge(a,b)"<< complex.to_string();

	cerr <<endl<<endl;

	// there should be one blocker (a,c,d,e) in the complex
	bool test ;
	test = complex.contains_blocker(Simplex_handle(a,x,y,z));
	test = test && complex.num_blockers()==1;
	return test;
}




bool test_collapse1(){
	// xxx implement remove_star(simplex) before
	/*
	Complex complex(5);
	buildComplete(4,complex);
	complex.add_vertex();
	complex.add_edge(2,4);
	complex.add_edge(3,4);
	// Print result
	cerr << "complex complex"<< complex.to_string();
	cerr <<endl<<endl;

	Simplex *sigma= new Simplex(1,2,3);
	complex.remove_star(*sigma);
	cerr << "complex.remove_star(1,2,3)"<< complex.to_string();
	cerr <<endl<<endl;

	// verification
	assert_blocker(complex,1,2,3);
	cerr <<"----> Ocomplex \n";
	*/
	return true;

}

bool test_collapse2(){
	// Make convenient labels for the vertices
	Complex complex(5);
	build_complete(4,complex);
	complex.add_vertex();
	complex.add_edge(Vertex_handle(1),Vertex_handle(4));
	complex.add_edge(Vertex_handle(2),Vertex_handle(4));
	complex.add_edge(Vertex_handle(3),Vertex_handle(4));
	complex.add_blocker(Vertex_handle(1),Vertex_handle(2),Vertex_handle(3),Vertex_handle(4));
	// Print result
	cerr << "complex complex"<< complex.to_string();
	cerr <<endl<<endl;

	Simplex_handle *sigma= new Simplex_handle(Vertex_handle(1),Vertex_handle(2),Vertex_handle(3));
	complex.remove_star(*sigma);
	cerr << "complex.remove_star(1,2,3)"<< complex.to_string();
	cerr <<endl<<endl;

	// verification
	assert_blocker(complex,Root_vertex_handle(1),Root_vertex_handle(2),Root_vertex_handle(3));
	assert(!complex.contains_blocker(Simplex_handle(Vertex_handle(1),Vertex_handle(2),Vertex_handle(3),Vertex_handle(4))));
	cerr <<"----> Ocomplex \n";
	return true;
}

bool test_collapse3(){
	// Make convenient labels for the vertices
	Complex complex(5);
	build_complete(4,complex);
	complex.add_vertex();
	complex.add_edge(Vertex_handle(1),Vertex_handle(4));
	complex.add_edge(Vertex_handle(2),Vertex_handle(4));
	complex.add_edge(Vertex_handle(3),Vertex_handle(4));
	complex.add_blocker(Vertex_handle(1),Vertex_handle(2),Vertex_handle(3),Vertex_handle(4));
	// Print result
	cerr << "complex complex"<< complex.to_string();
	cerr <<endl<<endl;

	complex.remove_star(Vertex_handle(2));
	cerr << "complex complex"<< complex.to_string();
	// verification
	//	assert_blocker(complex,1,2,3);
	//	assert(!complex.ContainsBlocker(new AddressSimplex(1,2,3,4)));
	cerr <<"----> Ocomplex \n";
	return true;
}

void test_collapse4(){
	// Make convenient labels for the vertices
	Complex complex(6);
	build_complete(4,complex);
	complex.add_vertex();
	complex.add_vertex();
	complex.add_edge(Vertex_handle(1),Vertex_handle(4));
	complex.add_edge(Vertex_handle(2),Vertex_handle(4));
	complex.add_edge(Vertex_handle(3),Vertex_handle(4));
	complex.add_edge(Vertex_handle(3),Vertex_handle(4));
	complex.add_edge(Vertex_handle(2),Vertex_handle(5));
	complex.add_edge(Vertex_handle(3),Vertex_handle(5));
	complex.add_edge(Vertex_handle(4),Vertex_handle(5));
	// Print result
	cerr << "complex complex"<< complex.to_string();
	cerr <<endl<<endl;
	complex.remove_star(Vertex_handle(1));
	cerr << "complex complex"<< complex.to_string();
	// verification
	//	assert_blocker(complex,1,2,3);
	//	assert(!complex.ContainsBlocker(new AddressSimplex(1,2,3,4)));
	cerr <<"----> Ocomplex \n";
}


bool test_remove_popable_blockers(){
	Complex complex;
	build_complete(4,complex);
	complex.add_vertex();
	complex.add_edge(Vertex_handle(3),Vertex_handle(4));
	complex.add_edge(Vertex_handle(2),Vertex_handle(4));
	Simplex_handle sigma1=Simplex_handle(Vertex_handle(1),Vertex_handle(2),Vertex_handle(3));
	Simplex_handle sigma2=Simplex_handle(Vertex_handle(2),Vertex_handle(3),Vertex_handle(4));

	complex.add_blocker(sigma1);
	complex.add_blocker(sigma2);
	cerr << "complex complex"<< complex.to_string();
	cerr <<endl<<endl;
	cerr << "complex.RemovePopableBlockers();"<<endl;
	complex.remove_popable_blockers();
	cerr << "complex complex"<< complex.to_string();
	cerr <<endl<<endl;

	bool test1 = (complex.num_blockers()==1);


	// test 2
	complex.clear();
	build_complete(4,complex);
	complex.add_vertex();
	complex.add_vertex();
	complex.add_edge(Vertex_handle(3),Vertex_handle(4));
	complex.add_edge(Vertex_handle(2),Vertex_handle(4));
	complex.add_edge(Vertex_handle(2),Vertex_handle(5));
	complex.add_edge(Vertex_handle(3),Vertex_handle(5));
	complex.add_edge(Vertex_handle(4),Vertex_handle(5));
	sigma1=Simplex_handle(Vertex_handle(1),Vertex_handle(2),Vertex_handle(3));
	sigma2=Simplex_handle(Vertex_handle(2),Vertex_handle(3),Vertex_handle(4));

	complex.add_blocker(sigma1);
	complex.add_blocker(sigma2);
	cerr << "complex complex"<< complex.to_string();
	cerr <<endl<<endl;	cerr << "complex.RemovePopableBlockers();"<<endl;
	complex.remove_popable_blockers();
	cerr << "complex complex"<< complex.to_string();

	cerr <<endl<<endl;
	bool test2 = (complex.num_blockers()==0);
	return test1&&test2;
}



int main (int argc, char *argv[])
{
	Tests tests_simplifiable_complex;
	tests_simplifiable_complex.add("Test contraction 1",test_contraction1);
	tests_simplifiable_complex.add("Test contraction 2",test_contraction2);
	tests_simplifiable_complex.add("Test remove popable blockers",test_remove_popable_blockers);


	tests_simplifiable_complex.run();

	if(tests_simplifiable_complex.run())
		return EXIT_SUCCESS;
	else
		return EXIT_FAILURE;
}
