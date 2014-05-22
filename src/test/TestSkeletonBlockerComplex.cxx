#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <fstream>
#include <sstream>
#include "Utils.h"

#include "Test.h"
//#include "Skeleton_blocker/Simplex.h"
#include "Skeleton_blocker/Skeleton_blocker_complex.h"
#include "Skeleton_blocker/Skeleton_blocker_complex_iterators.h"
#include "Skeleton_blocker/Skeleton_blocker_link_complex.h"
#include "Skeleton_blocker/Skeleton_blocker_link_superior.h"
#include "Skeleton_blocker/Skeleton_blocker_simple_traits.h"
//#include "Simple_vertex.h"
//#include "Simple_edge.h"

using namespace std;

template<typename ComplexType> class Skeleton_blocker_link_complex;
typedef Skeleton_blocker_complex<Skeleton_blocker_simple_traits> Complex;
typedef typename Complex::Vertex_handle Vertex_handle;
typedef typename Complex::Root_vertex_handle Root_vertex_handle;
typedef typename Complex::Simplex_handle Simplex_handle;
typedef typename Complex::Root_simplex_handle Root_simplex_handle;
typedef Simplex_handle::Simplex_vertex_const_iterator Simplex_vertex_const_iterator;
typedef Complex::Const_complex_blocker_around_vertex_iterator Const_blocker_iterator;
typedef Complex::Edge_handle Edge_handle;

// true iff v \in complex
bool assert_vertex(Complex &complex,Vertex_handle v){
	//assert(complex.contains(v));
	return complex.contains(v);
}

bool assert_simplex(Complex &complex,Root_vertex_handle a,Root_vertex_handle b,Root_vertex_handle c){
	return true;
//	AddressSimplex simplex((a),(b),(c));
//	return complex.contains(&simplex);
}

// true iff the blocker (a,b,c) is in complex
bool assert_blocker(Complex &complex,Root_vertex_handle a,Root_vertex_handle b,Root_vertex_handle c){
	return true;
	//return complex.contains_blocker((a),(b),(c));
}

// true iff the blocker (a,b,c,d) is in complex
bool assert_blocker(Complex &complex,Root_vertex_handle a,Root_vertex_handle b,Root_vertex_handle c,Root_vertex_handle d){
	return true;
	//Simplex blocker (a,b,c,d);
	//return complex.contains_blocker(&blocker);
}


void build_complete(int n,Complex& complex){
	complex.clear();
	for(int i=0;i<n;i++)
		complex.add_vertex();
	for(int i=0;i<n;i++)
		for(int j=0;j<i;j++)
			complex.add_edge(Vertex_handle(i),Vertex_handle(j));
}


bool test_simplex(){
	PRINT("test simplex");
	Simplex_handle simplex(Vertex_handle(0),Vertex_handle(1),Vertex_handle(2),Vertex_handle(3));
	for (auto i = simplex.begin() ; i != simplex.end() ; ++i){
		PRINT(*i);
		auto j = i;
		for (++j ;
				j != simplex.end() ;
				++j){
			PRINT(*j);
		}
	}
	return simplex.dimension()==3;
}


bool test_iterator_vertices(){
	int n = 10;
	Complex complex(10);
	cerr << "complex.num_vertices():"<<complex.num_vertices()<<endl;
	int num_vertex_seen = 0;
	for(auto vi :complex.vertex_range()){
		cerr << "vertex:"<<vi<<endl;
		++num_vertex_seen;
	}
	return num_vertex_seen == n;
}


bool test_iterator_edge(){
	const int n = 10;
	Complex complex(n);
	for(int i=0;i<n;i++)
		for(int j=0;j<i;j++)
			complex.add_edge(Vertex_handle(i),Vertex_handle(j));
	complex.remove_edge(Vertex_handle(2),Vertex_handle(3));
	complex.remove_edge(Vertex_handle(3),Vertex_handle(5));
	cerr << "complex.num_edges():"<<complex.num_edges()<<endl;
	int num_edges_seen = 0;
	for(auto edge : complex.edge_range()){
		cerr << "edge :"<<complex[edge]<<endl;
		++num_edges_seen;
	}

	return num_edges_seen == n*(n-1)/2-2;
}

bool test_iterator_edge2(){
	const int n = 10;
	Complex complex(n);
	for(int i=0;i<n;i++)
		for(int j=0;j<i;j++)
			complex.add_edge(Vertex_handle(i),Vertex_handle(j));
	complex.remove_edge(Vertex_handle(2),Vertex_handle(3));
	complex.remove_edge(Vertex_handle(3),Vertex_handle(5));
	cerr << "complex.num_edges():"<<complex.num_edges()<<endl;
	int num_neigbors_seen = 0;
	for(auto neighbor : complex.vertex_range(Vertex_handle(2))){
		cerr << "neighbor"<<neighbor<<endl;
		++num_neigbors_seen;
	}
	return num_neigbors_seen==8;
}



bool test_iterator_edge3(){
	const int n = 10;
	Complex complex(n);
	for(int i=0;i<n;i++)
		for(int j=0;j<i;j++)
			complex.add_edge(Vertex_handle(i),Vertex_handle(j));
	complex.remove_edge(Vertex_handle(2),Vertex_handle(3));
	complex.remove_edge(Vertex_handle(3),Vertex_handle(5));
	cerr << "complex.num_edges():"<<complex.num_edges()<<endl;
	int num_neigbors_seen = 0;
	for(auto edge : complex.edge_range(Vertex_handle(2))){
		std::cerr << edge<< std::endl;
		++num_neigbors_seen;
	}
	return num_neigbors_seen==8;
}



bool test_iterator_triangles(){
	const int n = 7;
	Complex complex(n);
	//create a "ring" around '0'
	for(int i=1;i<n;i++)
		complex.add_edge(Vertex_handle(0),Vertex_handle(i));
	for(int i=1;i<n-1;i++)
		complex.add_edge(Vertex_handle(i),Vertex_handle(i+1));
	complex.add_edge(Vertex_handle(1),Vertex_handle(6));

	PRINT(complex.to_string());

	int num_triangles_seen=0;
	//for (auto t : complex.triangle_range(5)){
	TEST("triangles around 5 (should be 2 of them):");
	for (auto t : complex.triangle_range(5)){
		PRINT(t);
		++num_triangles_seen;
	}
	bool test = (num_triangles_seen==2);

	num_triangles_seen=0;
	TEST("triangles around 0 (should be 6 of them):");
	for (auto t : complex.triangle_range(0)){
		PRINT(t);
		++num_triangles_seen;
	}
	test = test&&(num_triangles_seen==6);

	TEST("superior triangles around 0 (should be 6 of them):");
	num_triangles_seen=0;
	for (auto t : complex.superior_triangle_range(0)){
		PRINT(t);
		++num_triangles_seen;
	}
	test = test&&(num_triangles_seen==6);

	// bug
	TEST("superior triangles around 1 (should be 0 of them):");
	num_triangles_seen=0;
	for (auto t : complex.superior_triangle_range(1)){
		PRINT(t);
		++num_triangles_seen;
	}
	test = test&&(num_triangles_seen==0);

	// we now add another triangle
	complex.add_vertex();
	complex.add_edge(Vertex_handle(4),Vertex_handle(7));
	complex.add_edge(Vertex_handle(3),Vertex_handle(7));
	complex.add_blocker(Vertex_handle(0),Vertex_handle(1),Vertex_handle(6));
	num_triangles_seen=0;
	TEST("superior triangles around 3 (should be 1 of them):");
	for (auto t : complex.superior_triangle_range(3)){
		PRINT(t);
		++num_triangles_seen;
	}
	test = test&&(num_triangles_seen==1);

	TEST("triangles (should be 6 of them):");
	num_triangles_seen=0;
	for (auto t : complex.triangle_range()){
		PRINT(t);
		++num_triangles_seen;
	}
	test = test&&(num_triangles_seen==6);
	PRINT(num_triangles_seen);

	return test;

}

bool test_iterator_blockers(){
	Complex complex;
	Simplex_handle alpha;
	Simplex_handle vertex_set_expected;
	// Build the complexes
	for (int i=0;i<20;i++){
		complex.add_vertex();
	}
	for (int i=10;i<15;i++){
		for (int j=i+1;j<15;j++)
			complex.add_edge(Vertex_handle(i),Vertex_handle(j));
	}

	complex.add_blocker(Simplex_handle(Vertex_handle(10),Vertex_handle(11),Vertex_handle(12)));
	complex.add_blocker(Simplex_handle(Vertex_handle(2),Vertex_handle(1),Vertex_handle(10)));
	complex.add_blocker(Simplex_handle(Vertex_handle(10),Vertex_handle(9),Vertex_handle(15)));
	complex.add_blocker(Simplex_handle(Vertex_handle(1),Vertex_handle(9),Vertex_handle(8)));

	// Print result
	int num_blockers=0;
	for(auto blockers : complex.blocker_range(Vertex_handle(10))){
		TESTVALUE(*blockers) ;
		num_blockers++;
	}
	bool test = (num_blockers==3);

	num_blockers=0;
	for (auto blockers : complex.blocker_range()){
		TESTVALUE(*blockers) ;
		num_blockers++;
	}
	test = test && (num_blockers==4) ;
	return test;
}


bool test_link0(){

	enum { a, b, c, d, n };
	Complex complex(n);
	complex.add_edge(Vertex_handle(b),Vertex_handle(c));complex.add_edge(Vertex_handle(c),Vertex_handle(d));
	Simplex_handle alpha = Simplex_handle(Vertex_handle(c));
	Skeleton_blocker_link_complex<Complex> L(complex,alpha);
	PRINT(L.num_vertices());
	PRINT(L.to_string());

	bool test1 = L.contains_vertex(*L.get_address(Root_vertex_handle(b)));
	bool test2 = L.contains_vertex(*L.get_address(Root_vertex_handle(d)));
	bool test3 = L.num_edges()==0;
	bool test4 = L.num_blockers()==0;
	return test1&&test2&&test3&&test4;

}

bool test_link1(){
	Complex complex;


	// Build the complexes
	for (int i=0;i<20;i++){
		complex.add_vertex();
	}
	for (int i=10;i<15;i++){
		for (int j=i+1;j<15;j++)
			complex.add_edge(Vertex_handle(i),Vertex_handle(j));
	}
	Simplex_handle alpha(Vertex_handle(12),Vertex_handle(14));
	Skeleton_blocker_link_complex<Complex> L(complex,alpha);
	// Complexes built

	// verification
	bool test1 = L.contains_vertex(*L.get_address(Root_vertex_handle(10)));
	bool test2 = L.contains_vertex(*L.get_address(Root_vertex_handle(11)));
	bool test3 = L.contains_vertex(*L.get_address(Root_vertex_handle(13)));
	bool test4 = L.num_edges()==3;
	bool test5 = L.num_blockers()==0;
	Root_simplex_handle simplex;
	simplex.add_vertex(Root_vertex_handle(10));
	simplex.add_vertex(Root_vertex_handle(11));
	simplex.add_vertex(Root_vertex_handle(13));
	bool test6 = L.get_simplex_address(simplex);
	bool test7 = L.contains(*(L.get_simplex_address(simplex)));
	cerr <<"----> Ocomplex \n";
	return test1&&test2&&test3&&test4&&test5&&test6&&test7 ;

}


bool test_link2(){
	Complex complex;

	Simplex_handle alpha;
	Simplex_handle vertex_set_expected;
	// Build the complexes
	for (int i=0;i<20;i++){
		complex.add_vertex();
	}
	for (int i=10;i<15;i++){
		for (int j=i+1;j<15;j++)
			complex.add_edge(Vertex_handle(i),Vertex_handle(j));
	}
	complex.add_blocker(Vertex_handle(10),Vertex_handle(11),Vertex_handle(13));
	alpha = Simplex_handle(Vertex_handle(12),Vertex_handle(14));
	Skeleton_blocker_link_complex<Complex> L(complex,alpha);
	// Complexes built

	// Print result
	cerr << "complex complex"<< complex.to_string();
	cerr <<endl<<endl;
	cerr << "L= Link_complex("<<alpha<<") : \n"<<L.to_string();

	// verification
	bool test1 = L.contains_vertex(*L.get_address(Root_vertex_handle(10)));
	bool test2 = L.contains_vertex(*L.get_address(Root_vertex_handle(11)));
	bool test3 = L.contains_vertex(*L.get_address(Root_vertex_handle(13)));
	bool test4 = L.num_edges()==3;
	bool test5 = L.num_blockers()==1;
	Root_simplex_handle simplex;
	simplex.add_vertex(Root_vertex_handle(10));
	simplex.add_vertex(Root_vertex_handle(11));
	simplex.add_vertex(Root_vertex_handle(13));
	bool test6 = L.contains_blocker(*(L.get_simplex_address(simplex)));
	cerr <<"----> Ocomplex \n";
	return test1&&test2&&test3&&test4&&test5&&test6 ;
}

bool test_link3(){
	Complex complex;

	Simplex_handle alpha;
	Simplex_handle vertex_set_expected;
	// Build the complexes
	for (int i=0;i<20;i++){
		complex.add_vertex();
	}
	for (int i=10;i<15;i++){
		for (int j=i+1;j<15;j++)
			complex.add_edge(Vertex_handle(i),Vertex_handle(j));
	}
	complex.add_blocker(Vertex_handle(10),Vertex_handle(11),Vertex_handle(12));
	alpha = Simplex_handle(Vertex_handle(12),Vertex_handle(14));
	Skeleton_blocker_link_complex<Complex> L(complex,alpha);
	// Complexes built

	// Print result
	cerr << "complex complex"<< complex.to_string();
	cerr <<endl<<endl;
	cerr << "L= Link_complex("<<alpha<<") : \n"<<L.to_string();

	// verification
	bool test = assert_vertex(L,*L.get_address(Root_vertex_handle(10)));
	test = test&& assert_vertex(L,*L.get_address(Root_vertex_handle(11)));
	test = test&& assert_vertex(L,*L.get_address(Root_vertex_handle(13)));
	test = test&& L.num_edges()==2;
	test = test&&L.contains_edge(*L.get_address(Root_vertex_handle(10)),*L.get_address(Root_vertex_handle(13)));
	test=test&&L.contains_edge(*L.get_address(Root_vertex_handle(13)),*L.get_address(Root_vertex_handle(11)));
	test=test&&L.num_blockers()==0;
	return test;
}

bool test_link4(){
	Complex complex;

	// Build the complexes
	for (int i=0;i<20;i++){
		complex.add_vertex();
	}
	for (int i=10;i<15;i++){
		for (int j=i+1;j<15;j++)
			complex.add_edge(Vertex_handle(i),Vertex_handle(j));
	}
	complex.add_blocker(Vertex_handle(10),Vertex_handle(11),Vertex_handle(12),Vertex_handle(13));
	Simplex_handle alpha(Vertex_handle(12),Vertex_handle(14));
	Skeleton_blocker_link_complex<Complex> L(complex,alpha);
	// Complexes built

	// verification
	bool test1 = L.contains_vertex(*L.get_address(Root_vertex_handle(10)));
	bool test2 = L.contains_vertex(*L.get_address(Root_vertex_handle(11)));
	bool test3 = L.contains_vertex(*L.get_address(Root_vertex_handle(13)));
	bool test4 = L.num_edges()==3;
	bool test5 = L.num_blockers()==1;
	Root_simplex_handle simplex;
	simplex.add_vertex(Root_vertex_handle(10));
	simplex.add_vertex(Root_vertex_handle(11));
	simplex.add_vertex(Root_vertex_handle(13));
	bool test6 = L.contains_blocker(*(L.get_simplex_address(simplex)));
	cerr <<"----> Ocomplex \n";
	return test1&&test2&&test3&&test4&&test5&&test6 ;

}

bool test_link5(){
	Complex complex(0,new Print_complex_visitor<Vertex_handle>());
	// Build the complexes
	build_complete(4,complex);
	complex.add_blocker(Vertex_handle(0),Vertex_handle(1),Vertex_handle(2),Vertex_handle(3));

	Simplex_handle alpha(Vertex_handle(0),Vertex_handle(1),Vertex_handle(2));


	Skeleton_blocker_link_complex<Complex> L(complex,alpha);	// Complexes built

	// Print result
	PRINT(complex.to_string());
	cerr <<endl<<endl;
	PRINT(L.to_string());

	// verification
	return L.num_vertices()==0;
}




template<typename SimplexHandle>
void add_triangle_edges(int a,int b,int c,list<SimplexHandle>& simplices){
	typedef SimplexHandle Simplex_handle;
	typedef typename SimplexHandle::Vertex_handle Vertex_handle;

	simplices.push_back(Simplex_handle(Vertex_handle(a),Vertex_handle(b) ));
	simplices.push_back(Simplex_handle(Vertex_handle(b),Vertex_handle(c) ));
	simplices.push_back(Simplex_handle(Vertex_handle(c),Vertex_handle(a) ));
}

template<typename SimplexHandle>
void add_triangle(int a,int b,int c,list<SimplexHandle>& simplices){
	typedef SimplexHandle Simplex_handle;
	typedef typename SimplexHandle::Vertex_handle Vertex_handle;
	simplices.push_back(Simplex_handle(Vertex_handle(a),Vertex_handle(b),Vertex_handle(c)));
}

bool test_constructor(){
	list <Simplex_handle> simplices;

	simplices.push_back(Simplex_handle(Vertex_handle(0)));
	simplices.push_back(Simplex_handle(Vertex_handle(1)));
	simplices.push_back(Simplex_handle(Vertex_handle(2)));
	simplices.push_back(Simplex_handle(Vertex_handle(3)));
	simplices.push_back(Simplex_handle(Vertex_handle(4)));
	simplices.push_back(Simplex_handle(Vertex_handle(5)));

	simplices.push_back(Simplex_handle(Vertex_handle(3),Vertex_handle(5) ));

	add_triangle_edges(0,1,5,simplices);
	add_triangle_edges(1,2,3,simplices);
	add_triangle_edges(2,3,4,simplices);
	add_triangle_edges(1,3,4,simplices);
	add_triangle_edges(1,2,4,simplices);


	add_triangle(0,1,5,simplices);
	add_triangle(1,2,3,simplices);
	add_triangle(2,3,4,simplices);
	add_triangle(1,3,4,simplices);
	add_triangle(1,2,4,simplices);


	Complex complex(simplices);

	PRINT(complex.to_string());

	return ( complex.num_vertices()==6&&complex.num_edges()==10&&  complex.num_blockers()==2);
}



int main (int argc, char *argv[])
{
	Tests tests_complex;
	tests_complex.add("test simplex",test_simplex);

	tests_complex.add("test_link0",test_link0);
	tests_complex.add("test_link1",test_link1);
	tests_complex.add("test_link2",test_link2);
	tests_complex.add("test_link3",test_link3);
	tests_complex.add("test_link4",test_link4);
	tests_complex.add("test_link5",test_link5);

	tests_complex.add("test iterator vertices",test_iterator_vertices);
	tests_complex.add("test iterator edges",test_iterator_edge);
	tests_complex.add("test iterator edges2",test_iterator_edge2);
	tests_complex.add("test iterator edges3",test_iterator_edge3);

	tests_complex.add("test iterator blockers",test_iterator_blockers);
	tests_complex.add("test_iterator_triangles",test_iterator_triangles);

	tests_complex.add("test_constructor_list_simplices",test_constructor);

	if(tests_complex.run()){
		return EXIT_SUCCESS;
	}
	else{
		test_iterator_triangles();
		return EXIT_FAILURE;
	}
}

