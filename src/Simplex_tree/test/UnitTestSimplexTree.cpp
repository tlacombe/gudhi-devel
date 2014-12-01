#define BOOST_TEST_MODULE const_string test
#include <boost/test/included/unit_test.hpp>
#include <boost/system/error_code.hpp>
#include <boost/chrono/thread_clock.hpp>
#include <iostream>
#include <string>

#include <utility> // std::pair, std::make_pair

#include <cmath> // float comparison
#include <limits>

#include "gudhi/graph_simplicial_complex.h"
#include "gudhi/reader_utils.h"
#include "gudhi/Simplex_tree.h"



using namespace Gudhi;

typedef Simplex_tree<> typeST;
typedef std::pair< typeST::Simplex_handle, bool > typePairSimplexBool;
typedef std::vector< Vertex_handle > typeVectorVertex;
typedef std::pair<typeVectorVertex, Filtration_value> typeSimplex;

const Vertex_handle DEFAULT_VERTEX_HANDLE = (const Vertex_handle)-1;

void test_empty_simplex_tree(typeST& tst)
{
	BOOST_CHECK( tst.null_vertex() == DEFAULT_VERTEX_HANDLE );
	BOOST_CHECK( tst.filtration() == Filtration_value(0) );
	BOOST_CHECK( tst.num_vertices() == (size_t)0 );
	BOOST_CHECK( tst.num_simplices() == (size_t)0 );
	typeST::Siblings* STRoot = tst.root();
	BOOST_CHECK( STRoot != NULL );
	BOOST_CHECK( STRoot->oncles() == NULL );
	BOOST_CHECK( STRoot->parent() == DEFAULT_VERTEX_HANDLE );
	BOOST_CHECK( tst.dimension() == -1 );
}

void test_iterators_on_empty_simplex_tree(typeST& tst)
{
	std::cout << "Iterator on vertices: " << std::endl;
	for( auto vertex : tst.complex_vertex_range() )
	{
		std::cout << "vertice:" << vertex << std::endl;
		BOOST_CHECK( false); // shall be empty
	}
	std::cout << "Iterator on simplices: " << std::endl;
	for( auto simplex : tst.complex_simplex_range() )
	{
		BOOST_CHECK( false); // shall be empty
	}

	std::cout << "Iterator on Simplices in the filtration, with [filtration value]:" << std::endl;
	for( auto f_simplex : tst.filtration_simplex_range() )
	{
		BOOST_CHECK( false); // shall be empty
	}
}

BOOST_AUTO_TEST_CASE( simplex_tree_when_empty )
{
	const Filtration_value DEFAULT_FILTRATION_VALUE = 0;

	// TEST OF DEFAULT CONSTRUCTOR
	std::cout << "TEST OF DEFAULT CONSTRUCTOR" << std::endl;
	typeST st;

	test_empty_simplex_tree(st);

	test_iterators_on_empty_simplex_tree(st);
	// TEST OF EMPTY INSERTION
	std::cout << "TEST OF EMPTY INSERTION" << std::endl;
	typeVectorVertex simplexVectorEmpty;
	BOOST_CHECK(simplexVectorEmpty.empty() == true);
	typePairSimplexBool returnEmptyValue = st.insert ( simplexVectorEmpty, DEFAULT_FILTRATION_VALUE );
	BOOST_CHECK(returnEmptyValue.first == typeST::Simplex_handle(NULL));
	BOOST_CHECK(returnEmptyValue.second == true);

	test_empty_simplex_tree(st);

	test_iterators_on_empty_simplex_tree(st);
}

bool AreAlmostTheSame(float a, float b) {
	return std::fabs(a - b) < std::numeric_limits<float>::epsilon();
}

BOOST_AUTO_TEST_CASE( simplex_tree_from_file )
{
	// TEST OF INSERTION
	std::cout << "TEST OF SIMPLEX TREE FROM A FILE" << std::endl;
	typeST st;

	std::ifstream simplex_tree_stream("simplex_tree.txt");
	simplex_tree_stream >> st;

	// Display the Simplex_tree
	std::cout << "The complex contains " << st.num_simplices() << " simplices" << std::endl;
	std::cout << "   - dimension " << st.dimension() << "   - filtration " << st.filtration() << std::endl;
	std::cout << std::endl << std::endl << "Iterator on Simplices in the filtration, with [filtration value]:" << std::endl;
	for( auto f_simplex : st.filtration_simplex_range() )
	{
		std::cout << "   " << "[" << st.filtration(f_simplex) << "] ";
		for( auto vertex : st.simplex_vertex_range(f_simplex) )
		{
			std::cout << vertex << " ";
		}
		std::cout << std::endl;
	}

	//    1
	//    o
	//   / \
	//  o---o---o
	//  2   0   3

	//   [0.1] 0
	//   [0.1] 1
	//   [0.1] 2
	//   [0.1] 3
	//   [0.2] 1 0
	//   [0.2] 2 0
	//   [0.2] 2 1
	//   [0.2] 3 0
	//   [0.3] 2 1 0

	// Check
	BOOST_CHECK(st.num_simplices() == 9);
	BOOST_CHECK(st.dimension() == 2);
	BOOST_CHECK(st.filtration() == 0.3);

	int previous_size=0;
	for( auto f_simplex : st.filtration_simplex_range() )
	{
		// Size of simplex
		int size=0;
		for( auto vertex : st.simplex_vertex_range(f_simplex) )
		{
			size++;
		}
		std::cout << "size=" << size << " - filtration=" << st.filtration(f_simplex) << " - diff=" << (st.filtration(f_simplex)-(0.1* size)) << std::endl;
		BOOST_CHECK(AreAlmostTheSame(st.filtration(f_simplex),(0.1* size))); // Specific test: filtration = 0.1 * simplex_size
		BOOST_CHECK(previous_size <= size); // Check list is sorted (because of sorted filtrations in simplex_tree.txt)
		previous_size = size;
	}
}

void test_simplex_tree_contains(typeST& simplexTree, typeSimplex simplex, int pos)
{
	auto f_simplex = simplexTree.filtration_simplex_range().begin();
	f_simplex += pos;
	std::cout << "test_simplex_tree_contains - filtration=" << simplexTree.filtration(*f_simplex) << "||" << simplex.second << std::endl;
	BOOST_CHECK( AreAlmostTheSame(simplexTree.filtration(*f_simplex),simplex.second) );

	typeVectorVertex::iterator simplexIter = simplex.first.end()-1;
	for( auto vertex : simplexTree.simplex_vertex_range(*f_simplex) )
	{
		std::cout << "test_simplex_tree_contains - vertex=" << vertex << "||" << *simplexIter << std::endl;
		BOOST_CHECK( vertex ==  *simplexIter);
		simplexIter--;
	}
}

BOOST_AUTO_TEST_CASE( simplex_tree_insertion )
{
	Vertex_handle FIRST_VERTEX_HANDLE =  (Vertex_handle)0;
	const Filtration_value FIRST_FILTRATION_VALUE = 0.1;
	Vertex_handle SECOND_VERTEX_HANDLE = (Vertex_handle)1;
	const Filtration_value SECOND_FILTRATION_VALUE = 0.1;
	const Filtration_value THIRD_FILTRATION_VALUE = 0.2;

	// TEST OF INSERTION
	std::cout << "TEST OF INSERTION" << std::endl;
	typeST st;

	std::cout << "   - INSERT 0" << std::endl;
	typeVectorVertex firstSimplexVector;
	firstSimplexVector.push_back(FIRST_VERTEX_HANDLE);
	BOOST_CHECK( firstSimplexVector.size() == 1 );
	typeSimplex firstSimplex = std::make_pair(firstSimplexVector, Filtration_value(FIRST_FILTRATION_VALUE));
	typePairSimplexBool returnValue =
			st.insert ( firstSimplex.first, firstSimplex.second );

	BOOST_CHECK(returnValue.second == true);
	typeST::Simplex_handle shReturned = returnValue.first; // Simplex_handle = boost::container::flat_map< Vertex_handle, Node >::iterator
	BOOST_CHECK(shReturned != typeST::Simplex_handle(NULL));

	st.set_dimension(firstSimplexVector.size()-1);
	st.set_num_simplices(1);
	st.set_filtration(firstSimplex.second);

	BOOST_CHECK( st.dimension() == (size_t)(firstSimplexVector.size()-1) );
	BOOST_CHECK( st.num_simplices() == (size_t)1 );
	BOOST_CHECK( st.filtration() == firstSimplex.second );
	BOOST_CHECK( st.num_vertices() == (size_t)1 );

	// Display the Simplex_tree
	std::cout << "The complex contains " << st.num_simplices() << " simplices" << std::endl;
	std::cout << "   - dimension " << st.dimension() << "   - filtration " << st.filtration() << std::endl;
	std::cout << std::endl << std::endl << "Iterator on Simplices in the filtration, with [filtration value]:" << std::endl;

	//test_simplex_tree_contains(st, firstSimplex, 0);

	std::cout << "   - INSERT 1" << std::endl;
	typeVectorVertex secondSimplexVector;
	secondSimplexVector.push_back(SECOND_VERTEX_HANDLE);
	BOOST_CHECK( secondSimplexVector.size() == 1 );
	typeSimplex secondSimplex = std::make_pair(secondSimplexVector, Filtration_value(SECOND_FILTRATION_VALUE));
	returnValue =
			st.insert ( secondSimplex.first, secondSimplex.second );

	BOOST_CHECK(returnValue.second == true);
	shReturned = returnValue.first; // Simplex_handle = boost::container::flat_map< Vertex_handle, Node >::iterator
	BOOST_CHECK(shReturned != typeST::Simplex_handle(NULL));

	st.set_dimension(secondSimplexVector.size()-1);
	st.set_num_simplices(2);
	st.set_filtration(secondSimplex.second);

	BOOST_CHECK( st.dimension() == (size_t)(secondSimplexVector.size()-1) );
	BOOST_CHECK( st.num_simplices() == (size_t)2 );
	BOOST_CHECK( st.filtration() == secondSimplex.second );
	BOOST_CHECK( st.num_vertices() == (size_t)2 );

	std::cout << "   - INSERT (0,1)" << std::endl;
	typeVectorVertex thirdSimplexVector;
	thirdSimplexVector.push_back(FIRST_VERTEX_HANDLE);
	thirdSimplexVector.push_back(SECOND_VERTEX_HANDLE);
	BOOST_CHECK( thirdSimplexVector.size() == 2 );
	typeSimplex thirdSimplex = std::make_pair(thirdSimplexVector, Filtration_value(THIRD_FILTRATION_VALUE));
	returnValue =
			st.insert ( thirdSimplex.first, thirdSimplex.second );

	BOOST_CHECK(returnValue.second == true);
	shReturned = returnValue.first; // Simplex_handle = boost::container::flat_map< Vertex_handle, Node >::iterator
	BOOST_CHECK(shReturned != typeST::Simplex_handle(NULL));

	st.set_dimension(thirdSimplexVector.size()-1);
	st.set_num_simplices(3);
	st.set_filtration(thirdSimplex.second);


	BOOST_CHECK( st.dimension() == (size_t)(thirdSimplexVector.size()-1) );
	BOOST_CHECK( st.num_simplices() == (size_t)3 );
	BOOST_CHECK( st.filtration() == thirdSimplex.second );
	BOOST_CHECK( st.num_vertices() == (size_t)2 );

	test_simplex_tree_contains(st, firstSimplex, 0);
	test_simplex_tree_contains(st, secondSimplex, 1);
	test_simplex_tree_contains(st, thirdSimplex, 2);

	// Display the Simplex_tree - Can not be done in the middle of 2 inserts
	std::cout << "The complex contains " << st.num_simplices() << " simplices" << std::endl;
	std::cout << "   - dimension " << st.dimension() << "   - filtration " << st.filtration() << std::endl;
	std::cout << std::endl << std::endl << "Iterator on Simplices in the filtration, with [filtration value]:" << std::endl;
	for( auto f_simplex : st.filtration_simplex_range() )
	{
		std::cout << "   " << "[" << st.filtration(f_simplex) << "] ";
		for( auto vertex : st.simplex_vertex_range(f_simplex) )
		{
			std::cout << (int)vertex << " ";
		}
		std::cout << std::endl;
	}

}
