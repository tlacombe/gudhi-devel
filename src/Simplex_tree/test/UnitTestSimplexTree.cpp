#define BOOST_TEST_MODULE const_string test
#include <boost/test/included/unit_test.hpp>
#include <boost/system/error_code.hpp>
#include <boost/chrono/thread_clock.hpp>
#include <iostream>
#include <string>

#include <ctime>
#include "gudhi/graph_simplicial_complex.h"
//#include "gudhi/distance_functions.h"
#include "gudhi/Simplex_tree.h"

//using namespace boost::posix_time;
using namespace Gudhi;

//typedef int        Vertex_handle;
//typedef double     Filtration_value;

typedef Simplex_tree<> typeST;
typedef std::pair< typeST::Simplex_handle, bool > typePairSimplexBool;
typedef std::vector< Vertex_handle > typeVectorVertex;
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

BOOST_AUTO_TEST_CASE( simplex_tree_test )
{
   // TEST OF DEFAULT CONSTRUCTOR
   std::cout << "TEST OF DEFAULT CONSTRUCTOR" << std::endl;
   typeST st;

   test_empty_simplex_tree(st);

   // TEST OF EMPTY INSERTION
   std::cout << "TEST OF EMPTY INSERTION" << std::endl;
   typeVectorVertex simplexVectorEmpty;
   BOOST_CHECK(simplexVectorEmpty.empty() == true);
   const Filtration_value DEFAULT_FILTRATION_VALUE = 0;
   typePairSimplexBool returnEmptyValue = st.insert ( simplexVectorEmpty, DEFAULT_FILTRATION_VALUE );
   BOOST_CHECK(returnEmptyValue.first == typeST::Simplex_handle(NULL));
   BOOST_CHECK(returnEmptyValue.second == true);

   test_empty_simplex_tree(st);

   // TEST OF INSERTION
   std::cout << "TEST OF INSERTION" << std::endl;
   typeVectorVertex simplexVector;
   const Vertex_handle FIRST_VERTEX_HANDLE = (const Vertex_handle)1;
   simplexVector.push_back(FIRST_VERTEX_HANDLE);
   BOOST_CHECK(simplexVector.empty() == false);
   BOOST_CHECK(simplexVector.size() == 1);
   const Filtration_value FIRST_FILTRATION_VALUE = 1.5;
   typePairSimplexBool returnValue = st.insert ( simplexVector, Filtration_value(FIRST_FILTRATION_VALUE) );

   typeST::Simplex_handle shReturned = returnValue.first; // Simplex_handle = boost::container::flat_map< Vertex_handle, Node >::iterator
   BOOST_CHECK(shReturned != typeST::Simplex_handle(NULL));


   BOOST_CHECK( st.null_vertex() == DEFAULT_VERTEX_HANDLE );
   BOOST_CHECK( st.filtration() == DEFAULT_FILTRATION_VALUE );
   BOOST_CHECK( st.num_vertices() == (size_t)1 );
   BOOST_CHECK( st.num_simplices() == (size_t)0 );
   typeST::Siblings* STRoot = st.root();
   BOOST_CHECK( STRoot != NULL );
   BOOST_CHECK( STRoot->oncles() == NULL );
   BOOST_CHECK( STRoot->parent() == DEFAULT_VERTEX_HANDLE );
   BOOST_CHECK( st.dimension() == -1 );

   BOOST_CHECK(returnValue.second == true);
}
