#include <iostream>
#include <gudhi/LSAL.h>

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE "simplex array list"
#include <boost/test/unit_test.hpp>

using namespace Gudhi;


std::vector<int> sigma1 = {1, 2, 3, 4};
std::vector<int> sigma2 = {5, 2, 3, 6};
std::vector<int> sigma3 = {5};
std::vector<int> sigma4 = {5, 2, 3};
std::vector<int> sigma5 = {5, 2, 7};
std::vector<int> sigma6 = {4, 5, 3};
std::vector<int> sigma7 = {4, 5, 9};
std::vector<int> sigma8 = {1, 2, 3, 6};



BOOST_AUTO_TEST_CASE(sal) {
    SAL K;
    K.insert_simplex(sigma1);
    K.insert_simplex(sigma2);
    K.insert_simplex(sigma3);
    K.insert_simplex(sigma6);
    K.insert_simplex(sigma7);
    BOOST_CHECK(K.membership(sigma4));
    BOOST_CHECK(!K.criticality(sigma5));
    BOOST_CHECK(!K.membership(sigma5));
    K.contraction(4,5);
    BOOST_CHECK(!K.membership(sigma6));
}

BOOST_AUTO_TEST_CASE(lsal) {
    LSAL K;
    K.insert_simplex(sigma1);
    K.insert_simplex(sigma2);
    K.insert_simplex(sigma3);
    K.insert_simplex(sigma6);
    K.insert_simplex(sigma7);
    BOOST_CHECK(K.membership(sigma4));
    BOOST_CHECK(!K.membership(sigma5));
    K.contraction(4,5);
    BOOST_CHECK(!K.membership(sigma6));
}

BOOST_AUTO_TEST_CASE(sal_candidates) {
    SAL K;
    K.insert_simplex(sigma1);
    K.insert_simplex(sigma2);
    K.remove_simplex(sigma1);
    K.remove_simplex(sigma2);
    auto c = K.candidates();
    BOOST_CHECK(c.count(get_key(sigma1)));
    BOOST_CHECK(c.count(get_key(sigma2)));
    BOOST_CHECK(c.size()==2);
}
