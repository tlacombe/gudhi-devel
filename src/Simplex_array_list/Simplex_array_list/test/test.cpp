#include <iostream>
#include <gudhi/LSAL.h>
#include <gudhi/STW.h>

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE "simplex array list"
#include <boost/test/unit_test.hpp>

using namespace Gudhi;


std::unordered_set<Vertex> sigma1 = {1, 2, 3, 4};
std::unordered_set<Vertex> sigma2 = {5, 2, 3, 6};
std::unordered_set<Vertex> sigma3 = {5};
std::unordered_set<Vertex> sigma4 = {5, 2, 3};
std::unordered_set<Vertex> sigma5 = {5, 2, 7};
std::unordered_set<Vertex> sigma6 = {4, 5, 3};
std::unordered_set<Vertex> sigma7 = {4, 5, 9};


BOOST_AUTO_TEST_CASE(sal) {
    SAL K;
    K.add(sigma1);
    K.add(sigma2);
    K.add(sigma3);
    K.add(sigma6);
    K.add(sigma7);
    BOOST_CHECK(K.membership(sigma4));
    BOOST_CHECK(!K.membership(sigma5));
    K.contraction(4,5);
    /* for(auto s : K.max_cofaces(sigma3)){
        for(int p : s)
            std::cout << p << " ";
        std::cout<< std::endl;
    } */
}

BOOST_AUTO_TEST_CASE(lsal) {
    LSAL LK;
    LK.add(sigma1);
    LK.add(sigma2);
    LK.add(sigma3);
    LK.add(sigma6);
    LK.add(sigma7);
    BOOST_CHECK(LK.membership(sigma4));
    BOOST_CHECK(!LK.membership(sigma5));
    LK.contraction(4,5);
}

BOOST_AUTO_TEST_CASE(stw) {
    STW st;
    st.add(sigma1);
    st.add(sigma2);
    st.add(sigma3);
    st.add(sigma6);
    st.add(sigma7);
    BOOST_CHECK(st.membership(sigma4));
    BOOST_CHECK(!st.membership(sigma5));
    st.contraction(4,5);
}

