#include <iostream>
#include <gudhi/LSAL.h>
#include <gudhi/STW.h>
#include <gudhi/SALF.h>


#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE "simplex array list"
#include <boost/test/unit_test.hpp>

using namespace Gudhi;


Simplex sigma1 = {1, 2, 3, 4};
Simplex sigma2 = {5, 2, 3, 6};
Simplex sigma3 = {5};
Simplex sigma4 = {5, 2, 3};
Simplex sigma5 = {5, 2, 7};
Simplex sigma6 = {4, 5, 3};
Simplex sigma7 = {4, 5, 9};


BOOST_AUTO_TEST_CASE(sal) {
    SAL K;
    K.add(sigma1);
    K.add(sigma2);
    K.add(sigma3);
    K.add(sigma6);
    K.add(sigma7);
    BOOST_CHECK(K.membership(sigma4));
    BOOST_CHECK(!K.membership(sigma5));
    Vertex v = K.contraction(4,5);
    Simplex s = {2, 3, 6}; s.insert(v);
    K.remove(s);
    BOOST_CHECK(K.all_facets_inside(s));
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
    Vertex v = LK.contraction(4,5);
    Simplex s = {2, 3, 6}; s.insert(v);
    LK.remove(s);
    BOOST_CHECK(LK.all_facets_inside(s));
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

BOOST_AUTO_TEST_CASE(salf) {
    SALF FK;
    FK.add(sigma1, 2.);
    FK.add(sigma2, 2.);
    FK.add(sigma3, 2.);
    FK.add(sigma6, 2.);
    FK.add(sigma7, 2.);
    BOOST_CHECK(FK.membership(sigma4));
    BOOST_CHECK(!FK.membership(sigma5));
}
