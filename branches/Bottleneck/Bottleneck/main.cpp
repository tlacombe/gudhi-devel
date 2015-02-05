#include "Neighbors_finder.hpp"
#include "Multi_nf.hpp"
#include "Deque.hpp"
#include "Graph_matching.hpp"
#include <iostream>
#include <vector>
#include <list>


void print_point(Diagram_point q){
    std::cout << "(" << q.first << ";" << q.second << ")";
}

void print_point(Persistence_diagram& p, int point_index){
    std::cout << " {" << point_index << " : ";
    print_point(p.get_point(point_index));
    std::cout << "} " << std::endl;
}

void print_point(Persistence_diagram& p, std::list<int>* l){
    for(auto it= l->cbegin(); it != l->cend(); ++it)
        print_point(p,*it);
}

bool naive_nf_test()
{
    std::vector<Diagram_point> pv;
    pv.emplace_back(1.,2.);
    pv.emplace_back(5.,7.);
    pv.emplace_back(1.,3.5);
    pv.emplace_back(5.,7.);
    pv.emplace_back(5.,5.);
    pv.emplace_back(5.,7.);

    std::vector<Diagram_point> qv;
    qv.emplace_back(1.,3.);
    qv.emplace_back(5.,7.);
    qv.emplace_back(1.,3.);
    qv.emplace_back(4.,6.);
    qv.emplace_back(4.5,7.5);

    int n = 5;

    qv.emplace_back(1.,3.);
    qv.emplace_back(5.,6.);
    qv.emplace_back(10.,8.);

    double r = 1.;
    Persistence_diagram x(pv);
    Naive_nf snf(x,r);
    for(int i = 0; i< (int) x.size(); ++i)
        snf.add(i);
    for(int i = 0; i<n; ++i){
        Diagram_point q = qv.at(i);
        int ia = snf.pull_near(q);
        if(ia == null_point_index())
            return false;
        Diagram_point a = x.get_point(ia);
        if(distance(qv.at(i),a)>r)
            return false;
    }
    for(int i = 5; i< (int) qv.size(); ++i){
        Diagram_point q = qv.at(i);
        int ia = snf.pull_near(q);
        if(ia != null_point_index())
            return false;
    }
    std::list<int>* l = snf.pull_all_near(qv.at(0));
    if(l->size()!=0)
        return false;
    return true;
}

bool multi_nf_test(){
    std::vector< Diagram_point > v1, v2;
    int n = 40;
    for(int i = 0; i<100; ++i){
        int a = rand()%n;
        v1.emplace_back(a, a + rand()%(n-a));
        int b = rand()%n;
        v2.emplace_back(b, b + rand()%(n-b));
    }
    Persistence_diagram x(v1);
    Naive_nf snf(x, 2.6);
    Multi_nf mnf (x, 2.6);
    mnf.add_layer();
    mnf.add_layer();
    mnf.add_layer();
    for(int i = 0; i< (int) x.size(); ++i){
        snf.add(i);
        mnf.add(i,1);
    }
    for(auto it = v2.cbegin(); it != v2.cend(); ++it){
        std::list<int>* l1 = snf.pull_all_near(*it);
        std::list<int>* l2 = mnf.pull_all_near(*it,1);
        std::list<int>* l3 = mnf.pull_all_near(*it,0);
        if(l3->size()!=0)
            return false;
        l1->sort();
        l2->sort();
        auto it1 = l1->cbegin();
        auto it2 = l2->cbegin();
        for(int i = 0; i< (int) l1->size(); ++i){
            if(*it1 != *it2)
                return false;
            ++it1;
            ++it2;
        }
    }
    return true;
}

void all_tests()
{
    if(naive_nf_test())
        std::cout << "Naive_nf seems OK" << std::endl;
    else
        std::cout << "Naive_nf is not correct !" << std::endl;
    if(multi_nf_test())
        std::cout << "Multi_nf seems OK" << std::endl;
    else
        std::cout << "Multi_nf is not correct !" << std::endl;
}

int main()
{
    all_tests();
}
