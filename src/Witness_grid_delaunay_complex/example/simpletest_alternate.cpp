#define CGAL_EIGEN3_ENABLED
#include <gudhi/Simplex_tree.h>
#include <gudhi/distance_functions.h>
#include <gudhi/Point.h>
#include <gudhi/reader_utils.h>
#include <boost/math/constants/constants.hpp>
#include <algorithm>
#include <utility>
#include <vector>
#include <list>
#include <set>
#include <queue>
#include <limits>
#include <ctime>
#include <cmath>
#include <iostream>
#include <CGAL/point_generators_d.h>
#include <CGAL/Cartesian_d.h>
#include <CGAL/Epick_d.h>
#include <CGAL/constructions_d.h>
#include <CGAL/Kernel_d/Hyperplane_d.h>
#include <CGAL/Origin.h>
#include "output.h"
#include "sampling_radius.h"
#include <Delaunay_Triangulation_BDG15_alternate.h>
int main()
{
    typedef CGAL::Cartesian_d<double> Kd;
    typedef CGAL::Epick_d<CGAL::Dynamic_dimension_tag> Ed;
    typedef Ed::Point_d Point_ed;
    typedef Kd::Point_d Point;
    typedef std::vector<Point> Point_Vector;
    typedef Gudhi::Simplex_tree<> Simplicial_complex;
    typedef Gudhi::Simplex_tree<>::Simplex_handle Simplex_handle;
    typedef Gudhi::Simplex_tree<>::Complex_simplex_iterator Complex_simplex_iterator;
    Point_Vector L;
    int dim;
    std::cin>>dim;
    CGAL::Random_points_in_ball_d<Point> rp(dim, 1);//random points in unit ball well separated
int siz;std::cin>>siz;
    double sp;std::cin>>sp;
    //std::cin>>siz>>sp;
    L.push_back(*rp++);
    while(L.size()!=siz)
    {
        Point temp=*rp;
        double spar=sp;
        int j;
        for(j=0; j<L.size(); j++)
        {
            if(squared_distance(L[j],temp)<=spar*spar)
            {
                rp++;
                break;
            }
        }
        if(j==L.size())
        {
            L.push_back(*rp++);
        }
    }
    //std::cout<<"Points created"<<std::endl;
    write_points("point.mesh",L);
    std::vector<Ed::Point_d> points;
    for(int i=0; i<L.size(); i++)
    {
        //std::cout<<L[i]<<std::endl;
        for(int j=0; j<L[i].dimension(); j++)
        {
            auto st=L[i].cartesian_begin();
            auto en=L[i].cartesian_end();
            points.push_back(Point_ed(st,en));
        }
    }
    std::pair<double,Point_ed> lambdapair=sampling_radius<Ed>(points);
    //std::cout<<"lambda  "<<lambdapair.first<<" "<<lambdapair.second<<std::endl;
    double lambda=lambdapair.first;
    Gudhi::Simplex_tree<> simplex_tree;
    double rho=sp/(4*lambda);//std::cin>>rho;
    double eps;std::cin>>eps;
clock_t t=clock();
    Gudhi::Delaunay_Triangulation_BDG15::Delaunay_Triangulation_BDG15<Simplicial_complex> T(simplex_tree,L,dim,lambda,eps,rho);
t=clock()-t;
std::cout<<"Parameters(dim,size,eps)"<<dim<<" "<<siz<<" "<<eps<<"  "<<"Time taken : "<< double(t)/CLOCKS_PER_SEC<<std::endl;   
    std::vector<int> land;
    for(int y=0; y<L.size(); y++)
    {
        land.push_back(y);
    }
    write_witness_mesh(L,land, simplex_tree,simplex_tree.complex_simplex_range(),dim==2,true,"complex.mesh" );
    return 0;
}

