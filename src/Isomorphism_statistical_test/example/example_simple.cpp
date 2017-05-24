#include <gudhi/Distance_matrix.h>

#include <iostream>
#include <math.h>
#include <functional>
#include <vector>

using namespace Gudhi::isomorphism_test;


int main(){
    bool check = true; // to check is the matrices rmatP and rmatQ are symetrical with zeros on the diagonal.
    std::vector<double> rmatP = {0,1,1.41421,2.23607,1,0,1,2,1.41421,1,0,1,2.23607,2,1,0};
    std::vector<double> rmatQ = {0,1.1,1.5,2.23,1.1,0,1,2,1.5,1,0,1,2.23,2,1,0};

    Distance_matrix DP;
    DP.build_from_distance_matrix<std::vector<double>>(rmatP,check);

    Distance_matrix DQ;
    DQ.build_from_distance_matrix<std::vector<double>>(rmatQ,check);

    int n = 2; // 0<n<=N where N = 4 is the number of points in a sample
    double m = 0.4; // mass parameter in [0,1]
    double alpha = 0.05; // level of the test
    int Nboot = 1000; // number of steps for the bootstrap, in general 1000 is enough. 
    auto tst = test(DP,DQ,m,n,alpha,Nboot);

    std::cout<<"The p-value is equal to : "<<tst.first<<std::endl;
    std::cout<<"The hypothesis retained is H"<<tst.second<<std::endl;

}
