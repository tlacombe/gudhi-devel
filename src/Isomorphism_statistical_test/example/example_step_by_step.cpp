#include <gudhi/Distance_matrix.h>

#include <iostream>
#include <math.h>
#include <functional>
#include <vector>
#include <array>


using namespace Gudhi::isomorphism_test ;

typedef std::array<double,2> Point;
typedef std::function<double (Point,Point)> Dist;
typedef std::vector<Point> P_range;

int main(){

    int size_P = 100;
    int size_Q = 100;
    double alpha = 0.05;
    unsigned int Nboot = 1000;
    unsigned int nb_redo = 100;


/*  The points are generated on the neighborhood of a spiral according from the following formula :
 (R sin(v R) + sig N,R cos(v R) + sig N')
 with R uniform on [0,1], N and N' normal random variables.*/

  
    double v_P = 10;
    double sig_P = 0.03;
    double v_Q = 10;
    double sig_Q = 0.03;

    std::default_random_engine generator;
    std::normal_distribution<double> distribution_normal(0.0,1.0);
    std::uniform_real_distribution<double> distribution_uniform(0.0,1.0);

    Dist dist = [](Point p1,Point p2) {
        double temp_0 = p1[0]-p2[0];
        double temp_1 = p1[1]-p2[1];
        return sqrt(temp_0*temp_0+temp_1*temp_1);
    };// distance function.

    P_range rP(size_P);
    P_range rQ(size_Q);

    for (int i=0; i<size_P; ++i) {
        double normal_x = distribution_normal(generator);
        double normal_y = distribution_normal(generator);
        double unif = distribution_uniform(generator);
        rP[i][0] = unif * std::sin(v_P*unif) + sig_P * normal_x;
        rP[i][1] = unif * std::cos(v_P*unif) + sig_P * normal_y;
    }
    for (int i=0; i<size_Q; ++i) {
        double normal_x = distribution_normal(generator);
        double normal_y = distribution_normal(generator);
        double unif = distribution_uniform(generator);
        rQ[i][0] = unif * std::sin(v_Q*unif) + sig_Q * normal_x;
        rQ[i][1] = unif * std::cos(v_Q*unif) + sig_Q * normal_y;
    }

    Gudhi::isomorphism_test::Distance_matrix DP;
    DP.build_from_point_set<P_range,Dist>(rP,dist);

    Gudhi::isomorphism_test::Distance_matrix DQ;
    DQ.build_from_point_set<P_range,Dist>(rQ,dist);

/* m Choice */

    std::vector<double> vector_m = {0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8};
    std::vector<bool> signatures_to_plot = {false,true,true,true,true,true,true,true,true}; // plot signatures for all parameter m, except the first one, 0.05.

    auto m_cost = m_choice(DP,DQ,vector_m,true,signatures_to_plot); //plot the curve m->cost(m)

    // print all of the costs
    for(unsigned int it_m = 0; it_m<vector_m.size(); ++it_m){
        std::cout<<"{m,cost} = {"<<vector_m[it_m]<<","<<m_cost[it_m]<<"}."<<std::endl;
    }

    // The user can choose a m !
    std::string m_val;
    double m;
    bool contin = true;
    while(contin){
        std::cout<<"Which m do you want to choose ?"<<std::endl;
        std::cin>>m_val;
        m = stod(m_val);
        if(m<0||m>1){std::cout<<"m should be in [0,1] !"<<std::endl; contin = true;}
        else{contin = false;}
    }

/* Plot the signatures for the m which was chosen */

    plot_signatures(DP, DQ, m);

/* Approximation of the first type error for the sample P*/

    std::vector<std::pair<double,unsigned int>> pair_m_n={{m,5},{m,10},{m,20},{m,30}};

    auto mean_reject = DP.check_test_error(pair_m_n,alpha,Nboot,nb_redo);
    for(unsigned int i=0;i<mean_reject.size();++i){
        std::cout<<"{n,err type I} = {"<<pair_m_n[i].second<<","<<mean_reject[i]<<"} "<<std::endl;
    }

    // The user can choose a n !
    std::string n_val;
    int n;
    contin = true;
    while(contin){
        std::cout<<"Which n do you want to choose ?"<<std::endl;
        std::cin>>n_val;
        n = stoi(n_val);
        if(n<0||n>size_P||n>size_Q){std::cout<<"n should be in [0,min{size_P,size_Q}] !"<<std::endl; contin = true;}
        else{contin = false;}
    }

/* Make the test for the chosen parameters m and n */

    auto tst = test(DP,DQ,m,n,alpha,Nboot);
    std::cout<<"The p-value is equal to : "<<tst.first<<std::endl;
    std::cout<<"The hypothesis retained is H"<<tst.second<<std::endl;

return 1;

}
