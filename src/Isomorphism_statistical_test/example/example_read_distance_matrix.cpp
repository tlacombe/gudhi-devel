#include <gudhi/Distance_matrix.h>

#include <iostream>
#include <math.h>
#include <functional>
#include <vector>

using namespace Gudhi::isomorphism_test;

int main(){

    std::cout<<std::endl;
    std::cout<<"Different ways of reading a distance matrix."<<std::endl;
    std::cout<<std::endl;

/*******************************************************************************************************************************************************************
                                   build_from_point_set<typename P_range,typename Dist>(P_range const & rP, Dist const & dist)
*******************************************************************************************************************************************************************/

    typedef std::vector<double> Point;
    typedef std::function<double (Point,Point)> Dist;
    typedef std::vector<Point> P_range;

    const int dim = 2;
    P_range rP = {{1,0},{0,0},{0,1},{0,2}};// Here a vector of vectors of size dim.
    Dist dist = [dim](Point p1,Point p2) {
        double sum_sq = 0;
        double temp = 0;
        for (int i=0; i<dim; ++i){
            temp = p1[i]-p2[i];
            sum_sq = sum_sq + temp*temp;
        }
        return sqrt(sum_sq);
    };

    Distance_matrix D1;
    D1.build_from_point_set<P_range,Dist>(rP,dist,true,"matrix1.txt"); // The matrix will be stored in the file "matrix1.txt" if this file does not exist already.
    D1.printmat(); // To print the distance matrix.



/*******************************************************************************************************************************************************************
                              build_from_file_point_set<typename Dist>(std::string const & file_name, Dist const & dist, double dimension)
********************************************************************************************************************************************************************/

    typedef std::vector<double> Point2;
    typedef std::function<double (Point,Point)> Dist2;

    const int dimension = 2;
    Dist2 dist2 = [dimension](Point2 p1,Point2 p2) {
        double sum_sq = 0;
        double temp = 0;
        for (int i=0; i<dimension; ++i){
            temp = p1[i]-p2[i];
            sum_sq = sum_sq + temp*temp;
        }
        return sqrt(sum_sq);
    };

    Distance_matrix D2;
    D2.build_from_file_point_set("point_set.txt",dist2,dimension,true,"matrix2.txt"); // The matrix will be stored in the file "matrix2.txt" if this file does not exist already.
    D2.printmat(); 


/*******************************************************************************************************************************************************************
                             build_from_distance_matrix<typename Matrix_range>(Matrix_range & rmat, bool check_symmetry_and_null_diagonal);
*******************************************************************************************************************************************************************/

    bool check = true;
    std::vector<double> rmat = {0,1,1.41421,2.23607,1,0,1,2,1.41421,1,0,1,2.23607,2,1,0};

    Distance_matrix D3;
    D3.build_from_distance_matrix<std::vector<double>>(rmat,check);
    D3.printmat(); 


/*******************************************************************************************************************************************************************
                                    build_from_distance_lower_matrix<typename Lower_matrix_range(Lower_matrix_range & rmat);
*******************************************************************************************************************************************************************/

    std::vector<double> rmat2 = {1,1.41421,1,2.23607,2,1};

    Distance_matrix D4;
    D4.build_from_distance_lower_matrix<std::vector<double>>(rmat2);
    D4.printmat();


/*******************************************************************************************************************************************************************
                                        build_from_file_matrix(std::string const & file_name, bool check_symmetry_and_null_diagonal);
*******************************************************************************************************************************************************************/

    bool check2 = true;

    Distance_matrix D5;
    D5.build_from_file_matrix("distance_matrix.txt", check2);
    D5.printmat();


/*******************************************************************************************************************************************************************
                                                void build_from_file_lower_matrix(std::string const & file_name);
*******************************************************************************************************************************************************************/

    Distance_matrix D6;
    D6.build_from_file_lower_matrix("distance_lower_matrix.txt");
    D6.printmat();

}
