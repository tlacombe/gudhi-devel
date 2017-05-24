/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Claire Brécheteau
 *
 *    Copyright (C) 2017  INRIA Sophia-Saclay (France)
 *
 *    This program is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *
 *    This program is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef ISOMORPHISM_TEST_H_
#define ISOMORPHISM_TEST_H_

//-Wall-fopenmp

// Choisir un meilleur titre pour la branche. A statistical test of isomorphism with the DTM.

// Marc : ca compile avec gcc-6 ?


#include <iostream>
#include <string>
#include <fstream>
#include <ostream>
#include <vector>
#include <iterator>
#include <cmath>
#include <algorithm>
#include <math.h> 
#include <random>

#include <utility>
#include <GNUPLOT/gnuplot-iostream.h>



namespace Gudhi {
namespace isomorphism_test {

/**
* Check that the size of mat is a square : size*size.
* Check that the matrix \f$(mat[i*size+j])_{i,j}\f$ is symmetric with zeros on the diagonal.
* If some of these conditions is not satisfied, the program is aborted. 
**/
void check_symmetry_and_null_diagonal_of_matrix(std::vector<double> const & mat){
    double sqrt_size = std::sqrt(mat.size());
    int sqrt_size_int = int(sqrt_size);
    if (sqrt_size != sqrt_size_int){
        std::cerr<<"Warning : The matrix is not squared !\n";
        std::abort();
    }
    # pragma omp parallel for
    for (int i=0; i<sqrt_size_int; ++i){
        for (int j=0; j<i; ++j){
            if (mat[sqrt_size*i+j]!=mat[sqrt_size*j+i]){
                std::cerr<<"Warning : The matrix of distances is not symmetrical !\n";
                std::abort();
            }
        }
    }
    for (int i=0; i<sqrt_size_int; ++i){
        if (mat[sqrt_size*i+i]!=0){
            std::cerr<<"Warning : The diagonal of the matrix should be zero !\n";
            std::abort();
        }
    }
}

/**
* Read a file which contains doubles.
* Store these doubles in the vector data.
* It will be used to read files containing either the elements of a distance matrix, or the coordinates of points in \f$\mathbf{R}^d\f$.
**/
void read_file(std::vector<double> & data,std::string const & file_name){
    std::ifstream file(file_name, std::ios::in);
    if(file){
        double data_item;
        while(file>>data_item){
            data.push_back(data_item);
        };
        file.close();
    }
    else
    std::cerr << "\n Unable to open the file " << file_name <<" !\n" << std::endl;
}

/**
* Function to store the elements of the matrix mat in a file.
* Remark that mat[i*size+j]=d(X_i,X_j) but the elements writen to the file file_name_matrix are : d(X_2,X_1), d(X_3,X_2), d(X_3,X_1), d(X_4,X_1), d(X_4,X_2)...
**/
void write_file(std::vector<double> & mat, std::string const & file_name_matrix){
    std::string file_name_matrix_new = file_name_matrix;
    std::ifstream file_matrix_test(file_name_matrix_new, std::ios::out);
    bool continu = file_matrix_test.is_open();
    if(continu){file_matrix_test.close();}
    while(continu){
        std::string stop_loop = "";
        while(stop_loop!="y" && stop_loop!="n"){
            std::cout<<"The file "<<file_name_matrix_new<<" already exists and will be erased if you continue. Do you want to continue ? (y/n)"; 
            std::cin>>stop_loop;
        }
        if(stop_loop == "n"){
            std::cout<<"Please choose another file name for the matrix of distances :";
            std::cin>>file_name_matrix_new;
            std::ifstream file_matrix_test(file_name_matrix_new, std::ios::out);
            continu = file_matrix_test.is_open();
            if(continu){file_matrix_test.close();}
        }
        if(stop_loop == "y"){
            continu = false;
        }
    }
    std::ofstream file_matrix(file_name_matrix_new, std::ios::out);
    if(!file_matrix){std::cerr<<"Can't write on file "<<file_name_matrix_new<<" "; std::abort();}
    int size = std::sqrt(mat.size());
    if(size*size!=mat.size()){std::cerr<<"mat should be a squared matrix !"; std::abort();}
    for(unsigned int i=0; i<size; ++i){
        for (unsigned int j=0; j<i; ++j){
            file_matrix<<mat[i*size+j]<<" ";
        }
    }
    file_matrix.close();  
}

/*********************************************************************************************************************************************************************

                                                                 CLASS Distance_matrix

*********************************************************************************************************************************************************************/

/**
 * @brief Distance_matrix class for the isomorphism statistical test.
 * @ingroup isomorphism_test
 * @details The isomorphism statistical test is built from two samples coming from spaces equipped with a metric.
 * The input for the test is a pair of two distance matrices.
 * Each distance matrix contains the distances between all pairs of points within a sample.
 *
 * This class contains many methods to read or build distance matrices, from a C++ code or from a file.
 * It contains a method which is the implementation of the statistical test of isomorphism.
 * It also contains a method to help choosing the mass parameter m.
 * It finally contains a method based on a heuristic to check whether the test has correct first type error or not. By correct, we mean close to some fixed parameter alpha.
 * 
 * If this last method returns an error of type I too big, or too far from alpha,
 * then it means that the samples are not large enough to answer the question of isomorphism with the use of this statistical test.
 * We cannot expect too much of this test.
 * Indeed theoretical results justifying the performance of the test are asymptotic, the aforementioned method is a tool to check that it has good chances to work or not.
 **/

class Distance_matrix {//
public:
    Distance_matrix(): size_(0), mat_(0), mat_is_sorted_(0), mat_sorted_(0), vect_dtm_() {}

/**
* Prints the matrix of distances mat_, one raw per line.
**/
    void printmat(){
        if (mat_.size()==0){std::cerr<<"The matrix mat_ is still empty ! "; std::abort();}
        for(unsigned int i=0;i<size_;++i){
            for(unsigned int j=0;j<size_;++j){
                std::cout<<mat_[j*size_+i]<<" ";
            }
            std::cout<<std::endl;
        }
        std::cout<<std::endl;
    }

/**
* Matrix_range is a readable forward range of doubles.
* rmat should contain the distances between all pairs of points in a sample : \f$d(X_1,X_1)\f$, \f$d(X_1,X_2)\ldots\f$ \f$d(X_1,X_N)\f$, \f$d(X_2,X_1)\f$,
* \f$d(X_2,X_2)\ldots\f$ \f$d(X_2,X_N)\f$ etc.
* rmat is of size \f$N^2\f$.
*
* If check_symmetry_and_null_diagonal is TRUE, then the method checks that the matrix \f$(rmat[iN+j])_{i,j}\f$ is symmetrical with zeros on its diagonal.
* Uses the function check_symmetry_and_null_diagonal_of_matrix.
**/
    template<typename Matrix_range> void build_from_distance_matrix(Matrix_range const & rmat, bool check_symmetry_and_null_diagonal);

/**
* Lower_matrix_range is a readable forward range of doubles.
* rmat should contain the distances between all pairs of points in a sample : \f$d(X_2,X_1)\f$, \f$d(X_3,X_1)\f$, \f$d(X_3,X_2)\f$, \f$d(X_4,X_1)\f$, \f$d(X_4,X_2)\f$, \f$d(X_4,X_3)\f$ etc.
* rmat is of size \f$\frac{N(N-1)}{2}\f$.
* Note that \f$d(X_1,X_2)\f$ and \f$d(X_2,X_1)\f$ should appear only once, and \f$d(X_1,X_1)\f$ should not appear.
**/
    template<typename Lower_matrix_range> void build_from_distance_lower_matrix(Lower_matrix_range const & rmat);

/**
* Reads a file file_name containing all the distances between pairs of points in a sample : \f$d(X_1,X_1)\f$, \f$d(X_1,X_2)\ldots\f$ \f$d(X_1,X_N)\f$,
* \f$d(X_2,X_1)\f$, \f$d(X_2,X_2)\ldots\f$ \f$d(X_2,X_N)\f$ etc.
* The number of doubles in file_name is N^2.
*
* If check_symmetry_and_null_diagonal is TRUE, then the method checks that the matrix \f$(rmat[iN+j])_{i,j}\f$ is symmetrical with zeros on its diagonal.
* Uses the function check_symmetry_and_null_diagonal_of_matrix.
* Uses the function read_file.
**/
    void build_from_file_matrix(std::string const & file_name, bool check_symmetry_and_null_diagonal);

/**
* Reads a file file_name containing all the distances between pairs of points in a sample : \f$d(X_2,X_1)\f$, \f$d(X_3,X_1)\f$, \f$d(X_3,X_2)\f$, \f$d(X_4,X_1)\f$,
* \f$d(X_4,X_2)\f$, \f$d(X_4,X_3)\f$ etc.
* The number of doubles in file_name is \f$\frac{N(N-1)}{2}\f$.
* Note that \f$d(X_1,X_2)\f$ and \f$d(X_2,X_1)\f$ should appear only once, and \f$d(X_1,X_1)\f$ should not appear.
*
* Uses the function read_file.
**/
    void build_from_file_lower_matrix(std::string const & file_name);

/**
* P_range is a readable forward range of Points (any type).
* dist is a function with arguments two Points and which returns a double, dist is actually a distance on a metric space which elements are of type Points.
* If store_distance_matrix is true, one stores the computed distance matrix in the file file_name.
* The data stored are : \f$d(X_2,X_1)\f$, \f$d(X_3,X_1)\f$, \f$d(X_3,X_2)\f$, \f$d(X_4,X_1)\f$, \f$d(X_4,X_2)\f$, \f$d(X_4,X_3)\f$ etc.
**/
    template<typename P_range,typename Dist> void build_from_point_set(P_range const & rP, Dist const & dist, bool store_distance_matrix = false, std::string const & file_name = "distance_matrix.txt");

/**
* Reads a file file_name containing all the coordinates of \f$N\f$ points in \f$\mathbf{R}^{dimension}\f$ : \f$X_{1,1}\f$, \f$X_{1,2}\ldots\f$ \f$X_{1,dimension}\f$, \f$X_{2,1}\f$, \f$X_{2,2}\ldots\f$ \f$X_{2,dimension}\f$ etc.
* The number of doubles in file_name is \f$N dimension\f$.
* dist is a function with argument two vectors of \f$\mathbf{R}^{dimension}\f$.
*
* If store_distance_matrix is true, one stores the computed distance matrix in the file file_name.
* The data stored are : \f$d(X_2,X_1)\f$, \f$d(X_3,X_1)\f$, \f$d(X_3,X_2)\f$, \f$d(X_4,X_1)\f$, \f$d(X_4,X_2)\f$, \f$d(X_4,X_3)\f$ etc.
*
* Uses the function read_file.
**/
    void build_from_file_point_set(std::string const & file_name,  std::function<double (std::vector<double>,std::vector<double>)>  const & dist, unsigned int dimension,bool store_distance_matrix = false, std::string const & file_name_matrix = "distance_matrix.txt");


private:
    unsigned int size_; // size of the sample P.
    std::vector<double> mat_; // matrix of distances associated to some sample P.
    bool mat_is_sorted_; //=0 when the lines of the matrix mat_ are not sorted yet; =1 when they have been sorted.
    std::vector<double> mat_sorted_;
    static bool smaller2(std::pair<double,unsigned int> a, std::pair<double,unsigned int> b){return a.second<b.second;}   

/**********************************************************************************************************************************************************************
                                                                         SUB-CLASS Vector_dtm
**********************************************************************************************************************************************************************/ 

/**
* @brief Vector_dtm class in Distance_matrix
* @details This is a class which contains two attributes : the mass parameter m_ which is in [0,1] (initialized with -1, and taking this value whenever dtm_ is empty),
* and a vector dtm_ which contains the values of the distance to the uniform measure on the sample (for which mat_ is associated) at all points of the sample,
* with mass parameter m_.
* It contains a method test which returns a p-value and the hypothesis retained : 0 or 1 ; 0 meaning isomorphic, 1 non-isomorphic.
* The method test uses the methods wasserstein_from_pair and smaller. 
**/
    class Vector_dtm {
    public:
    Vector_dtm():dtm_(0), m_(-1) {}
    std::vector<double> dtm_;
    double m_;
    private:
    template<typename T > friend bool smaller (std::pair<double,T> p ,std::pair<double,T> q) {return p.first<q.first;}
    static double wasserstein_from_pair (std::vector<std::pair<double,short>> & pT, const int size);
    public:
    static std::pair<double,bool> test(std::vector<double> & dtm_P,std::vector<double> & dtm_Q, unsigned int n, double alpha, unsigned int Nboot);
    };


// MARC : elle sert à rien cette classe ?  Il faudrait virer ses méthodes qui n'ont rien à voir avec les attribus ???


/**
* vect_dtm_ is an element of type Vector_dtm. It stores a vector dtm_ for some mass parameter m_. Thus, while the mass parameter considered is m,
* there is no need to compute the vector of dtm again.
**/
    Vector_dtm vect_dtm_;

public:

/**
* init_dtm computes and stores the dtm with parameter m, on all points of the sample, in vect_dtm_. 
**/
    void init_dtm(double m);

/**
* The method m_choice computes \f$\frac{W1(d_{\mathbf{1}_P,m}(\mathbf{1}_P),d_{\mathbf{1}_Q,m}(\mathbf{1}_Q))}{\max{Ed_{\mathbf{1}_P,m}(\mathbf{1}_P),Ed_{\mathbf{1}_Q,m}(\mathbf{1}_Q)}}\f$ for all values of m in vector_m.
* Here \f$\mathbf{1}_P\f$ and \f$\mathbf{1}_Q\f$ are empirical measures on the samples associated to the matrices \f$D_P\f$ and \f$D_Q\f$. 
* The best m to use for the test may probably be such that this quantity is maximal or is a local max. But m should not be too small.
**/  
    friend std::vector<double>  m_choice(Distance_matrix & D_P, Distance_matrix & D_Q, std::vector<double> const & vector_m, bool plot_signatures = false, std::vector<bool> const & signatures_to_plot = {});

/**
* The method plot_signatures plots the cumulative distribution functions associated to \f$d_{\mathbf{1}_P,m}(\mathbf{1}_P)\f$ (for \f$D_P\f$) and \f$d_{\mathbf{1}_Q,m}(\mathbf{1}_Q)\f$ (for \f$D_Q\f$).
**/
    friend void plot_signatures(Distance_matrix & D_P, Distance_matrix & D_Q, double m);

/**
* The method test is probably the most important method since it takes in argument the two matrices of distances and returns a p-value
* and the hypothesis retained : 0 for isomorphism, 1 for non-isomorphism.
* It depends on the mass parameter m in [0,1], n which should be much smaller than the sample size, alpha for the level of the test and Nboot,
* the numbers of computations of the random variables from the bootstrap distribution in order to approximate this bootstrap distribution.
*
* Uses the method test of the class Vector_dtm.
**/
    friend std::pair<double,bool> test(Distance_matrix & DP, Distance_matrix & DQ, double m, unsigned int n, double alpha=0.05, unsigned int Nboot=1000);

/**
* The method check_test_error gives an approximative upper-bound for the error of type I, for the parameters m and n, for all pairs (m,n).
* The method consists in spitting the sample associated to mat_ in two parts : two new samples on which we can make the test for the same parameters m and n.
* We do many independent splitting, and compute the mean number of times H_1 is retained, this is an approximation of the type I error.
* Note that for the test to be working, this approximate type I error should be around alpha.
* We prefer to choose n as big as possible but keeping the type I error below alpha.
* If the error is never around alpha, it means that you have chosen wrong parameters, that your data are not large enough or that your data are part of the few examples for which the test cannot work. (as for instance from a uniform measure on a sphere) 
**/
    std::vector<double> check_test_error(std::vector<std::pair<double,unsigned int>> pair_m_n, double alpha, unsigned int Nboot, unsigned int nb_redo);
}; // class Distance_matrix


/***********************************************************************************************************************************************************************
                                
                                      FONCTIONS IN THE CLASS Distance_matrix ---------- TO READ MATRIX OF DISTANCES

************************************************************************************************************************************************************************/

template<typename Matrix_range>
void Distance_matrix::build_from_distance_matrix(Matrix_range const & rmat, bool check_symmetry_and_null_diagonal ){
    this->mat_.clear();
    this->mat_.shrink_to_fit();
    this->mat_is_sorted_=false;
    this->vect_dtm_.m_=-1; // ensure latest vect_dtm_ is deleted if the matrix of distances is redefined.
    auto first = std::begin(rmat);
    auto last = std::end(rmat);
    int size = std::distance(first,last);
    double sqrt_size = std::sqrt(size);
    if(sqrt_size!=int(sqrt_size)){
        std::cerr<<"The size of rmat should be a square, since rmat should be a matrix of size sizeP*sizeP !";
        std::abort();
    }
    this->size_ = sqrt_size;
    if(size>0){
        this->mat_.resize(size);
        int i=0;
        for (auto& d : rmat) {
            this->mat_[i]=d;
            ++i;
        }
        if (check_symmetry_and_null_diagonal){check_symmetry_and_null_diagonal_of_matrix(this->mat_);}
    }
    else{std::cerr<<"Warning : rmat is empty, thus mat_ is still NULL !\n";}
}

template<typename Lower_matrix_range> // If the input rmat is empty, then mat_=[0].
void Distance_matrix::build_from_distance_lower_matrix(Lower_matrix_range const & rmat){
    this->mat_.clear();
    this->mat_.shrink_to_fit();
    this->mat_is_sorted_=false;
    this->vect_dtm_.m_=-1;
    auto first = std::begin(rmat);
    auto last = std::end(rmat);
    int size_rmat = std::distance(first,last);
    int size = int(std::sqrt(2*size_rmat))+1;    
    if(size_rmat != size*(size-1)/2.0){
        std::cerr<<"The size of rmat should be like n(n-1)/2 !";
        std::abort();
    }
    this->size_ = size;
    this->mat_.resize(size*size);
    for(int i=0; i<size; ++i){
        this->mat_[i*size+i]=0;
        for (int j=0; j<i; ++j){
            double temp = *first;
            this->mat_[i*size+j]=temp;
            this->mat_[j*size+i]=temp;
            ++first;
        }
    }
}

void Distance_matrix::build_from_file_matrix(std::string const & file_name, bool check_symmetry_and_null_diagonal){
    std::vector<double> lower_matrix;
    read_file(lower_matrix,file_name);
    build_from_distance_matrix(lower_matrix,check_symmetry_and_null_diagonal);
}



void Distance_matrix::build_from_file_lower_matrix(std::string const & file_name){
    std::vector<double> lower_matrix;
    read_file(lower_matrix,file_name);
    build_from_distance_lower_matrix(lower_matrix);
}


void Distance_matrix::build_from_file_point_set(std::string const & file_name, std::function<double (std::vector<double>,std::vector<double>)> const & dist, unsigned int dimension, bool store_distance_matrix, std::string const & file_name_matrix){

    std::vector<double> point_cloud(0);

    this->mat_.clear();
    this->mat_.shrink_to_fit();
    this->mat_is_sorted_=false;
    this->vect_dtm_.m_=-1;

    read_file(point_cloud,file_name);

    unsigned int size_file = point_cloud.size();
    if(size_file==0){
        std::cerr<<"Warning : There is no point in the file, thus mat_ is still NULL !\n";
    }
    unsigned int size = int(double(size_file)/dimension); 
    if(size_file != size*dimension){
        std::cerr<<"The size of the file should be like nxd !";
        std::abort();
    }
    this->mat_.resize(size*size);
    this->size_ = size;

    std::vector<double> point_i(dimension);
    std::vector<double> point_j(dimension);
    for(unsigned int i=0; i<size; ++i){
        this->mat_[i*size+i]=0;
        for (unsigned int k=0; k<dimension; ++k) {
            point_i[k]=point_cloud[i*dimension+k];
        }
        for (unsigned int j=0; j<i; ++j){
            for (unsigned int k=0; k<dimension; ++k){
                point_j[k]=point_cloud[j*dimension+k];
            }
            double distance_value = dist(point_i,point_j);
            this->mat_[i*size+j]=distance_value;
            this->mat_[j*size+i]=distance_value;
        }
    }
    if(store_distance_matrix){write_file(this->mat_,file_name_matrix);}
}

template<typename P_range,typename Dist>
void Distance_matrix::build_from_point_set(P_range const & rP, Dist const & dist, bool store_distance_matrix, std::string const & file_name){
    this->mat_.clear();
    this->mat_.shrink_to_fit();
    this->mat_is_sorted_=false;
    this->vect_dtm_.m_=-1;
    auto first = std::begin(rP);
    auto last = std::end(rP);
    auto it_i = first;
    auto it_j = first;
    int size = std::distance(first,last);
    this->size_ = size;
    if(size>0){
        this->mat_.resize(size*size);
        for(int i=0; i<size; ++i){
            it_j = first;
            for (int j=0; j<i; ++j){
                double temp = dist(*it_i,*it_j);
                this->mat_[i*size+j]=temp;
                this->mat_[j*size+i]=temp;
                ++it_j;
            }
            this->mat_[i*size+i]=0;
            ++it_i; 
        }
        if(store_distance_matrix){write_file(this->mat_,file_name);}
    }
    else{std::cerr<<"Warning : rmat is empty, thus mat_ is still NULL !\n";}
}



/***********************************************************************************************************************************************************************
                                
                                                                  FONCTIONS IN THE CLASS Vector_dtm 

************************************************************************************************************************************************************************/


/******************************************************           Wasserstein            *****************************************************************************/



/// Compute the Wasserstein distance from a table of pairs with data from the same size.


double Distance_matrix::Vector_dtm::wasserstein_from_pair (std::vector<std::pair<double,short>> & pT, const int size){
    std::sort(pT.begin(),pT.end(),smaller<short>);
    short C = pT[0].second;
    double tmin = pT[0].first;
    double tmax = pT[1].first;
    double W = (tmax - tmin)*abs(C);
    for (int i = 1 ; i < 2*size-1 ; ++i){
        tmin = tmax ;
        tmax = pT[i+1].first;
        C += pT[i].second;
        W += (tmax - tmin)*abs(C);
    }
    return W/size;
}

/****************************************************       Test       ************************************************************************************************/


std::pair<double,bool> Distance_matrix::Vector_dtm::test(std::vector<double> & dtm_P,std::vector<double> & dtm_Q, unsigned int n, double alpha, unsigned int Nboot){
    
    std::random_device rd;  
    std::mt19937 rng(rd()); 

    unsigned int sizeP = dtm_P.size();
    unsigned int sizeQ = dtm_Q.size();
    if(n>sizeP || n>sizeQ){std::cerr<<" n should be smaller than N."<<std::endl; std::abort();}

    for (unsigned int i=sizeP; i>sizeP-n; --i){
        std::uniform_int_distribution<int> uni(0,i-1); 
        int a = uni(rng);
        std::swap(dtm_P[a],dtm_P[i-1]);
    }
    for (unsigned int i=sizeQ; i>sizeQ-n; --i){
        std::uniform_int_distribution<int> uni(0,i-1);
        int a = uni(rng);
        std::swap(dtm_Q[a],dtm_Q[i-1]);
    }
    std::vector<std::pair<double,short> > dtm_by_subset (2*n);
    for (unsigned int l = 0 ; l < n ; ++l){
        dtm_by_subset[l].first = dtm_P[sizeP-l-1];
        dtm_by_subset[l].second = 1;  
    }
    for (unsigned int l = 0 ; l < n ; ++l){
        dtm_by_subset[l+n].first = dtm_Q[sizeQ-l-1];
        dtm_by_subset[l+n].second = -1;  
    }

    double stat = wasserstein_from_pair (dtm_by_subset,n);

  
    std::uniform_int_distribution<int> uniP(0,sizeP-1);
    std::uniform_int_distribution<int> uniQ(0,sizeQ-1);
    std::vector<double> vect_boot (Nboot); // Nboot realisations of the bootstrap law.
    for (unsigned int i=0 ; i<(Nboot/2) ; ++i){
        for (unsigned int l=0 ; l<n ; ++l){
            int a = uniP(rng);
            dtm_by_subset[l].first = dtm_P[a];
            dtm_by_subset[l].second = 1;
        }
        for (unsigned int l=n ; l<2*n ; ++l){
            int a = uniP(rng);
            dtm_by_subset[l].first = dtm_P[a];
            dtm_by_subset[l].second = -1;
        }
        vect_boot[i] = wasserstein_from_pair(dtm_by_subset,n);
    }
    for (unsigned int i=(Nboot/2) ; i<Nboot ; ++i){
        for (unsigned int l=0 ; l<n ; ++l){
            int a = uniQ(rng);
            dtm_by_subset[l].first = dtm_Q[a];
            dtm_by_subset[l].second = 1;
        }
        for (unsigned int l=n ; l<2*n ; ++l){
            int a = uniQ(rng);
            dtm_by_subset[l].first = dtm_Q[a];
            dtm_by_subset[l].second = -1;
        }
        vect_boot[i] = wasserstein_from_pair(dtm_by_subset,n);
    }

    double temp = 0;
    for(unsigned int i=0;i<Nboot;++i){
        temp += (vect_boot[i]>=stat);
    }
    double p_value = temp/Nboot;
    bool b = (p_value<=alpha);
    std::pair<double,bool> pair = {p_value,b};

return pair;
}


/***********************************************************************************************************************************************************************
                                
                                                FONCTIONS IN THE CLASS Distance_matrix ------------ FOR THE DTM

************************************************************************************************************************************************************************/


// From a distance matrix, build the vector of size "size" with inside the DTM of parameter m.
void Distance_matrix::init_dtm (double m){

    this->vect_dtm_.m_=m;
    unsigned int size = this->size_;
    if (size!=this->vect_dtm_.dtm_.size()){
        this->vect_dtm_.dtm_.clear();
        this->vect_dtm_.dtm_.shrink_to_fit();
        this->vect_dtm_.dtm_.resize(size);
    }
    double* pDTM = this->vect_dtm_.dtm_.data();
    unsigned int k = size*m;
    if((k!=0) && (m!=1)){
        if (this->mat_is_sorted_){
            double* pD = this->mat_sorted_.data();
            for (unsigned int j = 0 ; j < size ; ++j){
                double s = 0;
                for (unsigned int i = 0 ; i < k ; ++i){
                    s += *(pD+j*size+i);
                }
                *(pDTM+j) = (s+(m*size-k)**(pD+j*size+k))/(m*size);
            }
        }
        else{
            mat_sorted_.resize(size*size);
            double* pD = this->mat_sorted_.data();
            for (unsigned int j = 0 ; j < size ; ++j){
                for (unsigned int l = 0 ; l < size ; ++l){ 
                    this->mat_sorted_[j*size + l] = this->mat_[j*size + l];
                }
                std::sort(pD+j*size,pD+(j+1)*size);
                double s = 0;
                for (unsigned int i = 0 ; i < k ; ++i){
                    s += *(pD+j*size+i);
                }
                *(pDTM+j) = (s+(m*size-k)**(pD+j*size+k))/(m*size);
            }
            this->mat_is_sorted_=true;
        }
    }
    else{
        if (k==0){
            for (unsigned int j = 0 ; j < size ; ++j){
                *(pDTM+j) = 0;
            }
        }
        else{
            double* pD = this->mat_.data();
            for (unsigned int j = 0 ; j < size ; ++j){
                double s = 0;
                for (unsigned int i = 0 ; i < size ; ++i){
                    s += *(pD+j*size+i);
                }
                *(pDTM+j) = s/size;
            }
        }
    }
}


/*******************************************************      Test        *******************************************************************************************/

std::pair<double,bool> test(Distance_matrix & DP, Distance_matrix & DQ, double m, unsigned int n, double alpha, unsigned int Nboot){
    if(DP.size_<=0){
        std::cerr<<"\n The number of points in the set P should be positive.\n"<<std::endl;
        std::abort();
    }
    if(DQ.size_<=0){
        std::cerr<<"\n The number of points in the set Q should be positive.\n"<<std::endl;
        std::abort();
    }    
    if(m<0 || m>1){
        std::cerr<<"\n m should belong to the interval [0,1].\n"<<std::endl;
        std::abort();
    }
    if(n<=0 || n>DP.size_ || n>DQ.size_){
        std::cerr<<"\n n should be an integer in {1,2,...,min("<<DP.size_<<","<<DQ.size_<<")}.\n"<<std::endl;
        std::abort();
    }
    if(alpha<=0 || alpha>1){
        std::cerr<<"\n alpha should belong to the interval (0,1].\n"<<std::endl;
        std::abort();
    }
    if(Nboot<=0){
        std::cerr<<"\n Nboot should be positive.\n"<<std::endl;
        std::abort();
    }

    if (DP.vect_dtm_.m_!=m){DP.init_dtm(m);}
    if (DQ.vect_dtm_.m_!=m){DQ.init_dtm(m);}

    std::pair<double,bool>  pair = DP.vect_dtm_.test(DP.vect_dtm_.dtm_, DQ.vect_dtm_.dtm_, n, alpha, Nboot);
    return pair;  
}

/******************************************************** Choice of m ***********************************************************************************************/

// If plot_signatures=1, the vector signatures_to_plot should be of the same size than vector_m. We will plot the signatures for all m=vector_m[i] such that signatures_to_plot[i] = 1.
std::vector<double>  m_choice(Distance_matrix & D_P, Distance_matrix & D_Q, std::vector<double> const & vector_m, bool plot_signatures, std::vector<bool> const & signatures_to_plot){

    std::vector<std::pair<double,bool>> vect_m_plot(vector_m.size());
    if(plot_signatures == 0){
        for(int i = 0 ; i<vector_m.size() ; ++i){
            vect_m_plot[i]={vector_m[i],false};
        }
    }
    else{
        if(signatures_to_plot.size()!=vector_m.size()){std::cerr<<"signatures_to_plot and vector_m should be of the same size"; std::abort;}
        else{
            for(int i = 0 ; i<vector_m.size() ; ++i){
                vect_m_plot[i]={vector_m[i],signatures_to_plot[i]};
            }
        }
    }

    std::sort(vect_m_plot.begin(),vect_m_plot.end(),smaller<bool>);

    std::vector<double> m_loss(vector_m.size());

    unsigned int size_P = D_P.size_;
    std::vector<double> cumsum_P(size_P);

    if(D_P.mat_is_sorted_==false){
        D_P.mat_sorted_.resize(size_P*size_P);
        for (unsigned int j = 0 ; j < size_P ; ++j){
            for (unsigned int l = 0 ; l < size_P ; ++l){ 
                D_P.mat_sorted_[j*size_P + l] = D_P.mat_[j*size_P + l];
            }
            std::sort(D_P.mat_sorted_.data() + j*size_P,D_P.mat_sorted_.data() + (j+1)*size_P);
        }
        D_P.mat_is_sorted_ = true;
    }
    D_P.vect_dtm_.dtm_.clear();
    D_P.vect_dtm_.dtm_.shrink_to_fit();
    D_P.vect_dtm_.dtm_.resize(size_P);
    double* pDTM_P = D_P.vect_dtm_.dtm_.data();
    unsigned int temp_k_P = 0;
    unsigned int k_P;

    unsigned int size_Q = D_Q.size_;
    std::vector<double> cumsum_Q(size_Q);

    if(D_Q.mat_is_sorted_==false){
        D_Q.mat_sorted_.resize(size_Q*size_Q);
        for (unsigned int j = 0 ; j < size_Q ; ++j){
            for (unsigned int l = 0 ; l < size_Q ; ++l){ 
                D_Q.mat_sorted_[j*size_Q + l] = D_Q.mat_[j*size_Q + l];
            }
            std::sort(D_Q.mat_sorted_.data() + j*size_Q,D_Q.mat_sorted_.data() + (j+1)*size_Q);
        }
        D_Q.mat_is_sorted_ = true;
    }
    D_Q.vect_dtm_.dtm_.clear();
    D_Q.vect_dtm_.dtm_.shrink_to_fit();
    D_Q.vect_dtm_.dtm_.resize(size_Q);
    double* pDTM_Q = D_Q.vect_dtm_.dtm_.data();
    unsigned int temp_k_Q = 0;
    unsigned int k_Q;

    std::vector<std::pair<double,short>> vect_pair(size_P+size_Q);

    double* pD_P = D_P.mat_sorted_.data();
    double* pD_Q = D_Q.mat_sorted_.data();

    for(unsigned int it_m = 0; it_m<vector_m.size(); ++it_m){

        double m = vect_m_plot[it_m].first;
        D_P.vect_dtm_.m_ = m;
        D_Q.vect_dtm_.m_ = m;

        if(m!=1){
            k_P = size_P*m;
            k_Q = size_Q*m;

            for (unsigned int j = 0 ; j < size_P ; ++j){
                double s = cumsum_P[j];
                for (unsigned int i = temp_k_P ; i < k_P ; ++i){
                    s += *(pD_P+j*size_P+i);
                }
                cumsum_P[j] = s;
                *(pDTM_P+j) = (s+(m*size_P-k_P)**(pD_P+j*size_P+k_P))/(m*size_P);
            }
            
            for (unsigned int j = 0 ; j < size_Q ; ++j){
                double s = cumsum_Q[j];
                for (unsigned int i = temp_k_Q ; i < k_Q ; ++i){
                    s += *(pD_Q+j*size_Q+i);
                }
                cumsum_Q[j] = s;
                *(pDTM_Q+j) = (s+(m*size_Q-k_Q)**(pD_Q+j*size_Q+k_Q))/(m*size_Q);
            }
        }
     
        else{
            k_P = size_P*m;
            k_Q = size_Q*m;

            for (unsigned int j = 0 ; j < size_P ; ++j){
                double s = cumsum_P[j];
                for (unsigned int i = temp_k_P ; i < k_P ; ++i){
                    s += *(pD_P+j*size_P+i);
                }
                cumsum_P[j] = s;
                *(pDTM_P+j) = s/(m*size_P);
            }
            
            for (unsigned int j = 0 ; j < size_Q ; ++j){
                double s = cumsum_Q[j];
                for (unsigned int i = temp_k_Q ; i < k_Q ; ++i){
                    s += *(pD_Q+j*size_Q+i);
                }
                cumsum_Q[j] = s;
                *(pDTM_Q+j) = s/(m*size_Q);
            }
        }

        if(vect_m_plot[it_m].second==true){
            std::vector<double> sorted_dtm_P(size_P);
            std::vector<double> stairs_P(size_P);
            for(unsigned int j = 0; j<size_P; ++j){sorted_dtm_P[j] = *(pDTM_P+j);}
            std::sort(sorted_dtm_P.begin(),sorted_dtm_P.end());
            for(unsigned int j = 0; j<size_P; ++j){stairs_P[j] = j/double(size_P);}

            std::vector<double> sorted_dtm_Q(size_Q);
            std::vector<double> stairs_Q(size_Q);
            for(unsigned int j = 0; j<size_Q; ++j){sorted_dtm_Q[j] = *(pDTM_Q+j);}
            std::sort(sorted_dtm_Q.begin(),sorted_dtm_Q.end());
            for(unsigned int j = 0; j<size_Q; ++j){stairs_Q[j] = j/double(size_Q);}

            double xmin = std::min(sorted_dtm_P[0],sorted_dtm_Q[0]);
            double xmax = std::max(sorted_dtm_P[size_P-1],sorted_dtm_Q[size_Q-1]);
            Gnuplot gp;
            gp << "set title 'Empirical signatures, for m = "<<std::setprecision(5)<<m<<"'\n";
            gp << "set yrange [0:1]\n";
            gp << "set xrange ["<<xmin<<":"<<xmax<<"]\n";
            gp << "plot '-' with steps tit 'sample P', '-' with steps tit 'sample Q'\n ";
            gp.send1d(boost::make_tuple(sorted_dtm_P,stairs_P));
            gp.send1d(boost::make_tuple(sorted_dtm_Q,stairs_Q));
        }

    // Wasserstein between the dtm


        for (unsigned int i = 0 ; i < size_P ; ++i){
            vect_pair[i].first = *(pDTM_P+i);
            vect_pair[i].second = size_Q;  
        }
        for (unsigned int i = 0 ; i < size_Q ; ++i){
            vect_pair[i+size_P].first = *(pDTM_Q+i);
            vect_pair[i+size_P].second = -size_P;  
        }

        std::sort(vect_pair.begin(),vect_pair.end(),smaller<int>);
        double C = vect_pair[0].second;
        double tmin = vect_pair[0].first;
        double tmax = vect_pair[1].first;
        double W = (tmax - tmin)*abs(C);
        for (unsigned int i = 1 ; i < size_P+size_Q-1 ; ++i){
            tmin = tmax ;
            tmax = vect_pair[i+1].first;
            C += vect_pair[i].second;
            W += (tmax - tmin)*abs(C);
        }
        double w = W/(size_P*size_Q);

        double s = 0;
        for(unsigned int j = 0 ; j<size_P ; ++j){s +=*(pDTM_P+j);}
        double E_P = s/size_P;// mean(DP.vect_dtm_.dtm_)

        s = 0;
        for(unsigned int j = 0 ; j<size_Q ; ++j){s +=*(pDTM_Q+j);}
        double E_Q = s/size_Q;// mean(DQ.vect_dtm_.dtm_)

        m_loss[it_m] = w*w/std::max(E_P,E_Q);

        temp_k_P = k_P;
        temp_k_Q = k_Q;
    }
    Gnuplot gp;
    std::vector<std::array<double,2>> cost_plot(vector_m.size());
    for(int i = 0; i<vector_m.size(); ++i){cost_plot[i][0] = vector_m[i]; cost_plot[i][1] = m_loss[i];}
    gp << "plot [0:1] '-' with lines tit 'cost'\n";
    gp.send1d(cost_plot);
    return m_loss;
}


/**************************************************** Plot signatures *******************************************************************************************/

void plot_signatures(Distance_matrix & D_P, Distance_matrix & D_Q, double m){
    unsigned int size_P = D_P.size_;
    std::vector<double> cumsum_P(size_P);

    if(D_P.mat_is_sorted_==false){
        D_P.mat_sorted_.resize(size_P*size_P);
        for (unsigned int j = 0 ; j < size_P ; ++j){
            for (unsigned int l = 0 ; l < size_P ; ++l){ 
                D_P.mat_sorted_[j*size_P + l] = D_P.mat_[j*size_P + l];
            }
            std::sort(D_P.mat_sorted_.data() + j*size_P,D_P.mat_sorted_.data() + (j+1)*size_P);
        }
        D_P.mat_is_sorted_ = true;
    }
    D_P.vect_dtm_.dtm_.clear();
    D_P.vect_dtm_.dtm_.shrink_to_fit();
    D_P.vect_dtm_.dtm_.resize(size_P);
    double* pDTM_P = D_P.vect_dtm_.dtm_.data();
    unsigned int temp_k_P = 0;
    unsigned int k_P;

    unsigned int size_Q = D_Q.size_;
    std::vector<double> cumsum_Q(size_Q);

    if(D_Q.mat_is_sorted_==false){
        D_Q.mat_sorted_.resize(size_Q*size_Q);
        for (unsigned int j = 0 ; j < size_Q ; ++j){
            for (unsigned int l = 0 ; l < size_Q ; ++l){ 
                D_Q.mat_sorted_[j*size_Q + l] = D_Q.mat_[j*size_Q + l];
            }
            std::sort(D_Q.mat_sorted_.data() + j*size_Q,D_Q.mat_sorted_.data() + (j+1)*size_Q);
        }
        D_Q.mat_is_sorted_ = true;
    }
    D_Q.vect_dtm_.dtm_.clear();
    D_Q.vect_dtm_.dtm_.shrink_to_fit();
    D_Q.vect_dtm_.dtm_.resize(size_Q);
    double* pDTM_Q = D_Q.vect_dtm_.dtm_.data();
    unsigned int temp_k_Q = 0;
    unsigned int k_Q;

    double* pD_P = D_P.mat_sorted_.data();
    double* pD_Q = D_Q.mat_sorted_.data();

    D_P.vect_dtm_.m_ = m;
    D_Q.vect_dtm_.m_ = m;

    if(m!=1){
        k_P = size_P*m;
        k_Q = size_Q*m;

        for (unsigned int j = 0 ; j < size_P ; ++j){
            double s = cumsum_P[j];
            for (unsigned int i = temp_k_P ; i < k_P ; ++i){
                s += *(pD_P+j*size_P+i);
            }
            cumsum_P[j] = s;
            *(pDTM_P+j) = (s+(m*size_P-k_P)**(pD_P+j*size_P+k_P))/(m*size_P);
        }
        
        for (unsigned int j = 0 ; j < size_Q ; ++j){
            double s = cumsum_Q[j];
            for (unsigned int i = temp_k_Q ; i < k_Q ; ++i){
                s += *(pD_Q+j*size_Q+i);
            }
            cumsum_Q[j] = s;
            *(pDTM_Q+j) = (s+(m*size_Q-k_Q)**(pD_Q+j*size_Q+k_Q))/(m*size_Q);
        }
    }
     
    else{
        k_P = size_P*m;
        k_Q = size_Q*m;
        for (unsigned int j = 0 ; j < size_P ; ++j){
            double s = cumsum_P[j];
            for (unsigned int i = temp_k_P ; i < k_P ; ++i){
                s += *(pD_P+j*size_P+i);
            }
            cumsum_P[j] = s;
            *(pDTM_P+j) = s/(m*size_P);
        }
            
        for (unsigned int j = 0 ; j < size_Q ; ++j){
            double s = cumsum_Q[j];
            for (unsigned int i = temp_k_Q ; i < k_Q ; ++i){
                s += *(pD_Q+j*size_Q+i);
            }
            cumsum_Q[j] = s;
            *(pDTM_Q+j) = s/(m*size_Q);
        }
    }

    std::vector<double> sorted_dtm_P(size_P);
    std::vector<double> stairs_P(size_P);
    for(unsigned int j = 0; j<size_P; ++j){sorted_dtm_P[j] = *(pDTM_P+j);}
    std::sort(sorted_dtm_P.begin(),sorted_dtm_P.end());
    for(unsigned int j = 0; j<size_P; ++j){stairs_P[j] = j/double(size_P);}

    std::vector<double> sorted_dtm_Q(size_Q);
    std::vector<double> stairs_Q(size_Q);
    for(unsigned int j = 0; j<size_Q; ++j){sorted_dtm_Q[j] = *(pDTM_Q+j);}
    std::sort(sorted_dtm_Q.begin(),sorted_dtm_Q.end());
    for(unsigned int j = 0; j<size_Q; ++j){stairs_Q[j] = j/double(size_Q);}

    double xmin = std::min(sorted_dtm_P[0],sorted_dtm_Q[0]);
    double xmax = std::max(sorted_dtm_P[size_P-1],sorted_dtm_Q[size_Q-1]);
    Gnuplot gp;
    gp << "set title 'Empirical signatures, for m = "<<std::setprecision(5)<<m<<"'\n";
    gp << "set yrange [0:1]\n";
    gp << "set xrange ["<<xmin<<":"<<xmax<<"]\n";
    gp << "plot '-' with steps tit 'sample P', '-' with steps tit 'sample Q'\n ";
    gp.send1d(boost::make_tuple(sorted_dtm_P,stairs_P));
    gp.send1d(boost::make_tuple(sorted_dtm_Q,stairs_Q));
}

/***************************************************** Check test error *****************************************************************************************/

std::vector<double> Distance_matrix::check_test_error(std::vector<std::pair<double,unsigned int>> pair_m_n, double alpha, unsigned int Nboot, unsigned int nb_redo){
    // On split la matrice this->mat_ en deux nb_redo fois. Attention, il faut qu'elle n'ai pas été triée avant !!!

    unsigned int size  = this->size_;

    //for(unsigned int ll=0;ll<size*size;++ll){std::cout<<mat_[ll]<<" ";}

    unsigned int size_P = size/2;
    unsigned int size_Q = size-size/2;

    Distance_matrix D_P;
    Distance_matrix D_Q;

    D_P.size_ = size_P;
    D_Q.size_ = size_Q;

    D_P.mat_.clear();
    D_P.mat_.shrink_to_fit();
    D_P.mat_.resize(size_P*size_P);
    double* pmat_P = D_P.mat_.data();
    D_P.vect_dtm_.dtm_.clear();
    D_P.vect_dtm_.dtm_.shrink_to_fit();
    D_P.vect_dtm_.dtm_.resize(size_P);

    D_Q.mat_.clear();
    D_Q.mat_.shrink_to_fit();
    D_Q.mat_.resize(size_Q*size_Q);
    double* pmat_Q = D_Q.mat_.data();
    D_Q.vect_dtm_.dtm_.clear();
    D_Q.vect_dtm_.dtm_.shrink_to_fit();
    D_Q.vect_dtm_.dtm_.resize(size_Q);

    double* pmat = this->mat_.data();

    // sort pair_m_n with m increasing.
    std::sort(pair_m_n.begin(),pair_m_n.end(),smaller<unsigned int>);

    // preprocessing
    if(size_P<=0){
        std::cerr<<"\n The number of points in the set P should be positive.\n"<<std::endl;
        std::abort();
    }
    if(size_Q<=0){
        std::cerr<<"\n The number of points in the set Q should be positive.\n"<<std::endl;
        std::abort();
    }    
    if(alpha<=0 || alpha>1){
        std::cerr<<"\n alpha should belong to the interval (0,1].\n"<<std::endl;
        std::abort();
    }
    if(Nboot<=0){
        std::cerr<<"\n Nboot should be positive.\n"<<std::endl;
        std::abort();
    }
    if(nb_redo<=0){
        std::cerr<<"\n nb_redo should be positive.\n"<<std::endl;
        std::abort();
    }
    for(unsigned int i=0;i<pair_m_n.size(); ++i){
        if(pair_m_n[i].first<0 || pair_m_n[i].first>1){
            std::cerr<<"\n m should belong to the interval [0,1].\n"<<std::endl;
            std::abort();
        }
        if(pair_m_n[i].second<=0 || pair_m_n[i].second>size_P || pair_m_n[i].second>size_Q){
            std::cerr<<"\n n should be an integer in {1,2,...,min("<<size_P<<","<<size_Q<<")}.\n"<<std::endl;
            std::abort();
        }
    }


    std::vector<int> permutation(size);
    for(unsigned int i = 0; i<size ; ++i){permutation[i]=i;}

    std::pair<double,bool> test_result;
    std::vector<double> nb_rejected_tests(pair_m_n.size());

    for(unsigned int redo=0; redo<nb_redo ; ++redo){
        // split the matrix in two matrices of distances.
        std::random_shuffle(permutation.begin(),permutation.end());

        //D_P.vect_dtm_.dtm_.clear();
        //D_Q.vect_dtm_.dtm_.clear();

        //std::cout<<"permutation";
        //for(unsigned int ll=0;ll<permutation.size();++ll){std::cout<<permutation[ll]<<" ";}
        //std::cout<<std::endl;

        for(unsigned int i = 0; i<size_P ; ++i){
            pmat_P[i*size_P+i] = 0;
            unsigned int ind_i = permutation[i];
            for(unsigned int j = 0; j<i ; ++j){
                double temp = pmat[ind_i*size+permutation[j]];
                pmat_P[i*size_P+j] = temp;
                pmat_P[j*size_P+i] = temp;
            }
        }

        //for(unsigned int ll=0;ll<size_P*size_P;++ll){std::cout<<D_P.mat_[ll]<<" ";}

        for(unsigned int i = 0; i<size_Q ; ++i){
            pmat_Q[i*size_Q+i] = 0;
            unsigned int ind_i = permutation[size_P+i];
            for(unsigned int j = 0; j<i ; ++j){
                double temp = pmat[ind_i*size+permutation[size_P+j]];
                pmat_Q[i*size_Q+j] = temp;
                pmat_Q[j*size_Q+i] = temp;
            }
        }

        //for(unsigned int ll=0;ll<size_P*size_P;++ll){std::cout<<D_Q.mat_[ll]<<" ";}
        
        // initialise vect_dtm_.m_ for each matrix.
        D_P.vect_dtm_.m_=-1;
        D_Q.vect_dtm_.m_=-1;

        // the matrices are not sorted yet.
        D_P.mat_is_sorted_ = false;
        D_Q.mat_is_sorted_ = false;

        // compute the test result
        for(unsigned int i=0;i<pair_m_n.size();++i){
            double m = pair_m_n[i].first;
            unsigned int n = pair_m_n[i].second;
            if (D_P.vect_dtm_.m_!=m){D_P.init_dtm(m);}
            if (D_Q.vect_dtm_.m_!=m){D_Q.init_dtm(m);}

            //std::cout<<"dtmP : ";
            //for(unsigned int ll=0;ll<D_P.size_;++ll){std::cout<<" "<<D_P.vect_dtm_.dtm_[ll];}
            //std::cout<<std::endl;

            //std::cout<<"dtmQ : ";
            //for(unsigned int ll=0;ll<D_Q.size_;++ll){std::cout<<" "<<D_Q.vect_dtm_.dtm_[ll];}
            //std::cout<<std::endl;

            test_result = Distance_matrix::Vector_dtm::test(D_P.vect_dtm_.dtm_, D_Q.vect_dtm_.dtm_, n, alpha, Nboot);

            nb_rejected_tests[i] += test_result.second;
            //std::cout<<" is_rej = "<<test_result.second;
            //std::cout<<std::endl;
            //std::cout<<std::endl;                   
        }
    }
    for(auto& d : nb_rejected_tests){d /= nb_redo;}

    Gnuplot gp;
    std::vector<std::array<double,2>> error_plot(pair_m_n.size());
    for(int i = 0; i<pair_m_n.size(); ++i){error_plot[i][0] = pair_m_n[i].second; error_plot[i][1] = nb_rejected_tests[i];} 
    auto minn = std::min_element(pair_m_n.begin(),pair_m_n.end(),smaller2);
    auto maxn = std::max_element(pair_m_n.begin(),pair_m_n.end(),smaller2);
    int min_bound = (*minn).second-1;
    int max_bound = (*maxn).second+1;
    gp << "plot ["<<min_bound<<":"<<max_bound<<"] '-' with lines tit 'cost'\n";
    gp.send1d(error_plot);

    return nb_rejected_tests;

}


}//isomorphism_test
}// Gudhi

#endif  // ISOMORPHISM_TEST_H_
