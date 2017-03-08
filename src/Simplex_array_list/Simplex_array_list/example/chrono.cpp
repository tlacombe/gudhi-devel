#include <iostream>
#include <random>
#include <chrono>

#include <gudhi/LSAL.h>
#include <gudhi/STW.h>

using namespace Gudhi;

int n = 300;
//int d = 40;

int nb_add1 = 3000;
int nb_membership1 = 4000;
int nb_contraction = 300;
int nb_add2 = 3000;
int nb_membership2 = 400000;

Simplex random_simplex(int n, int d){
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dis(1, n);
    Simplex s;
    while(s.size()!=d)
        s.insert(dis(gen));
    return s;
}

std::vector<Simplex> r_vector_simplices(int n, int max_d, int m){
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dis(1, max_d);
    std::vector<Simplex> v;
    for(int i=0; i<m; i++)
        v.push_back(random_simplex(n,dis(gen)));
    return v;
}

template<typename complex_type>
void chrono(int n, int d){
    complex_type K;
    std::vector<Simplex> simplices_add1 = r_vector_simplices(n,d,nb_add1);
    std::vector<Simplex> simplices_membership1 = r_vector_simplices(n,d,nb_membership1);
    std::vector<Simplex> simplices_add2 = r_vector_simplices(n - 2*nb_contraction,d,nb_add2);
    std::vector<Simplex> simplices_membership2 = r_vector_simplices(n - 2*nb_contraction,d,nb_membership2);
    std::chrono::time_point<std::chrono::system_clock> start, end;

    //start = std::chrono::system_clock::now();
    for(const Simplex& s : simplices_add1)
        K.add(s);
    //end = std::chrono::system_clock::now();
    //std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count() << std::endl;

    //start = std::chrono::system_clock::now();
    for(const Simplex& s : simplices_membership1)
        K.membership(s);
    //end = std::chrono::system_clock::now();
    //std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count() << std::endl;

    start = std::chrono::system_clock::now();
    for(int i = 0; i<=nb_contraction; i++)
        K.contraction(n-2*i,n-2*i-1);
    end = std::chrono::system_clock::now();
    auto c3 = std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count();

    start = std::chrono::system_clock::now();
    for(const Simplex& s : simplices_add2)
        K.add(s);
    end = std::chrono::system_clock::now();
    auto c1 = std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count();

    start = std::chrono::system_clock::now();
    for(const Simplex& s : simplices_membership2)
        K.membership(s);
    end = std::chrono::system_clock::now();
    auto c2 = std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count();

    std::cout << c1 << "\t \t" << c2 << "\t \t" << c3 << "\t \t" << K.size() << std::endl;
}

int main(){
    for(int d=5;d<=40;d+=5){
        std::cout << "d=" << d << " \t  Insertions \t Membership \t   Contractions \t Mem" << std::endl;
        std::cout << "LSAL \t \t";
        chrono<LSAL>(n,d);
        std::cout << "SAL \t \t";
        chrono<SAL>(n,d);
        if(d<=15){
            std::cout << "STW \t \t";
            chrono<STW>(n,d);
        }
        std::cout << std::endl;
    }
}
