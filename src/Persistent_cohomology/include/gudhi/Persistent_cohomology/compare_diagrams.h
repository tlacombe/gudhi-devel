
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cstdio>
#include <cstring>

#include <boost/tuple/tuple.hpp>


/** Compare two intervals by length.*/
struct cmp_bares_lex {
  bool operator() (  boost::tuple< int, double, double > & bare1
                  ,  boost::tuple< int, double, double > & bare2 )
  {
    if(bare1.get<1>() != bare2.get<1>()) return bare1.get<1>() < bare2.get<1>();
    if(bare1.get<2>() != bare2.get<2>()) return bare1.get<2>() < bare2.get<2>();
    if(bare1.get<0>() != bare2.get<0>()) return bare1.get<0>() < bare2.get<0>();

    return false;
  }
};

void
read_diagram ( std::string file_name
             , std::vector< boost::tuple< int, double, double > > & diag)
{  
  std::ifstream in_file ( file_name.c_str(), std::ios::in);
  if(!in_file.is_open()) {
    std::cerr << "Unable to open file " << file_name << std::endl;
    return;}

  std::string line;
  int dim;
  double b,d;
  double x;
  while( getline ( in_file , line ) )
  {
    std::istringstream iss( line );
    iss >> dim; iss >> b; iss >> d;
    boost::tuple< int, double, double > bare(dim,b,d);
    if( iss >> x ) { std::cerr << "More than 3 values in a line. \n"; }
    diag.push_back(bare);
  }
  in_file.close();
}

void compare_diagrams ( std::string d_file1
                      , std::string d_file2 )
{
  std::vector< boost::tuple< int, double, double > > diag1, diag2;
  read_diagram( d_file1, diag1 );
  read_diagram( d_file2, diag2 );

  cmp_bares_lex cmp;

  sort(diag1.begin(),diag1.end(),cmp);
  sort(diag2.begin(),diag2.end(),cmp);

  auto it1 = diag1.begin(); auto it2 = diag2.begin();
  while ( (it1 != diag1.end()) && (it2 != diag2.end()) )
  {
    // std::cout << it1->get<0>() << " "
    //                            << it1->get<1>() << " "
    //                            << it1->get<2>() << "     ";
    // std::cout << it2->get<0>() << " "
    //                            << it2->get<1>() << " "
    //                            << it2->get<2>() << std::endl; 
   
    // ++it1; ++it2;

    if( !(cmp(*it1,*it2)) && !(cmp(*it2,*it1)) ) { ++it1; ++it2; } //==
    else {
      if( cmp(*it1,*it2) ) { std::cout << it1->get<0>() << " "
                                       << it1->get<1>() << " "
                                       << it1->get<2>() << std::endl; 
                             ++it1; }                                   // <
      else                 { std::cout << "                         "
                                       << it2->get<0>() << " "
                                       << it2->get<1>() << " "
                                       << it2->get<2>() << std::endl; 
                             ++it2; }                                    // >
    }
  }
  while(it1 != diag1.end()) { std::cout << it1->get<0>() << " "
                                        << it1->get<1>() << " "
                                        << it1->get<2>() << std::endl; 
                             ++it1; } 
  while(it2 != diag2.end()) { std::cout << "                         "
                                        << it2->get<0>() << " "
                                        << it2->get<1>() << " "
                                        << it2->get<2>() << std::endl; 
                             ++it2; }
}

