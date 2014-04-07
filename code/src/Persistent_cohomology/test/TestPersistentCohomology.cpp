/*
 *  TestPersistentCohomology.h
 *  Gudhi
 *
 *  Created by Clément Maria on 02/28/14.
 *  Copyright 2013 INRIA. All rights reserved.
 *
 */

#include <ctime>
#include "io.h"
#include "Nearest_neighbors.h"
#include "Simplex_tree.h"
#include "Hasse_complex.h"
#include "Persistent_cohomology.h"
#include "Persistent_cohomology/Field_Zp.h"
#include "Persistent_cohomology/Multi_field.h"
#include "Persistent_cohomology/Multi_field_bench.h"


// typedef Persistent_cohomology< Simplex_tree<> , Field_Zp<2> >       PcoH2;
// typedef Persistent_cohomology< Simplex_tree<> , Field_Zp<3> >       PcoH3;
// typedef Persistent_cohomology< Simplex_tree<> , Multi_field<2,30000> >  PcoH23;


bool test_Klein_bottle()
{

  std::string filecomplex = "../Persistent_cohomology/test/Klein_bottle_complex.txt";
  std::ifstream in_ (filecomplex.c_str(),std::ios::in);
  if(!in_.is_open()) { std::cerr << "Unable to open file " << filecomplex << std::endl; }
 //Construct the Simplex Tree
  Simplex_tree<>  st;     in_ >> st;     in_.close();
  st.initialize_filtration();

  std::string file_hasse = "/Users/cmaria/Klein_hasse.txt";
  std::ofstream out_ (file_hasse.c_str(),std::ios::trunc);
  if(!out_.is_open()) { std::cerr << "Unable to open file " << file_hasse << std::endl; }
//Output the Hasse complex
  st.output_hasse(out_);  out_.close();

  std::ifstream in2_ (file_hasse.c_str(),std::ios::in);
  Hasse_complex<> hcpx;    in2_ >> hcpx;   in2_.close();

  Persistent_cohomology< Hasse_complex<> , Field_Zp<2> > Hpcoh (hcpx);
  Hpcoh.compute_persistent_cohomology(0);
  Hpcoh.output_diagram();

std::cout << "\n \n \n";

  Persistent_cohomology< Simplex_tree<> , Field_Zp<11> > Spcoh (st);
  Spcoh.compute_persistent_cohomology(0);
  Spcoh.output_diagram();

  return true;


  for(int i = 2 ;i < 3; i++)
  {
    st.initialize_filtration();
    Persistent_cohomology< Simplex_tree<> , Field_Zp<11> > pcoh (st);

    // clock_t start, end;
    // start = clock();
    pcoh.compute_persistent_cohomology (1000);
    // end = clock();
    // std::cout << (double)(end-start)/CLOCKS_PER_SEC << " s. \n";
  }

}


int main (int argc, char * const argv[]) 
{
test_Klein_bottle();
return 0;

  // int nb_of_primes = 1000000;

  // std::string file_name = "/Users/cmaria/prime_numbers.txt";
  // std::ofstream out_ (file_name.c_str(),std::ios::trunc);
  // if(!out_.is_open()) { std::cerr << "Unable to open file " << file_name << std::endl; }

  // unsigned int curr_prime = 2;
  // mpz_t tmp_prime; mpz_init_set_ui(tmp_prime,2);

  // int count = 1;
  // while (count <= nb_of_primes)
  // {
  //   ++count;
  //   out_ << curr_prime << " ";
  //   mpz_nextprime(tmp_prime,tmp_prime);
  //   curr_prime = mpz_get_ui(tmp_prime);
  // }
  // out_.close();
  // return 0;







  std::string filepoints = "/Users/cmaria/MF_datasets/Cy8.txt";
  std::string filegraph  = "/Users/cmaria/MF_datasets/Cy8_graph.txt";
  double threshold = 0.5;
  compute_rips_graph(filepoints,filegraph,threshold);
  std::cout << "Graph computed \n";

  // return 0;

  auto g = read_graph(filegraph); 

  Simplex_tree<> st;
  std::cout << "Start ST construction \n";
  st.insert_graph (g); //insert the graph in the simplex tree as 1-skeleton
  st.expansion(3);

  std::cout << " Number of simplices: " << st.num_simplices() << std::endl;

  st.initialize_filtration();
  Persistent_cohomology< Simplex_tree<> , Field_Zp<11> > pcoh (st);

  clock_t start, end;
  start = clock();
  pcoh.compute_persistent_cohomology (1000);
  end = clock();
  std::cout << (double)(end-start)/CLOCKS_PER_SEC << " s. \n";


  return 0;
}
