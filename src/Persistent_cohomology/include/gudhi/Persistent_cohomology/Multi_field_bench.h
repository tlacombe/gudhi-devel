/*
 *  Multi_field.h
 *  Gudhi
 *
 *  Created by Clément Maria on 15/02/14.
 *  Copyright 2014 INRIA. All rights reserved.
 *
 */

#ifndef GUDHI_MULTI_FIELD_BENCH_H
#define GUDHI_MULTI_FIELD_BENCH_H

#include <iostream>
#include <vector>
#include <gmpxx.h>
/**
  *
  */
class Multi_field_bench {
public:
  typedef mpz_class      Element;

  Multi_field_bench (int num_of_primes)
  { 

std:: cout << "\n \n \n \n INIT \n \n";


    std::string fileprimes = "/Users/cmaria/prime_numbers.txt";
    std::ifstream in_ (fileprimes.c_str(),std::ios::in);
    if(!in_.is_open()) {
      std::cerr << "Unable to open file " << fileprimes << std::endl;
      return;}

    unsigned int curr_prime;
    for(int i=0; i<num_of_primes; ++i) {
      if(!(in_ >> curr_prime)) { std::cerr << "Not enough primes in file \n";}
      primes_.push_back(curr_prime);
    }
    in_.close();

    //set m to primorial(bound_prime)
    mpz_primorial_ui(prod_characteristics_.get_mpz_t(),curr_prime);
    num_primes_ = primes_.size();
    //Uvect_ 
    Element Ui; Element tmp_elem;
    for(auto p : primes_)
    {
      tmp_elem = prod_characteristics_ / p;
      mpz_powm_ui ( tmp_elem.get_mpz_t()
                  , tmp_elem.get_mpz_t()
                  , p - 1
                  , prod_characteristics_.get_mpz_t() );
      Uvect_.push_back(tmp_elem);
    }
    mult_id_all = 0;
    for(int idx = 0; idx < num_primes_; ++idx) 
      { mult_id_all = (mult_id_all + Uvect_[idx]) % prod_characteristics_; }



    std::cout << "Initialization Multi-Field: number of primes = " << primes_.size() << "\n";
    // std::cout << prod_characteristics_ << "    ";
    // for(auto p : primes_) { std::cout << p << " ";}
    //   std::cout << std::endl;


std:: cout << "\n \n ENDINIT \n \n \n \n";

  }



  /** \brief Returns the additive idendity \f$0_{\F}\f$ of the field.*/
  Element additive_identity () { return 0; }
  /** \brief Returns the multiplicative identity \f$1_{\F}\f$ of the field.*/
  Element multiplicative_identity () { return mult_id_all; }// 1 everywhere 

  Element multiplicative_identity (Element Q) 
  {
    Element mult_id = 0;
    for(int idx = 0; idx < num_primes_; ++idx) {
      if( (Q % primes_[idx]) == 0 ) 
        { mult_id = (mult_id + Uvect_[idx]) % prod_characteristics_; }
    }
    return mult_id;
  } 

  /** Returns y * w */
  Element times ( Element y, int w ) { 
    Element tmp = (y*w) % prod_characteristics_;
    if(tmp < 0) return prod_characteristics_ + tmp;
    return tmp;
  }

  void plus_equal(Element & x, Element y) 
  { x += y; x %= prod_characteristics_; }

  /** \brief Returns the characteristic \f$p\f$ of the field.*/
  Element characteristic() { return prod_characteristics_; }

  /** Returns the inverse in the field. Modifies P.*/
  std::pair<Element,Element> inverse ( Element x
                                     , Element QS ) 
  { 
    Element QR;
    mpz_gcd( QR.get_mpz_t(), x.get_mpz_t(), QS.get_mpz_t() ); // QR <- gcd(x,QS) 
    if( QR == QS ) return std::pair<Element,Element>(additive_identity()
                                                    , multiplicative_identity() );   //partial inverse is 0
    Element QT = QS / QR;
    Element inv_qt;
    mpz_invert(inv_qt.get_mpz_t(), x.get_mpz_t(), QT.get_mpz_t());

    return std::pair<Element,Element>(
                (inv_qt * multiplicative_identity(QT)) % prod_characteristics_
                , QT                 );
  }
  /** Returns -x * y.*/
  Element times_minus ( Element x, Element y ) 
  { return prod_characteristics_ - ((x*y)%prod_characteristics_); }

  /** Set x <- x + w * y*/
  void plus_times_equal ( Element & x, Element y, Element w ) 
  { x = (x + w * y) % prod_characteristics_; }

  Element               prod_characteristics_;
  std::vector<int>      primes_;
  std::vector<Element>  Uvect_;
  size_t                num_primes_;
  Element               mult_id_all;
  gmp_randstate_t       state; //setup for generate rand numbers with gmp

};
/**
  *
  */
// class Field_Q {

// };

#endif // GUDHI_MULTI_FIELD_BENCH_H 
