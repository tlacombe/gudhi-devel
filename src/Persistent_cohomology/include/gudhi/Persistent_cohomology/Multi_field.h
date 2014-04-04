/*
 *  Multi_field.h
 *  Gudhi
 *
 *  Created by Clément Maria on 15/02/14.
 *  Copyright 2014 INRIA. All rights reserved.
 *
 */

#ifndef GUDHI_MULTI_FIELD_H
#define GUDHI_MULTI_FIELD_H

#include <iostream>
#include <vector>
#include <gmpxx.h>
/**
  *
  */
template< unsigned int min_prime = 2
        , unsigned int max_prime = 17 >
class Multi_field {
public:
  typedef mpz_class      Element;

  Multi_field ()
  { 
    if(max_prime<2) 
      { std::cerr << "There is no prime less than " << max_prime << std::endl; }
    if(min_prime > max_prime) 
      { std::cerr << "No prime in ["<<min_prime<<":"<<max_prime<<"]"<<std::endl; }
    //initialize the random generation settings with Mersenne Twister algo
    gmp_randinit_mt(state);
    // fill the list of prime numbers
    unsigned int curr_prime = min_prime;
    mpz_t tmp_prime; mpz_init_set_ui(tmp_prime,min_prime);
    //test if min_prime is prime
    int is_prime = mpz_probab_prime_p(tmp_prime,25); //probabilistic primality test
    while(is_prime == 1) //uncertainty result
    { is_prime = mpz_probab_prime_p(tmp_prime,25); } //test again

    if(is_prime == 0) //min_prime is composite
    {
      mpz_nextprime(tmp_prime,tmp_prime);
      curr_prime = mpz_get_ui(tmp_prime);
    }
    
    
    while (curr_prime <= max_prime)
    {
      primes_.push_back(curr_prime);
      mpz_nextprime(tmp_prime,tmp_prime);
      curr_prime = mpz_get_ui(tmp_prime);
    }
    //set m to primorial(bound_prime)
    mpz_primorial_ui(prod_characteristics_.get_mpz_t(),max_prime);

    num_primes_ = primes_.size();

    //Uvect_ 
    Element Ui;
    Element tmp_elem;
    for(auto p : primes_)
    {
      tmp_elem = prod_characteristics_ / p;
      Element tmp_elem_bis = 10;
      mpz_powm_ui(tmp_elem.get_mpz_t(),tmp_elem_bis.get_mpz_t(),
        p,
        prod_characteristics_.get_mpz_t());
      Uvect_.push_back(tmp_elem);
    }
  }



  /** \brief Returns the additive idendity \f$0_{\F}\f$ of the field.*/
  Element additive_identity () { return 0; }
  /** \brief Returns the multiplicative identity \f$1_{\F}\f$ of the field.*/
  Element multiplicative_identity () { return prod_characteristics_; } 

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
    Element tmp;
    if(w < 0) { tmp = prod_characteristics_ - w; }
    else { tmp = w; }
    return (y*tmp)%prod_characteristics_;
  }

  void plus_equal(Element & x, Element y) 
  { x += y; x %= prod_characteristics_; }

  /** \brief Returns the characteristic \f$p\f$ of the field.*/
  Element characteristic() { return prod_characteristics_; }

  /** Returns the inverse in the field. Modifies P.*/
  Element inverse ( Element x
                  , Element & QS ) 
  { 
    Element QR;
    mpz_gcd( QR.get_mpz_t(), x.get_mpz_t(), QS.get_mpz_t() ); // QR <- gcd(x,QS) 
    if( QR == QS ) return 0;   //partial inverse is 0

    Element QT = QS / QR;
    Element inv_qt;
    mpz_invert(inv_qt.get_mpz_t(), x.get_mpz_t(), QT.get_mpz_t());

    QS = QR;
    return (inv_qt * multiplicative_identity(QT)) % prod_characteristics_;
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
  gmp_randstate_t       state; //setup for generate rand numbers with gmp

};
/**
  *
  */
// class Field_Q {

// };

#endif // GUDHI_MULTI_FIELD_H 
