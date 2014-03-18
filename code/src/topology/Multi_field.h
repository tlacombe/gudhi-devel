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


//#include "boost/pool.hpp" 


/**
*
*/
template < int P , int Q >
class Multi_field {
  typedef uint64_t    Element;

  Multi_field ()
  { /* ... */ }

  /** Set x <- x + w * y*/
  void plus_equal(Element & x, Element y, int w);
  /** Returns 0.*/
  Element additive_identity() { return 0; };
  /** Returns the partial identity w.r.t. all fields.*/
  Element multiplicative_identity(); // as if prod == product of all primes
  /** Returns the partial identity w.r.t. prod.*/
  Element multiplicative_identity(Element prod);
  /** Computes the partial inverse of x w.r.t. prod. <- modifies prod */
  Element inverse ( Element x, Element prod ); // <-- modifies prod
    
  bool is_one ( Element x );
  bool is_zero ( Element x );
  /** Returns the product of the characteristic of the fields in the multi-field
  * structure.*/
  Element characteristic() { return prod_primes; }

private:
  Element prod_primes; //product of all characteristics

};




/** \brief */
template < int Prime = 11 >
class Field_Zp {
public:
typedef int Element;

Field_Zp() 
: inverse_(Prime)
{ 
  inverse_[0] = 0;
  for(int i=1 ; i<Prime ; ++i)
  { 
    int inv = 1; 
    while(((inv * i) % Prime) != 1) ++inv;
    inverse_[i] = inv; 
  }
}

/** Set x <- x + w * y*/
void plus_times_equal ( Element & x, Element y, Element w ) 
{ x = (x + w * y) % Prime; }

// operator= defined on Element


/** Returns y * w */
Element times ( Element y, int w ) { 
  Element res = (y * w) % Prime;
  if(res > 0) return res;
  else return res+Prime; 
}

void plus_equal(Element & x, Element y) { x = ((x+y)%Prime); }

/** \brief Returns the additive idendity \f$0_{\F}\f$ of the field.*/
Element additive_identity () { return 0; }
/** \brief Returns the multiplicative identity \f$1_{\F}\f$ of the field.*/
Element multiplicative_identity ( Element P = Prime ) { return 1; }
/** Returns the inverse in the field. Modifies P.*/
Element inverse ( Element x
                , Element & P ) 
{
 //std::cout << "   inverse: " << x << " " << inverse_[x] << std::endl;
 P = 1; return inverse_[ x ]; 

}
/** Returns -x * y.*/
Element times_minus ( Element x, Element y ) 
{ 
  Element out = (-x * y) % Prime;
  return (out < 0) ? out + Prime : out; 
}


bool is_one ( Element x ) { return x == 1; }
bool is_zero ( Element x ) { return x == 0; }

//bool is_null()

/** \brief Returns the characteristic \f$p\f$ of the field.*/
Element characteristic() { return Prime; }

private:
/** Property map Element -> Element, which associate to an element its inverse in the field.*/
std::vector< Element > inverse_;
//boost::object_pool< Column_cell > * cell_pool_;
};

/**
*
*/
class Field_Q {

};













#endif // GUDHI_MULTI_FIELD_H 








