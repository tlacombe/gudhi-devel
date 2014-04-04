/*
 *  Field_Zp.h
 *  Gudhi
 *
 *  Created by Clément Maria on 02/28/14.
 *  Copyright 2013 INRIA. All rights reserved.
 *
 */

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
  if(res < 0) return res+Prime;
  else return res; 
}

void plus_equal(Element & x, Element y) { x = ((x+y)%Prime); }

/** \brief Returns the additive idendity \f$0_{\F}\f$ of the field.*/
Element additive_identity () { return 0; }
/** \brief Returns the multiplicative identity \f$1_{\F}\f$ of the field.*/
Element multiplicative_identity ( Element P = Prime ) { return 1; }
/** Returns the inverse in the field. Modifies P.*/
Element inverse ( Element x
                , Element & P ) 
{ P = 1; return inverse_[ x ]; 
}  // <------ return the product of field characteristic for which x is invertible




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
