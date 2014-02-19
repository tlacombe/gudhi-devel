/*
*  FilteredSimplicialComplexDS.h
*  Gudhi
*
*  Created by Clément Maria on 12/10/13.
*  Copyright 2013 INRIA. All rights reserved.
*
*/

/** \brief Data structure for representing a filtered simplicial complex.
  *
  *  It is designed to be used by any algorithm to compute persistent homology
  *  and persistent cohomology as long as no additional data needs to be stored 
  *  in the simplices by the persistence algorithm.
  *
  *  \extends SimplicialComplexDS
  */
struct FilteredSimplicialComplexDS
{
/** \brief Iterator on the simplices belonging to the
  * boundary of a simplex.
  *
  * 'value_type' must be 'Simplex_handle'.
  */
typedef unspecified Boundary_simplex_iterator;
/** \brief Range over the simplices in the boundary of 
  * a simplex.
  *
  * .begin() and .end() return type Boundary_simplex_iterator.
  */
typedef unspecified Boundary_simplex_range;

/** \brief Range over all simplices of the boundary of a simplex, i.e.
  *  the set of codimension $1$ subsimplices of the Simplex.
  *
  * If the simplex is \f$[v_0, \cdots ,v_d]\f$, with canonical orientation
  * induced by \f$ v_0 < \cdots < v_d \f$, the iterator enumerates the 
  * simplices of the boundary in the order: 
  * \f$[v_0,\cdots,\widehat{v_i},\cdots,v_d]\f$ for \f$i\f$ from 0 to d
  *
  * We note that the alternate sum of the simplices given by the iterator
  * gives the chains corresponding to the boundary of the simplex.*/
Boundary_simplex_range boundary_simplex_range(Simplex_handle sh);

/** \brief Iterator over all simplices of the complex 
  * in the order of the filtration.
  *
  * 'value_type' must be 'Simplex_handle'.
  */
typedef unspecified Filtration_simplex_iterator;
/** \brief Range over the simplices of the complex
  * in the order of the filtration.
  *
  * .begin() and .end() return type Filtration_simplex_iterator.*/
typedef unspecified Filtration_simplex_range;
/** \brief Returns a range over the simplices of the complex
  * in the order of the filtration.
  *
  * .begin() and .end() return type Filtration_simplex_iterator.*/
Filtration_simplex_range filtration_simplex_range();


/** \brief Iterator over the simplices of the complex,
  * in an arbitrary order.
  *
  * 'value_type' must be 'Simplex_handle'.*/
//typedef unspecified Complex_simplex_iterator;
//typedef unspecified Complex_simplex_range;

/**
* Returns a range over all the simplices of a
* complex.
*/
//Complex_simplex_range complex_simplex_range();

/*************************************************/		
/**
* @details Type of the filtration values
*
* @note OPTIONAL
*/
//typedef unspecified Filtration_value;
/**
* @details Returns the filtration value of a simplex s
*
* @note OPTIONAL
*/
//Filtration_value filtration_value(Simplex_handle s);
/**
* @details Compares the filtration values of simplices s and t
*
* @return -1 if s comes before t in the filtration, +1
* if it comes later, and 0 if they come at the same time
*
* @note OPTIONAL
* @todo use an enum? Just a bool?
*/
//int is_before_in_filtration(Simplex_handle s, Simplex_handle t);
/*************************************************/		

};
