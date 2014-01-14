/*
 *  FilteredSimplicialComplexDS.h
 *  GUDHI
 *
 *  Created by Clément Maria on 12/10/13.
 *  Copyright 2013 INRIA. All rights reserved.
 *
 */

/**
 *  \brief Data structure for representing a filtered simplicial complex.
 *
 *  It is designed to be used by any algorithm to compute persistent homology
 *  and persistent cohomology as long as no additional data needs to be stored 
 *  in the simplices by the persistence algorithm.
 *
 *  \extends SimplicialComplexDS
 */

struct FilteredSimplicialComplexDS
{
	/*************************************************/		
	/// \name Boundary Iterator
	/// @{
	/**
	 *	 Iterator over all simplices of the boundary of a simplex, i.e.
	 *   the set of codimension $1$ subsimplices of the Simplex.
	 *
	 *	 `value_type` must be a `Simplex`.
	 */
	typedef unspecified Boundary_simplex_iterator;
	/**
	 * Returns an iterator to the beginning of the sequence of
	 * simplices of the boundary of a simplex
	 *
	 * If the simplex is \f$[v_0, \cdots ,v_d]\f$, with canonical orientation
	 * induced by \f$ v_0 < \cdots < v_d \f$, the iterator enumerates the 
	 * simplices of the boundary in the order: 
	 * \f$[v_0,\cdots,\widehat{v_i},\cdots,v_d]\f$ for \f$i\f$ from 0 to d
	 *
	 * We note that the alternate sum of the simplices given by the iterator
	 * gives the chains corresponding to the boundary of the simplex.
	 */
	Boundary_simplex_iterator boundary_simplex_begin(Simplex s);
	/**
	 * Returns an iterator to the end of the sequence of 
	 * simplices of the boundary of a simplex
	 */
	Boundary_simplex_iterator boundary_simplex_end(Simplex s);
	/// @}
	/*************************************************/		
	
	/*************************************************/		
	/// \name Filtration Iterator
	/// @{
	/**
	 * Iterator over all simplices of the complex in the order of the filtration
	 *
	 * `value_type` must be a `Simplex`
	 */
	typedef unspecified Filtration_simplex_iterator;
	/**
	 *  Returns an iterator to the beginning of the sequence of
	 *  simplices of a complex in the order of the filtration
	 */
	Filtration_simplex_iterator filtration_simplex_begin();
	/**
	 *  Returns an iterator to the end of the sequence of 
	 *  simplices of a complex in the order of the filtration
	 */
	Filtration_simplex_iterator filtration_simplex_end();
	/// @}
	/*************************************************/		
	
	
	
	/*************************************************/		
	/**
	 * @details Type of the filtration values
	 *
	 * @note OPTIONAL
	 */
	typedef unspecified Filtration_value;
	/**
	 * @details Returns the filtration value of a simplex s
	 *
	 * @note OPTIONAL
	 */
	Filtration_value filtration_value(Simplex s);
	/**
	 * @details Compares the filtration values of simplices s and t
	 *
	 * @return -1 if s comes before t in the filtration, +1
	 * if it comes later, and 0 if they come at the same time
	 *
	 * @note OPTIONAL
	 * @todo use an enum? Just a bool?
	 */
	int compare_in_filtration(Simplex s, Simplex t);
	/*************************************************/		
	
};
