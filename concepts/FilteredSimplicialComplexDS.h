/*
 *  FilteredSimplicialComplexDS.h
 *  GUDHI
 *
 *  Created by Cl??ment Maria on 12/10/13.
 *  Copyright 2013 __MyCompanyName__. All rights reserved.
 *
 */

/**
 *  \brief Data structure for representing a filtered simplicial complex.
 *  \extends SimplicialComplexDS
 */

struct FilteredSimplicialComplexDS
{
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

  /**
   * Iterator over all simplices of a complex in the order of the filtration
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
};
