/*
 *  CAMFilteredSimplicialComplexDS.h
 *  GUDHI
 *
 *  Created by Cl??ment Maria on 12/10/13.
 *  Copyright 2013 __MyCompanyName__. All rights reserved.
 *
 */

/**
 * \brief Data structure for representing a filtered simplicial complex that
 * can be used with the Compressed Annotation Matrix algorithm.
 * \extends FilteredSimplicialComplexDS
 * @todo Replace Simplex_data with exactly what CAM needs.
 */
struct CAMFilteredSimplicialComplexDS
{
	/**
	 * Type of the extra space provided for each simplex
	 * to store additional information
	 */
	typedef unspecified Simplex_data;
	
	/**
	 * Returns a pointer to the extra information provided for
	 * the simplex s.  If not `Simplex_data` is not `void`, it must be
	 * persistent, i.e. it cannot be part of a temporary Simplex object.
	 */
	Simplex_data *simplex_data_pointer(Simplex s);
};
