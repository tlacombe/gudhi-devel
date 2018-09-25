 /*    This file is part of the Gudhi Library. The Gudhi library 
  *    (Geometric Understanding in Higher Dimensions) is a generic C++ 
  *    library for computational topology.
  *
  *    Author(s):       Clément Maria
  *
  *    Copyright (C) 2018 Inria
  *
  *    This program is free software: you can redistribute it and/or modify
  *    it under the terms of the GNU General Public License as published by
  *    the Free Software Foundation, either version 3 of the License, or
  *    (at your option) any later version.
  *
  *    This program is distributed in the hope that it will be useful,
  *    but WITHOUT ANY WARRANTY; without even the implied warranty of
  *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  *    GNU General Public License for more details.
  *
  *    You should have received a copy of the GNU General Public License
  *    along with this program.  If not, see <http://www.gnu.org/licenses/>.
  */

/** \brief Iterator over a zigzag filtration (insertion and removal of cells). 
  * The dynamic nature of computations requires the iterator to encode methods, 
  * that are usually encoded by the data structure storing the cell complex.
  */
struct ZigzagfiltrationSimplexIterator {
/** Handle to specify a cell. */
  typedef unspecified      Simplex_handle;  
/** Star operator of the iterator, of value_type Simplex_handle.
  */
  Simplex_handle & operator*();
/** Comparison of iterators.
  */
  bool operator!=();
/** \brief Type for the value of the filtration function.
  *
  * Must be comparable with <. */
  typedef unspecified      Filtration_value;
/** \brief Returns the filtration value of the cell pointed at.
  */
  Filtration_value filtration();
/** \brief Returns an upper bound on the dimension of the cells iterated over.
  *
  * In particular used in in order to not record maximal dimension intervals 
  * in the zigzag persistence barcode.
  */
  int dim_max();
/** \brief Returns the direction (insertion or deletion) of the corresponding map.
  *
  * true for insertion, false for deletion.
  */
  bool arrow_direction();
};


/** \brief The concept ZigzagFilteredComplex describes the requirements 
  * for a type to implement a filtered cell complex, from which 
  * one can compute zigzag persistent homology via the class 
  * Zigzag_persistence< ZigzagFilteredComplex >
  */
struct ZigzagFilteredComplex
{
/** \brief Data stored for each simplex. 
  *
  * Must be an integer type. */
  typedef unspecified      Simplex_key;
/** Handle to specify a simplex. */
  typedef unspecified      Simplex_handle;
/** \brief Type for the value of the filtration function.
  *
  * Must be comparable with <. */
  typedef unspecified      Filtration_value;
/** \brief Returns the number stored for a simplex by `assign_key`.
  *
  * This is never called on null_simplex(). */
  Simplex_key              key      ( Simplex_handle sh );
/** \brief Store a number for a simplex, which can later be retrieved with 
  * `key(sh)`.
  *
  * This is never called on null_simplex(). */
  void                     assign_key(Simplex_handle sh, Simplex_key n);
/** \brief Iterator over all simplices of the complex 
  * in the order of the indexing scheme. Must be a model of
  * ZigzagfiltrationSimplexIterator
  *
  * 'value_type' must be 'Simplex_handle'.
  */
typedef unspecified Zigzagfiltration_simplex_iterator;
/** \brief Range over the simplices of the complex
  * in the order of the filtration.
  *
  * .begin() and .end() return type Zigzagfiltration_simplex_iterator.*/
typedef unspecified Zigzagfiltration_simplex_range;
/** \brief Returns a range over the simplices of the complex
  * in the order of the filtration.
  *
  * .begin() and .end() return type Filtration_simplex_iterator.*/
Zigzagfiltration_simplex_range filtration_simplex_range();
/** \brief Returns the dimension of a simplex. */
  int dimension(Simplex_handle sh);
/** \brief In a Morse theoretical framework, returns true if the simplex is critical
  * with regards to a Morse matching, false if it is paired.
  */
  bool critical(Simplex_handle sh);
/** \brief If the simplex is paired in a Morse matching, returns a Simplex_hanbdle 
  * to the simplex it is paired with. If sh is not paired (i.e., critical), 
  * morse_pair returns sh itself.
  */
  Simplex_handle morse_pair(Simplex_handle sh);
/** \brief Iterator over the simplices belonging to the
  * boundary of a simplex.
  *
  * <CODE>value_type</CODE> must be 'Simplex_handle'.
  */
typedef unspecified Boundary_simplex_iterator;
/** \brief Range giving access to the simplices in the boundary of 
  * a simplex.
  *
  * .begin() and .end() return type Boundary_simplex_iterator.
  */
typedef unspecified Boundary_simplex_range;

/** \brief Returns a range giving access to all simplices of the 
  * boundary of a simplex, i.e.
  * the set of codimension 1 subsimplices of the Simplex.
  *
  * If the simplex is \f$[v_0, \cdots ,v_d]\f$, with canonical orientation
  * induced by \f$ v_0 < \cdots < v_d \f$, the iterator enumerates the 
  * simplices of the boundary in the order: 
  * \f$[v_0,\cdots,\widehat{v_i},\cdots,v_d]\f$ for \f$i\f$ from 0 to d
  *
  * We note that the alternate sum of the simplices given by the iterator
  * gives the chains corresponding to the boundary of the simplex.*/
Boundary_simplex_range boundary_simplex_range(Simplex_handle sh);
/** \brief Iterator over the simplices belonging to the
  * coboundary of a simplex (cofaces of codimension 1).
  *
  * <CODE>value_type</CODE> must be 'Simplex_handle'.
  */
typedef unspecified Coboundary_simplex_iterator;
/** \brief Range giving access to the simplices in the coboundary of 
  * a simplex.
  *
  * .begin() and .end() return type Coboundary_simplex_iterator.
  */
typedef unspecified Coboundary_simplex_range;

/** \brief Returns a range giving access to all simplices of the 
  * coboundary of a simplex, i.e.
  * the set of codimension 1 cofaces of the Simplex.
  */
Coboundary_simplex_range coboundary_simplex_range(Simplex_handle sh);
/** \brief Specifies the nature of the indexing scheme. 
  * 
  * is model of IndexingTag. */
  typedef unspecified      Indexing_tag;

/** Returns a Simplex_handle that is different from all simplex handles 
  * of the simplices. */
  Simplex_handle           null_simplex();
/** \brief Returns a constant dummy number that is either negative,
  * or at least as large as `num_simplices()`.  Suggested value: -1.  */
  Simplex_key              null_key ();
};
