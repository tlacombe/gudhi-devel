/*    This file is part of the Gudhi Library. The Gudhi library
*    (Geometric Understanding in Higher Dimensions) is a generic C++
*    library for computational topology.
*
*    Author(s):       Clément Jamin
*
*    Copyright (C) 2017  INRIA
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

#ifndef DOC_COMMON_FILE_FORMAT_H_
#define DOC_COMMON_FILE_FORMAT_H_

namespace Gudhi {

/*! \page fileformats File formats

 \tableofcontents

 \section FileFormatsPers Persistence Diagram

 Such a file, whose extension is usually `.pers`, contains a list of persistence intervals.<br>
 Lines starting with `#` are ignored (comments).<br>
 Other lines might contain 2, 3 or 4 values (the number of values on each line must be the same for all lines):
 \verbatim
   [[field] dimension] birth death
 \endverbatim

 Here is a simple sample file:
 \verbatim
   # Persistence diagram example
   2 2.7 3.7
   2 9.6 14.
   # Some comments
   3 34.2 34.974
   4 3. inf
 \endverbatim

 Other sample files can be found in the `data/persistence_diagram` folder.

 Such files can be generated with `Gudhi::persistent_cohomology::Persistent_cohomology::output_diagram()` and read with
 `Gudhi::read_persistence_intervals_and_dimension()`, `Gudhi::read_persistence_intervals_grouped_by_dimension()` or
 `Gudhi::read_persistence_intervals_in_dimension()`.
  
  
 \section FileFormatHasseDiagram Format of file based on which Hasse diagram can be created. 
  Such a Hasse diagram can be read by Hasse_diagram/utilities/Hasse_diagram_from_file.cpp
  Lines starting with `#` are ignored (comments).
  We assume that the file do not contain empty lines.
  
  The cells stored in that file are assumed to be enumerated with integers starting from zero.
  
  The first line contains a non negative integer N determining a number of cells in the Hasse diagram.  
  Subsequently the file contains N blocks. Each block represent a cell in the chain complex. Below a description 
  of a block is given.
   
  First line of a block consist of two or three numbers: number of a cell (between zero and N-1), dimension of cell (nonnegative integer). 
  The third (optional) number is a filtration of a cell.
  The next line contains sequence of ids of boundary elements of a given cell alternated by the incidence coefficient between the given cell and the boundary element.
 
  For an exampe of a file, please consult Hasse_diagram/test/cw_decomposition_of_torus.hasse    
*/
}  // namespace Gudhi

#endif  // DOC_COMMON_FILE_FORMAT_H_
