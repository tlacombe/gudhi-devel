import gudhi

"""This file is part of the Gudhi Library. The Gudhi library
   (Geometric Understanding in Higher Dimensions) is a generic C++
   library for computational topology.

   Author(s):       Pawel Dlotko

   Copyright (C) 2018 Swansea University.

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

__author__ = "Pawel Dlotko"
__copyright__ = "Copyright (C) 2018 Swansea University"
__license__ = "GPL v3"


def test_persistence_diagrams_range():
	diagram = [ (1,3),(2,4),(3,5) ]
	diag1 = gudhi.PersistenceIntervals(diagram)
	assert ( diag1.get_x_range() == (1.0, 5.0) )
	assert ( diag1.get_y_range() == (3.0, 5.0) )
	
	
	
def test_persistence_diagrams_characteristic_functions():
	diagram = [ (1,3),(2,4),(3,5) ]
	diag1 = gudhi.PersistenceIntervals(diagram)	
	assert ( diag1.characteristic_function_of_diagram(1,5,20) == [0.4, 0.4, 0.4, 0.4, 0.4, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.4, 0.4, 0.4, 0.4, 0.4] )
	assert (diag1.cumulative_characteristic_function_of_diagram(1,5,20) == [0.4, 0.8, 1.2000000000000002, 1.6, 2.0, 2.8, 3.5999999999999996, 4.3999999999999995, 5.199999999999999, 5.999999999999999, 6.799999999999999, 7.599999999999999, 8.399999999999999, 9.2, 10.0, 10.4, 10.8, 11.200000000000001, 11.600000000000001, 12.000000000000002])
	diagram2 = [ (1,2),(1,3),(1,4),(1,5) ]
	diag2 = gudhi.PersistenceIntervals(diagram2)	
	assert (diag2.histogram_of_lengths(11) == [0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1])
	assert (diag2.cumulative_histogram_of_lengths(11) == [0, 0, 1, 1, 1, 2, 2, 2, 3, 3, 3, 4] )	
	

def test_persistence_diagrams_persistence_betti_numners():
	diagram = [ (1,3),(2,4),(3,5) ]
	diag1 = gudhi.PersistenceIntervals(diagram)	
	assert ( diag1.compute_persistent_betti_numbers() == [(1.0, 1), (2.0, 2), (3.0, 1), (3.0, 2), (4.0, 1), (5.0, 0)] )	

def test_domnant_intervals():	
	diagram2 = [ (1,2),(1,3),(1,4),(1,5) ]
	diag2 = gudhi.PersistenceIntervals(diagram2)
	assert( diag2.dominant_intervals(1) == [(1.0, 5.0)] )  
	assert( diag2.dominant_intervals(2) == [(1.0, 5.0), (1.0, 4.0)] )  
	assert( diag2.dominant_intervals(3) == [(1.0, 5.0), (1.0, 4.0), (1.0, 3.0)] )  
	assert( diag2.dominant_intervals(4) == [(1.0, 5.0), (1.0, 4.0), (1.0, 3.0), (1.0, 2.0)] )  
	assert( diag2.dominant_intervals(5) == [(1.0, 5.0), (1.0, 4.0), (1.0, 3.0), (1.0, 2.0)] )  
	assert ( diag2.length_of_dominant_intervals(7) == [4.0, 3.0, 2.0, 1.0] )
	
	
