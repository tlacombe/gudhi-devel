import gudhi	
import math
import sys

"""
    This file is part of the Gudhi Library. The Gudhi library
    (Geometric Understanding in Higher Dimensions) is a generic C++
    library for computational topology.

    Author(s):       Pawel Dlotko

    Copyright (C) 2017  Swansea University

    This program is free software: you can redistribute it and/or modify
     it under the terms of the GNU General Public License as published by
     the Free Software Foundation, either version 3 of the License, or
     (at your option) any later version.

    This program is distributed in the hope that it will be useful,
     but WITHOUT ANY WARRANTY without even the implied warranty of
     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
     GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
     along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

__author__ = "Pawel Dlotko"
__copyright__ = "Copyright (C) 2018 Swansea University"
__license__ = "GPL v3"

epsilon = 0.0000005



def test_read_write_to_files():
	intervals = [ (2, 3),(4, 7),(9, 10),(3, 11) ]
	p = gudhi.PersistenceVectors(vector_of_intervals=intervals)
	p.print_to_file("test_vector_representation_write_read")

	q = gudhi.PersistenceVectors()
	q.load_from_file("test_vector_representation_write_read")

	assert p.compare(q)


def test_sortev_vector_distances_template(): 
	p = gudhi.PersistenceVectors(file_with_intervals="data/file_with_diagram", where_to_cut=100)
	sortev_vector_distances_template = [  0.609968,  0.566317,  0.538858,  0.534927,  0.515741,  0.507828,  0.500911,  0.496986,  0.495306,  0.439945,  0.424097,  0.413891,  0.413891,  0.413891,  0.412621,  0.410613,  0.407853,  0.407853,  0.402306,  0.401937,  0.377605,  0.363859,  0.357765,  0.357765,  0.357765,  0.357765,  0.357765,  0.353401,  0.348004,  0.345124,  0.345124,  0.345124,  0.345124,  0.345124,  0.345124,  0.345124,  0.345124,  0.345124,  0.345124,  0.345124,  0.345124,  0.345124,  0.345124,  0.34469,  0.339466,  0.33935,  0.32834,  0.327276,  0.318626,  0.318082,  0.30603,  0.30525,  0.297308,  0.296333,  0.296333,  0.296333,  0.296333,  0.293372,  0.292666,  0.292666,  0.292666,  0.292666,  0.292666,  0.292666,  0.292666,  0.292666,  0.292666,  0.292666,  0.292666,  0.292666,  0.292666,  0.292666,  0.292666,  0.29029,  0.290218,  0.289782,  0.288128,  0.286416,  0.285969,  0.282046,  0.28154,  0.281085,  0.280227,  0.279273,  0.278936,  0.278706,  0.278507,  0.278097,  0.276293,  0.276293,  0.276293,  0.276293,  0.276293,  0.276293,  0.276293,  0.276293,  0.276293,  0.276169,  0.270563,  0.264009]
	proj_no = p.number_of_vectorize_functions()
	p_vect = p.vectorize(proj_no)

	for i in range(0, len(p_vect)):
		assert math.fabs(sortev_vector_distances_template[i] - p_vect[i]) < epsilon
  


def test_projections_to_R():
	p = gudhi.PersistenceVectors(file_with_intervals="data/file_with_diagram", where_to_cut=100)
	proj = [0,0.6099679993,1.176284775,1.715142954,2.25006986,2.765810506,3.273638431,3.774549309,4.271535042,4.766840875,5.206786149,5.63088295,6.04477433,6.45866571,6.87255709,7.285177939,7.695791381,8.103643945,8.511496508,8.913802775,9.315740229,9.693344927,10.0572035,10.41496899,10.77273448,11.13049996,11.48826545,11.84603094,12.19943233,12.5474364,12.89256042,13.23768444,13.58280846,13.92793248,14.2730565,14.61818051,14.96330453,15.30842855,15.65355257,15.99867659,16.34380061,16.68892462,17.03404864,17.37917266,17.7238622,18.06332781,18.40267754,18.73101759,19.05829313,19.3769189,19.69500045,20.0010306,20.30628026,20.60358868,20.89992192,21.19625516,21.4925884,21.78892164,22.08229394,22.37495987,22.66762581,22.96029174,23.25295768,23.54562361,23.83828955,24.13095549,24.42362142,24.71628736,25.00895329,25.30161923,25.59428516,25.8869511,26.17961703,26.47228297,26.76257262,27.05279049,27.34257265,27.63070097,27.91711687,28.20308566,28.48513176,28.76667161,29.04775635,29.32798359,29.60725702,29.88619335,30.16489915,30.44340655,30.72150329,30.99779604,31.27408878,31.55038153,31.82667427,32.10296702,32.37925976,32.6555525,32.93184525,33.20813799,33.48430662,33.7548692]
	for proj_no in range(0,len(proj)):		
		assert math.fabs(p.project_to_R(proj_no) - proj[proj_no]) < epsilon
  


def test_distance_computations():
	p = gudhi.PersistenceVectors(file_with_intervals="data/file_with_diagram", where_to_cut=100)  
	intervals = [(1, 2) , (3, 4), (5, 6), (7, 8), (9, 10), (11, 12), (13, 14), (15, 16), (17, 18), (19, 20)]
	p_bis = gudhi.PersistenceVectors(vector_of_intervals=intervals, where_to_cut=10)  
	assert math.fabs(p.distance(p_bis, 1)- 30.676373638271652) < epsilon


def test_default_parameters_of_distances():
	diag = gudhi.read_persistence_intervals_in_dimension("data/file_with_diagram")
	p = gudhi.PersistenceVectors(vector_of_intervals=diag, where_to_cut=100)

	diag1 = gudhi.read_persistence_intervals_in_dimension("data/file_with_diagram_1")
	q = gudhi.PersistenceVectors(vector_of_intervals=diag1, where_to_cut=100)

	dist_numeric_limit_max = p.distance(q, sys.maxint)
	dist_infinity = p.distance(q, sys.maxint)

	assert dist_numeric_limit_max == dist_infinity


def test_compute_average():
	i1 = [ (1, 2),(3, 8),(1, 6) ]
	i2 = [(2, 9), (2, 15),(6, 17)]

	A = gudhi.PersistenceVectors(vector_of_intervals=i1)
	B = gudhi.PersistenceVectors(vector_of_intervals=i1)

	average = gudhi.PersistenceVectors()
	average.compute_average([A, B])

	template_average = gudhi.PersistenceVectors()
	template_average.load_from_file("data/average_of_persistence_vectors")

	assert template_average.compare(average)
