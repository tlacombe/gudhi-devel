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
     but WITHOUT ANY WARRANTY; without even the implied warranty of
     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
     GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
     along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

__author__ = "Pawel Dlotko"
__copyright__ = "Copyright (C) 2018 Swansea University"
__license__ = "GPL v3"

epsilon = 0.0000005;

#I do not understand how come those tests are passing, since the module should be 
#called gudhi.PersistenceLandscapes....

def test_check_construction_of_landscape():
	p = gudhi.PersistenceLandscapes(file_with_intervals="data/file_with_diagram")
	q = gudhi.PersistenceLandscapes()
	q.load_landscape_from_file("data/file_with_landscape_from_file_with_diagram")
	assert p.compare(q)   


def test_check_construction_of_landscape_form_gudhi_style_file():
	p = gudhi.PersistenceLandscapes(file_with_intervals="data/persistence_file_with_four_entries_per_line", dimension=1)
	q = gudhi.PersistenceLandscapes()
	q.load_landscape_from_file("data/persistence_file_with_four_entries_per_line_landscape");  
	assert p.compare(q)

def test_check_computations_of_integrals():  
	p = gudhi.PersistenceLandscapes(file_with_intervals="data/file_with_diagram")
	integral = p.compute_integral_of_landscape()
	assert math.fabs(integral - 2.34992) <= 0.00001
    
    
def test_check_construction_of_landscape():
	diag = gudhi.read_persistence_intervals_in_dimension(persistence_file='data/file_with_diagram')
	p = gudhi.PersistenceLandscapes(diag)
	q = gudhi.PersistenceLandscapes(file_with_intervals="data/file_with_diagram")    
	assert p.compare(q) 
	q = gudhi.PersistenceLandscapes()
	q.load_landscape_from_file("data/file_with_landscape_from_file_with_diagram")
	assert p.compare(q) 
	diag = [(0.1497983027, 0.9062991885), (0.7176230095, 0.730709507), (0.0544614334, 0.3223612854), (0.2674207884, 0.8895973614),
	(0.1234925873, 0.6294493901), (0.1769560101, 0.1832335867), (0.2473518713, 0.9767192878), (0.8728962245, 0.9756544009),
	(0.5930361075, 0.7011736326), (0.1045988346, 0.8666593796), (0.0920032915, 0.2323297802), (0.1554214768, 0.2199019471),
	(0.1910692175, 0.7908326657), (0.1918362495, 0.9002312066), (0.1926749829, 0.9108521216), (0.650780786, 0.7222979715),
	(0.3723952303, 0.649569716), (0.1173791633, 0.1377982685),  (0.1147177995, 0.9773428167), (0.1584430523, 0.4892093645),
	(0.3877435371, 0.4629111887), (0.2345582605, 0.7220228854), (0.2847216921, 0.8536686709), (0.2184254499, 0.5872268199),
	(0.2085643625, 0.2787338524), (0.2680930813, 0.6871715707), (0.0290362441, 0.5171153115), (0.8493046183, 0.9447571123),
	(0.7709131132, 0.9901270741), (0.3710208764, 0.9478115036), (0.4355365708, 0.5326650707), (0.0307675984, 0.2371657227),
	(0.2940688241, 0.9945370505), (0.545262265, 0.9591544003), (0.6834892731, 0.7583554041), (0.483568935, 0.5221970463),
	(0.3231522712, 0.9038473698), (0.6549508263, 0.9801285535), (0.1336375578, 0.9345304225), (0.605044048, 0.636718961),
	(0.2849982539, 0.9878422172), (0.1913332292, 0.7253469135), (0.6026467816, 0.8212106009), (0.0366311713, 0.6219619742),
	(0.3062930094, 0.6970299556)];
	p = gudhi.PersistenceLandscapes(diag)
	assert p.compare(q) 
    
    


def test_check_computations_of_integrals_for_each_level_separatelly():
	diag = gudhi.read_persistence_intervals_in_dimension(persistence_file='data/file_with_diagram')
	p = gudhi.PersistenceLandscapes(diag)    
	integrals_for_different_levels = [0.216432,0.204763,0.188793,0.178856,0.163142,0.155015,0.143046,0.133765,0.123531,0.117393,0.111269,0.104283,0.0941308,0.0811208,0.0679001,0.0580801,0.0489647,0.0407936,0.0342599,0.02896,0.0239881,0.0171792,0.0071511,0.00462067,0.00229033,0.000195296]
	for lv in range(0, len(integrals_for_different_levels)):
		integral = p.compute_integral_of_a_level_of_a_landscape(lv);        
		assert math.fabs(integral - integrals_for_different_levels[lv]) <= 0.00001
        
        

def test_check_computations_of_integrals_of_powers_of_landscape():
	diag = gudhi.read_persistence_intervals_in_dimension("data/file_with_diagram")
	p = gudhi.PersistenceLandscapes(diag)
	integrals_fir_different_powers = [17.1692,2.34992,0.49857,0.126405,0.0355235]
	for power in range(0,5):
		integral = p.compute_integral_of_power_of_landscape(power)
		assert math.fabs(integral - integrals_fir_different_powers[power]) <= 0.00005
  

def test_check_computations_of_values_on_different_points():
	diag = gudhi.read_persistence_intervals_in_dimension("data/file_with_diagram")
	p = gudhi.PersistenceLandscapes(diag);
	assert math.fabs(p.compute_value_at_a_given_point(1, 0.0)) <= 0.00001 
	assert math.fabs(p.compute_value_at_a_given_point(1, 0.1) - 0.0692324) <= 0.00001
	assert math.fabs(p.compute_value_at_a_given_point(1, 0.2) - 0.163369) <= 0.00001
	assert math.fabs(p.compute_value_at_a_given_point(1, 0.3) - 0.217115) <= 0.00001 
	assert math.fabs(p.compute_value_at_a_given_point(2, 0.0)) <= 0.00001
	assert math.fabs(p.compute_value_at_a_given_point(2, 0.1) - 0.0633688) <= 0.00001
	assert math.fabs(p.compute_value_at_a_given_point(2, 0.2) - 0.122361) <= 0.00001
	assert math.fabs(p.compute_value_at_a_given_point(2, 0.3) - 0.195401) <= 0.00001
	assert math.fabs(p.compute_value_at_a_given_point(3, 0.0)) <= 0.00001
	assert math.fabs(p.compute_value_at_a_given_point(3, 0.1) - 0.0455386) <= 0.00001
	assert math.fabs(p.compute_value_at_a_given_point(3, 0.2) - 0.0954012) <= 0.00001
	assert math.fabs(p.compute_value_at_a_given_point(3, 0.3) - 0.185282) <= 0.00001


def test_check_computations_of_maxima_and_norms():
	diag = gudhi.read_persistence_intervals_in_dimension("data/file_with_diagram")
	p = gudhi.PersistenceLandscapes(diag)
	second = gudhi.PersistenceLandscapes()
	second.load_landscape_from_file("data/file_with_landscape_from_file_with_diagram_1")
	assert math.fabs(p.compute_maximum() - 0.431313) <= 0.00001
	assert math.fabs(p.compute_norm_of_landscape(1) - 2.34992) <= 0.00001
	assert math.fabs(p.compute_norm_of_landscape(2) - 0.706095) <= 0.00001
	assert math.fabs(p.compute_norm_of_landscape(3) - 0.501867) <= 0.00001


def test_check_default_parameters_of_distances():
	diag = gudhi.read_persistence_intervals_in_dimension("data/file_with_diagram")
	p = gudhi.PersistenceLandscapes(diag)    
	q = gudhi.PersistenceLandscapes()
	q.load_landscape_from_file("data/file_with_landscape_from_file_with_diagram_1")
	dist_numeric_limit_max = p.distance(q, sys.float_info.max);
	dist_infinity = p.distance(q, sys.float_info.max);
	assert dist_numeric_limit_max == dist_infinity


def test_check_computations_of_averages():
	diag = gudhi.read_persistence_intervals_in_dimension("data/file_with_diagram")
	p = gudhi.PersistenceLandscapes(diag)
	q = gudhi.PersistenceLandscapes()
	q.load_landscape_from_file("data/file_with_landscape_from_file_with_diagram_1")
	av = gudhi.PersistenceLandscapes()
	av.compute_average([p, q])
	template_average = gudhi.PersistenceLandscapes()
	template_average.load_landscape_from_file("data/average")
	assert template_average.compare(av)


def test_check_computations_of_distances():
	diag = gudhi.read_persistence_intervals_in_dimension("data/file_with_diagram")
	p = gudhi.PersistenceLandscapes(diag)
	q = gudhi.PersistenceLandscapes()
	q.load_landscape_from_file("data/file_with_landscape_from_file_with_diagram_1")
	assert math.fabs(p.distance(q,1) - 25.5824) <= 0.00005
	assert math.fabs(p.distance(q, 2) - 2.1262155641322) <= 0.00001
	assert math.fabs(p.distance(q, sys.float_info.max) - 0.359068) <= 0.00001


def test_check_computations_of_scalar_product():
	diag = gudhi.read_persistence_intervals_in_dimension("data/file_with_diagram")
	p = gudhi.PersistenceLandscapes(diag)
	q = gudhi.PersistenceLandscapes()
	q.load_landscape_from_file("data/file_with_landscape_from_file_with_diagram_1")
	assert math.fabs(p.compute_scalar_product(q) - 0.754498) <= 0.00001

