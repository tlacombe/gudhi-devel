import gudhi	
import math


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

epsilon = 0.0000005



def test_check_construction_of_heat_map():    
	p = gudhi.PersistenceHeatMaps(file_with_intervals="data/file_with_diagram",how_many_pixels_raidus_of_Gausian_kernel=100,min_=0,max_=1)
	q = gudhi.PersistenceHeatMaps()
	q.load_from_file("data/persistence_heat_map_from_file_with_diagram")
	assert p.compare(q)
    
    
    
def test_construction_of_heat_maps():  
	p = gudhi.PersistenceHeatMaps(file_with_intervals="data/file_with_diagram", how_many_pixels_raidus_of_Gausian_kernel=100, min_=0, max_=1)
	p.print_to_file("data/persistence_heat_map_from_file_with_diagram")
	q = gudhi.PersistenceHeatMaps()
	q.load_from_file("data/persistence_heat_map_from_file_with_diagram")
	assert p.compare(q)


def test_averages_of_heat_maps():
	p = gudhi.PersistenceHeatMaps(file_with_intervals="data/file_with_diagram", how_many_pixels_raidus_of_Gausian_kernel=30, min_=0, max_=10)
	q = gudhi.PersistenceHeatMaps(file_with_intervals="data/file_with_diagram_1", how_many_pixels_raidus_of_Gausian_kernel=30, min_=0, max_=10)
	r = gudhi.PersistenceHeatMaps(file_with_intervals="data/file_with_diagram_2", how_many_pixels_raidus_of_Gausian_kernel=30, min_=0, max_=10)

	av = gudhi.PersistenceHeatMaps()
	av.compute_average([p,q,r])

	template_average = gudhi.PersistenceHeatMaps()
	template_average.load_from_file("data/template_average_of_heat_maps")

	assert av.compare(template_average)


def test_median_of_heat_maps():  
	p = gudhi.PersistenceHeatMaps(file_with_intervals="data/file_with_diagram",how_many_pixels_raidus_of_Gausian_kernel=30, min_=0, max_=1)
	q = gudhi.PersistenceHeatMaps(file_with_intervals="data/file_with_diagram_1", how_many_pixels_raidus_of_Gausian_kernel=30, min_=0, max_=1)
	r = gudhi.PersistenceHeatMaps(file_with_intervals="data/file_with_diagram_2", how_many_pixels_raidus_of_Gausian_kernel=30, min_=0, max_=1)

	to_compute_median = [p,q,r]
	median = gudhi.PersistenceHeatMaps()
	median.compute_median(to_compute_median)

	template_median = gudhi.PersistenceHeatMaps()
	template_median.load_from_file("data/template_median_of_heat_maps")

	assert median.compare(template_median)


def test_compute_percentage_of_active_of_heat_maps():  
	p = gudhi.PersistenceHeatMaps(file_with_intervals="data/file_with_diagram",how_many_pixels_raidus_of_Gausian_kernel=30, min_=0, max_=1)
	q = gudhi.PersistenceHeatMaps(file_with_intervals="data/file_with_diagram_1",how_many_pixels_raidus_of_Gausian_kernel=30, min_=0, max_=1)
	r = gudhi.PersistenceHeatMaps(file_with_intervals="data/file_with_diagram_2",how_many_pixels_raidus_of_Gausian_kernel=30, min_=0, max_=1)

	to_compute_percentage_of_active = [p,q,r]
	percentage_of_active = gudhi.PersistenceHeatMaps()
	percentage_of_active.compute_percentage_of_active(to_compute_percentage_of_active, 0.1)

	template_percentage_of_active = gudhi.PersistenceHeatMaps()
	template_percentage_of_active.load_from_file("data/template_percentage_of_active_of_heat_maps")

	assert percentage_of_active.compare(template_percentage_of_active)


def test_vectorize_for_heat_maps():  
	p = gudhi.PersistenceHeatMaps(file_with_intervals="data/file_with_diagram", how_many_pixels_raidus_of_Gausian_kernel=30, number_of_pixels=5, min_=0, max_=1)
	p_vect_template= [0.0606728,0.0610023,0.0607978,0.0600647,0.0588224,0.0619829,0.0623218  ,0.0621152  ,0.0613686  ,0.0601016  ,0.0627679  ,0.0631134  ,0.0629066  ,0.0621528  ,0.0608719  ,0.0630073  ,0.0633564  ,0.0631511  ,0.0623968  ,0.0611132  ,0.0626947  ,0.0630445  ,0.0628425  ,0.0620941  ,0.060819]
	p_vect = p.vectorize(0)
	for i in range(0, len(p_vect)):
		assert math.fabs(p_vect_template[i] - p_vect[i]) < 0.0000005
  


def test_distance_for_heat_maps():  
	p = gudhi.PersistenceHeatMaps(file_with_intervals="data/file_with_diagram", how_many_pixels_raidus_of_Gausian_kernel=30, min_=0, max_=1)
	q = gudhi.PersistenceHeatMaps(file_with_intervals="data/file_with_diagram_1", how_many_pixels_raidus_of_Gausian_kernel=30, min_=0, max_=1)
	r = gudhi.PersistenceHeatMaps(file_with_intervals="data/file_with_diagram_2", how_many_pixels_raidus_of_Gausian_kernel=30, min_=0, max_=1)

	epsilon = 0.0005
	assert math.fabs(p.distance(p) - 0.) < epsilon
	assert math.fabs(p.distance(q) - 624.183) < epsilon
	assert math.fabs(p.distance(r) - 415.815) < epsilon
	assert math.fabs(q.distance(p) - 624.183) < epsilon
	assert math.fabs(q.distance(q) - 0.) < epsilon
	assert math.fabs(q.distance(r) - 528.066) < epsilon
	assert math.fabs(r.distance(p) - 415.815) < epsilon
	assert math.fabs(r.distance(q) - 528.066) < epsilon
	assert math.fabs(r.distance(r) - 0.) < epsilon


def test_projections_to_R_for_heat_maps():
	p = gudhi.PersistenceHeatMaps(file_with_intervals="data/file_with_diagram",  how_many_pixels_raidus_of_Gausian_kernel=30, min_=0, max_=1)
	q = gudhi.PersistenceHeatMaps(file_with_intervals="data/file_with_diagram_1",  how_many_pixels_raidus_of_Gausian_kernel=30, min_=0, max_=1)
	r = gudhi.PersistenceHeatMaps(file_with_intervals="data/file_with_diagram_2",  how_many_pixels_raidus_of_Gausian_kernel=30, min_=0, max_=1)

	epsilon = 0.0005	
	assert math.fabs(p.project_to_R(0) - 44.3308) < epsilon
	assert math.fabs(q.project_to_R(0) - 650.568) < epsilon
	assert math.fabs(r.project_to_R(0) - 429.287) < epsilon


def test_scalar_products_for_heat_maps():
	p = gudhi.PersistenceHeatMaps(file_with_intervals="data/file_with_diagram",  how_many_pixels_raidus_of_Gausian_kernel=30, min_=0, max_=1)
	q = gudhi.PersistenceHeatMaps(file_with_intervals="data/file_with_diagram_1",  how_many_pixels_raidus_of_Gausian_kernel=30, min_=0, max_=1)
	r = gudhi.PersistenceHeatMaps(file_with_intervals="data/file_with_diagram_2",  how_many_pixels_raidus_of_Gausian_kernel=30, min_=0, max_=1)


	epsilon = 0.000005
	assert math.fabs(p.compute_scalar_product(p) - 0.0345687) < epsilon
	assert math.fabs(p.compute_scalar_product(q) - 0.0509357) < epsilon
	assert math.fabs(p.compute_scalar_product(r) - 0.0375608) < epsilon
	assert math.fabs(q.compute_scalar_product(p) - 0.0509357) < epsilon
	assert math.fabs(q.compute_scalar_product(q) - 1.31293) < epsilon
	assert math.fabs(q.compute_scalar_product(r) - 0.536799) < epsilon
	assert math.fabs(r.compute_scalar_product(p) - 0.0375608) < epsilon
	assert math.fabs(r.compute_scalar_product(q) - 0.536799) < epsilon
	assert math.fabs(r.compute_scalar_product(r) - 0.672907) < epsilon



    

