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
    
def test_check_construction_of_PSSK():
	p = gudhi.PSSK(file_with_intervals="data/file_with_diagram", how_many_pixels_raidus_of_Gausian_kernel=100, min_=0, max_=1)
	p.print_to_file("data/pssk_from_file_with_diagram")
	q = gudhi.PSSK()
	q.load_from_file("data/persistence_heat_map_from_file_with_diagram")
	assert p.compare(q)


def test_check_averages_of_heat_maps():
	p = gudhi.PSSK(file_with_intervals="data/file_with_diagram", how_many_pixels_raidus_of_Gausian_kernel=30, min_=0, max_=10);
	q = gudhi.PSSK(file_with_intervals="data/file_with_diagram_1", how_many_pixels_raidus_of_Gausian_kernel=30, min_=0, max_=10);
	r = gudhi.PSSK(file_with_intervals="data/file_with_diagram_2", how_many_pixels_raidus_of_Gausian_kernel=30, min_=0, max_=10);

	av = gudhi.PSSK();
	av.compute_average([p,q,r]);

	template_average = gudhi.PSSK();
	template_average.load_from_file("data/template_average_of_pssk");

	assert av.compare(template_average)


def test_check_median_of_heat_maps():
	p = gudhi.PSSK(file_with_intervals="data/file_with_diagram",how_many_pixels_raidus_of_Gausian_kernel=30, min_=0, max_=1);
	q = gudhi.PSSK(file_with_intervals="data/file_with_diagram_1", how_many_pixels_raidus_of_Gausian_kernel=30, min_=0, max_=1);
	r = gudhi.PSSK(file_with_intervals="data/file_with_diagram_2", how_many_pixels_raidus_of_Gausian_kernel=30, min_=0, max_=1);

	to_compute_median = [p,q,r];
	median = gudhi.PSSK();
	median.compute_median(to_compute_median);

	template_median = gudhi.PSSK();
	template_median.load_from_file("data/template_median_of_pssk");

	assert median.compare(template_median)


def test_check_compute_percentage_of_active_of_heat_maps():
	p = gudhi.PSSK(file_with_intervals="data/file_with_diagram",how_many_pixels_raidus_of_Gausian_kernel=30, min_=0, max_=1);
	q = gudhi.PSSK(file_with_intervals="data/file_with_diagram_1",how_many_pixels_raidus_of_Gausian_kernel=30, min_=0, max_=1);
	r = gudhi.PSSK(file_with_intervals="data/file_with_diagram_2",how_many_pixels_raidus_of_Gausian_kernel=30, min_=0, max_=1);

	to_compute_percentage_of_active = [p,q,r]
	percentage_of_active = gudhi.PSSK();
	percentage_of_active.compute_percentage_of_active(to_compute_percentage_of_active, 0.1);

	template_percentage_of_active = gudhi.PSSK();
	template_percentage_of_active.load_from_file("data/template_percentage_of_active_of_pssk");

	assert percentage_of_active.compare(template_percentage_of_active)


def test_check_vectorize_for_heat_maps():
	p = gudhi.PSSK(file_with_intervals="data/file_with_diagram", how_many_pixels_raidus_of_Gausian_kernel=30, number_of_pixels=5, min_=0, max_=1);
	p_vect_template= [0.0,  -0.000980622792375895, -0.0019701161569953343, -0.0029425369316000752, -0.0038723844540576, 0.000980622792375895, 0.0, -0.0009981974959707251, -0.0019877974973262287, -0.002942849271999715, 0.0019701161569953343, 0.0009981974959707251, 0.0, -0.0009983033008608634, -0.0019705338815308084, 0.0029425369316000752, 0.0019877974973262287, 0.0009983033008608634, 0.0, -0.0009809343151245984, 0.0038723844540576, 0.002942849271999715, 0.0019705338815308084, 0.0009809343151245984, 0.0]
	p_vect = p.vectorize(0);
	for i in range(0, len(p_vect)):
		assert math.fabs(p_vect_template[i] - p_vect[i]) < 0.0000005
  


def test_check_distance_for_heat_maps():
	p = gudhi.PSSK(file_with_intervals="data/file_with_diagram", how_many_pixels_raidus_of_Gausian_kernel=30, min_=0, max_=1);
	q = gudhi.PSSK(file_with_intervals="data/file_with_diagram_1", how_many_pixels_raidus_of_Gausian_kernel=30, min_=0, max_=1);
	r = gudhi.PSSK(file_with_intervals="data/file_with_diagram_2", how_many_pixels_raidus_of_Gausian_kernel=30, min_=0, max_=1);

	epsilon = 0.0005;
	assert math.fabs(p.distance(p) - 0.) < epsilon
	assert math.fabs(p.distance(q) - 1210.5692976797534) < epsilon
	assert math.fabs(p.distance(r) - 799.3897520109073) < epsilon
	assert math.fabs(q.distance(p) - 1210.5692976797534) < epsilon
	assert math.fabs(q.distance(q) - 0.) < epsilon
	assert math.fabs(q.distance(r) - 1022.5830723976615) < epsilon
	assert math.fabs(r.distance(p) - 799.3897520109073) < epsilon
	assert math.fabs(r.distance(q) - 1022.5830723976615) < epsilon
	assert math.fabs(r.distance(r) - 0.) < epsilon


def test_check_projections_to_R_for_heat_maps():
	p = gudhi.PSSK(file_with_intervals="data/file_with_diagram",  how_many_pixels_raidus_of_Gausian_kernel=30, min_=0, max_=1);
	q = gudhi.PSSK(file_with_intervals="data/file_with_diagram_1",  how_many_pixels_raidus_of_Gausian_kernel=30, min_=0, max_=1);
	r = gudhi.PSSK(file_with_intervals="data/file_with_diagram_2",  how_many_pixels_raidus_of_Gausian_kernel=30, min_=0, max_=1);

	epsilon = 0.0005;
	assert math.fabs(p.project_to_R(0) - 5.524932004884425e-13) < epsilon
	assert math.fabs(q.project_to_R(0) - 9.202995379043659e-13) < epsilon
	assert math.fabs(r.project_to_R(0) - -7.056164654707466e-12) < epsilon


def test_check_scalar_products_for_heat_maps():
	p = gudhi.PSSK(file_with_intervals="data/file_with_diagram",  how_many_pixels_raidus_of_Gausian_kernel=30, min_=0, max_=1);
	q = gudhi.PSSK(file_with_intervals="data/file_with_diagram_1",  how_many_pixels_raidus_of_Gausian_kernel=30, min_=0, max_=1);
	r = gudhi.PSSK(file_with_intervals="data/file_with_diagram_2",  how_many_pixels_raidus_of_Gausian_kernel=30, min_=0, max_=1);


	epsilon = 0.00000005;
	assert math.fabs(p.compute_scalar_product(p) - 0.06703940113854207) < epsilon
	assert math.fabs(p.compute_scalar_product(q) - 0.09969226448499101) < epsilon
	assert math.fabs(p.compute_scalar_product(r) - 0.07276588829788255) < epsilon
	assert math.fabs(q.compute_scalar_product(p) - 0.09969226448499101) < epsilon
	assert math.fabs(q.compute_scalar_product(q) - 2.5613914413826695) < epsilon
	assert math.fabs(q.compute_scalar_product(r) -  1.0352605244716466) < epsilon
	assert math.fabs(r.compute_scalar_product(p) - 0.07276588829788255) < epsilon
	assert math.fabs(r.compute_scalar_product(q) - 1.0352605244716466) < epsilon
	assert math.fabs(r.compute_scalar_product(r) -  1.2910989358886231) < epsilon



    

