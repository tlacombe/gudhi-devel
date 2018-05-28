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


def test_check_construction_of_landscape():
	import gudhi
	diagram = [(0.9062991885,0.1497983027),(0.730709507,0.7176230095),(0.3223612854,0.0544614334),
			   (0.8895973614,0.2674207884),(0.6294493901,0.1234925873),(0.1832335867,0.1769560101),
			   (0.9767192878,0.2473518713),(0.8728962245,0.9756544009),(0.5930361075,0.7011736326),
			   (0.1045988346,0.8666593796),(0.0920032915,0.2323297802),(0.2199019471,0.1554214768),
			   (0.1910692175,0.7908326657),(0.9002312066,0.1918362495),(0.1926749829,0.9108521216),
			   (0.650780786,0.7222979715),(0.3723952303,0.649569716),(0.1377982685,0.1173791633),
			   (0.1147177995,0.9773428167),(0.4892093645,0.1584430523),(0.4629111887,0.3877435371),
			   (0.7220228854,0.2345582605),(0.2847216921,0.8536686709),(0.2184254499,0.5872268199),
			   (0.2085643625,0.2787338524),(0.2680930813,0.6871715707),(0.0290362441,0.5171153115),
			   (0.8493046183,0.9447571123),(0.7709131132,0.9901270741),(0.3710208764,0.9478115036),
			   (0.5326650707,0.4355365708),(0.0307675984,0.2371657227),(0.2940688241,0.9945370505),
			   (0.9591544003,0.545262265),(0.7583554041,0.6834892731),(0.483568935,0.5221970463),
			   (0.3231522712,0.9038473698),(0.9801285535,0.6549508263),(0.1336375578,0.9345304225),
			   (0.636718961,0.605044048),(0.2849982539,0.9878422172),(0.1913332292,0.7253469135),
			   (0.6026467816,0.8212106009),(0.0366311713,0.6219619742),(0.3062930094,0.6970299556)]

	land = gudhi.PersistenceLandscapes(diagram)
	#TODO and here we should use the constructor that read the landscape from a file and compare those two landscapes.
	

def test_check_computations_of_integrals_of_landscapes():
	import gudhi
	diagram = [(0.9062991885,0.1497983027),(0.730709507,0.7176230095),(0.3223612854,0.0544614334),
			   (0.8895973614,0.2674207884),(0.6294493901,0.1234925873),(0.1832335867,0.1769560101),
			   (0.9767192878,0.2473518713),(0.8728962245,0.9756544009),(0.5930361075,0.7011736326),
			   (0.1045988346,0.8666593796),(0.0920032915,0.2323297802),(0.2199019471,0.1554214768),
			   (0.1910692175,0.7908326657),(0.9002312066,0.1918362495),(0.1926749829,0.9108521216),
			   (0.650780786,0.7222979715),(0.3723952303,0.649569716),(0.1377982685,0.1173791633),
			   (0.1147177995,0.9773428167),(0.4892093645,0.1584430523),(0.4629111887,0.3877435371),
			   (0.7220228854,0.2345582605),(0.2847216921,0.8536686709),(0.2184254499,0.5872268199),
			   (0.2085643625,0.2787338524),(0.2680930813,0.6871715707),(0.0290362441,0.5171153115),
			   (0.8493046183,0.9447571123),(0.7709131132,0.9901270741),(0.3710208764,0.9478115036),
			   (0.5326650707,0.4355365708),(0.0307675984,0.2371657227),(0.2940688241,0.9945370505),
			   (0.9591544003,0.545262265),(0.7583554041,0.6834892731),(0.483568935,0.5221970463),
			   (0.3231522712,0.9038473698),(0.9801285535,0.6549508263),(0.1336375578,0.9345304225),
			   (0.636718961,0.605044048),(0.2849982539,0.9878422172),(0.1913332292,0.7253469135),
			   (0.6026467816,0.8212106009),(0.0366311713,0.6219619742),(0.3062930094,0.6970299556)]	
			   
			   land = gudhi.PersistenceLandscapes(diagram)
			   			   
			   assert( land.compute_integral_of_a_level_of_a_landscape(1) == 2.3499233043206282 )
			   
			   
			   assert( land.compute_integral_of_a_level_of_a_landscape(0) == 0.21643215266245644 )
			   assert( land.compute_integral_of_a_level_of_a_landscape(1) == 0.20396934789953142 )
			   assert( land.compute_integral_of_a_level_of_a_landscape(2) == 0.18497070859639664 )
			   assert( land.compute_integral_of_a_level_of_a_landscape(3) == 0.15984838067576018 )
			   assert( land.compute_integral_of_a_level_of_a_landscape(4) == 0.1444026781903593 )
			   assert( land.compute_integral_of_a_level_of_a_landscape(5) == 0.12519460961637702 )
			   assert( land.compute_integral_of_a_level_of_a_landscape(6) == 0.11577501003153404 )
	
	
