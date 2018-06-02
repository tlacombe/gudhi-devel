import gudhi
import math

"""This file is part of the Gudhi Library. The Gudhi library
   (Geometric Understanding in Higher Dimensions) is a generic C++
   library for computational topology.

   Author(s):       Pawel Dlotko

   Copyright (C) 2017 Swansea University

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
__copyright__ = "Copyright (C) 2017 Swansea University"
__license__ = "GPL v3"

epsilon = 0.0000005;
 
def test_check_min_max_function():	
	p = gudhi.PersistenceIntervals(file_with_intervals="data/file_with_diagram")
	min_max_ = p.get_x_range()
	assert math.fabs(min_max_[0] - 0.0290362) <= epsilon
	assert math.fabs(min_max_[1] - 0.994537) <= epsilon

def test_check_length_of_dominant_intervals():
	p = gudhi.PersistenceIntervals(file_with_intervals="data/file_with_diagram")
	dominant_ten_intervals_length = p.length_of_dominant_intervals(10)
	dominant_intervals_length_ = [0.8626250172000001, 0.8008928647, 0.762060545, 0.7565008858, 0.7293674164999999, 0.7181771387, 0.7083949571, 0.7028439633, 0.7004682264, 0.622176573]
	assert dominant_ten_intervals_length == dominant_intervals_length_ 


def test_check_dominant_intervals(): 
	p = gudhi.PersistenceIntervals(file_with_intervals="data/file_with_diagram");
	ten_dominant_intervals = p.dominant_intervals(10);
	templ = [(0.1147177995, 0.9773428167), (0.1336375578, 0.9345304225), (0.1045988346, 0.8666593796), (0.1497983027, 0.9062991885), (0.2473518713, 0.9767192878), (0.1926749829, 0.9108521216), (0.1918362495, 0.9002312066), (0.2849982539, 0.9878422172), (0.2940688241, 0.9945370505), (0.2674207884, 0.8895973614)]	
	assert ten_dominant_intervals == templ
	
	
def test_check_histogram_of_lengths():
	p = gudhi.PersistenceIntervals(file_with_intervals="data/file_with_diagram")
	histogram = p.histogram_of_lengths(10);
	template_histogram = [10,5,3,4,4,3,6,1,7,1,1]	
	assert histogram == template_histogram


def test_check_cumulative_histograms_of_lengths():
	p = gudhi.PersistenceIntervals(file_with_intervals="data/file_with_diagram")
	cumulative_histogram = p.cumulative_histogram_of_lengths(10)
	template_cumulative_histogram  = [10,15,18,22,26,29,35,36,43,44,45]
	assert cumulative_histogram == template_cumulative_histogram


def test_check_characteristic_function_of_diagram():
	p = gudhi.PersistenceIntervals(file_with_intervals="data/file_with_diagram")
	min_max_ = p.get_x_range();
	char_funct_diag = p.characteristic_function_of_diagram(min_max_[0], min_max_[1]);
	template_char_funct_diag = [0.3706650096057107, 0.8405801460478484, 1.2464949441848387, 1.3664031938056187, 1.3403175486249868, 1.3190409371273135, 1.1407584208620511, 0.9912592030380287, 0.8007140228196504, 0.06763026374467777]
	assert char_funct_diag == template_char_funct_diag


def test_check_cumulative_characteristic_function_of_diagram():
	p = gudhi.PersistenceIntervals(file_with_intervals="data/file_with_diagram")
	min_max_ = p.get_x_range()
	cumul_char_funct_diag = p.cumulative_characteristic_function_of_diagram(min_max_[0], min_max_[1]);
	template_char_funct_diag_cumul = [0.3706650096057107, 1.211245155653559, 2.457740099838398, 3.8241432936440165, 5.164460842269003, 6.483501779396317, 7.624260200258368, 8.615519403296396, 9.416233426116047, 9.483863689860724]
	assert cumul_char_funct_diag == template_char_funct_diag_cumul


def test_check_compute_persistent_betti_numbers():
	p = gudhi.PersistenceIntervals(file_with_intervals="data/file_with_diagram")	
	pbns = [(0.0290362441,1),(0.0307675984,2),(0.0366311713,3),(0.0544614334,4),(0.0920032915,5),(0.1045988346,6),(0.1147177995,7),(0.1173791633,8),(0.1234925873,9),(0.1336375578,10),(0.1377982685,9),(0.1497983027,10),(0.1554214768,11),(0.1584430523,12),(0.1769560101,13),(0.1832335867,12),(0.1910692175,13),(0.1913332292,14),(0.1918362495,15),(0.1926749829,16),(0.2085643625,17),(0.2184254499,18),(0.2199019471,17),(0.2323297802,16),(0.2345582605,17),(0.2371657227,16),(0.2473518713,17),(0.2674207884,18),(0.2680930813,19),(0.2787338524,18),(0.2847216921,19),(0.2849982539,20),(0.2940688241,21),(0.3062930094,22),(0.3223612854,21),(0.3231522712,22),(0.3710208764,23),(0.3723952303,24),(0.3877435371,25),(0.4355365708,26),(0.4629111887,25),(0.483568935,26),(0.4892093645,25),(0.5171153115,24),(0.5221970463,23),(0.5326650707,22),(0.545262265,23),(0.5872268199,22),(0.5930361075,23),(0.6026467816,24),(0.605044048,25),(0.6219619742,24),(0.6294493901,23),(0.636718961,22),(0.649569716,21),(0.650780786,22),(0.6549508263,23),(0.6834892731,24),(0.6871715707,23),(0.6970299556,22),(0.7011736326,21),(0.7176230095,22),(0.7220228854,21),(0.7222979715,20),(0.7253469135,19),(0.730709507,18),(0.7583554041,17),(0.7709131132,18),(0.7908326657,17),(0.8212106009,16),(0.8493046183,17),(0.8536686709,16),(0.8666593796,15),(0.8728962245,16),(0.8895973614,15),(0.9002312066,14),(0.9038473698,13),(0.9062991885,12),(0.9108521216,11),(0.9345304225,10),(0.9447571123,9),(0.9478115036,8),(0.9591544003,7),(0.9756544009,6),(0.9767192878,5),(0.9773428167,4),(0.9801285535,3),(0.9878422172,2),(0.9901270741,1),(0.9945370505,0)]
	pbns_new = p.compute_persistent_betti_numbers();  
	assert pbns == pbns_new
    

