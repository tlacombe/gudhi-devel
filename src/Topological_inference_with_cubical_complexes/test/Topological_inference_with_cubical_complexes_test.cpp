/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Pawel Dlotko
 *
 *    Copyright (C) 2017  Swansea University, UK
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

/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Pawel Dlotko
 *
 *    Copyright (C) 2015  INRIA Saclay (France)
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

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE "cubical_complex"
#include <boost/test/unit_test.hpp>

#include <gudhi/reader_utils.h>
#include <gudhi/Off_reader.h>

#include <gudhi/functions_for_topological_inference/functions_for_topological_inference.h>
#include <gudhi/Topological_inference.h>
#include <gudhi/Bitmap_cubical_complex.h>
#include <gudhi/Persistent_cohomology.h>

// standard stuff
#include <iostream>
#include <sstream>
#include <vector>

//first a few tests for functions_for_topological_inference

double epsilon = 0.0000001;


BOOST_AUTO_TEST_CASE(check_distance_functions) 
{
  std::vector< double > point1(3);
  std::vector< double > point2(3);
  
  point1[0] = point1[1] = point1[2] = 0;
  point2[0] = point2[1] = point2[2] = 1; 
  
  Gudhi::Topological_inference_with_cubical_complexes::Euclidan_distance_squared eu;
  Gudhi::Topological_inference_with_cubical_complexes::Manhattan_distance manh;
  Gudhi::Topological_inference_with_cubical_complexes::Max_norm_distance max_norm_dist;
  
  BOOST_CHECK( fabs(eu(point1,point2)-3) <= epsilon );
  BOOST_CHECK( fabs(manh(point1,point2)-3) <= epsilon );  
  BOOST_CHECK( fabs(max_norm_dist(point1,point2)-1) <= epsilon );  
}

BOOST_AUTO_TEST_CASE(kernels_centerd_in_point_cloud_test)
{
	std::vector< double > point1(2);
    std::vector< double > point2(2);
    std::vector< double > point3(2);
    std::vector< double > point4(2);
    
    point1[0] = 0;point1[1] = 1;
    point2[0] = 0;point2[1] = -1;
    point3[0] = -1;point3[1] = 0;
    point4[0] = 1;point4[1] = 0;
    
    std::vector< std::vector< double > > point_cloud(4);
    point_cloud[0] = point1;
    point_cloud[1] = point2;
    point_cloud[2] = point3;
    point_cloud[3] = point4;
        
    Gudhi::Topological_inference_with_cubical_complexes::Euclidan_distance_squared eu;
    
	Gudhi::Topological_inference_with_cubical_complexes::kernels_centerd_in_point_cloud< Gudhi::Topological_inference_with_cubical_complexes::Euclidan_distance_squared > 
	kernel( eu , point_cloud );
	
	std::vector< double > test_point(2);
	test_point[0] = test_point[1] = 0;
	BOOST_CHECK( kernel( test_point ) == 4 );
	
	test_point[0] = -10;test_point[1] = 0;	
	BOOST_CHECK( kernel( test_point ) == 404 );
	
	test_point[0] = -1;test_point[1] = 0;	
	BOOST_CHECK( kernel( test_point ) == 8 );
}


BOOST_AUTO_TEST_CASE(Sum_of_distances_from_points_test)
{
	std::vector< double > point1(2);
    std::vector< double > point2(2);
    std::vector< double > point3(2);
    std::vector< double > point4(2);
    
    point1[0] = 0;point1[1] = 1;
    point2[0] = 0;point2[1] = -1;
    point3[0] = -1;point3[1] = 0;
    point4[0] = 1;point4[1] = 0;
    
    std::vector< std::vector< double > > point_cloud(4);
    point_cloud[0] = point1;
    point_cloud[1] = point2;
    point_cloud[2] = point3;
    point_cloud[3] = point4;
    
	Gudhi::Topological_inference_with_cubical_complexes::Sum_of_distances_from_points kernel( point_cloud );
	
	std::vector< double > test_point(2);
	test_point[0] = test_point[1] = 0;
	BOOST_CHECK( kernel( test_point ) == 4 );
	
	test_point[0] = -10;test_point[1] = 0;	
	BOOST_CHECK( kernel( test_point ) == 404 );
	
	test_point[0] = -1;test_point[1] = 0;	
	BOOST_CHECK( kernel( test_point ) == 8 );
}


 
BOOST_AUTO_TEST_CASE(Distance_to_k_th_closest_point_test)
{
	std::vector< double > point1(2);
    std::vector< double > point2(2);
    std::vector< double > point3(2);
    std::vector< double > point4(2);
    
    point1[0] = 0;point1[1] = 1;
    point2[0] = 0;point2[1] = -1;
    point3[0] = -1;point3[1] = 0;
    point4[0] = 1;point4[1] = 0;
    
    std::vector< std::vector< double > > point_cloud(4);
    point_cloud[0] = point1;
    point_cloud[1] = point2;
    point_cloud[2] = point3;
    point_cloud[3] = point4;
    
    Gudhi::Topological_inference_with_cubical_complexes::Euclidan_distance_squared eu;
    
	Gudhi::Topological_inference_with_cubical_complexes::Distance_to_k_th_closest_point<Gudhi::Topological_inference_with_cubical_complexes::Euclidan_distance_squared> kernel1( point_cloud,eu,1 );
	Gudhi::Topological_inference_with_cubical_complexes::Distance_to_k_th_closest_point<Gudhi::Topological_inference_with_cubical_complexes::Euclidan_distance_squared> kernel2( point_cloud,eu,2 );
	Gudhi::Topological_inference_with_cubical_complexes::Distance_to_k_th_closest_point<Gudhi::Topological_inference_with_cubical_complexes::Euclidan_distance_squared> kernel3( point_cloud,eu,3 );
	Gudhi::Topological_inference_with_cubical_complexes::Distance_to_k_th_closest_point<Gudhi::Topological_inference_with_cubical_complexes::Euclidan_distance_squared> kernel4( point_cloud,eu,4 );
	Gudhi::Topological_inference_with_cubical_complexes::Distance_to_k_th_closest_point<Gudhi::Topological_inference_with_cubical_complexes::Euclidan_distance_squared> kernel5( point_cloud,eu,5 );	
	Gudhi::Topological_inference_with_cubical_complexes::Distance_to_k_th_closest_point<Gudhi::Topological_inference_with_cubical_complexes::Euclidan_distance_squared> kernel6( point_cloud,eu,6 );
	
	
	std::vector< double > test_point(2);
	test_point[0] = test_point[1] = 0;
	BOOST_CHECK( kernel1( test_point ) == 1 );
	BOOST_CHECK( kernel2( test_point ) == 1 );//here we get  memory access violation at address: 0x00000001: no mapping at fault address error. 
	BOOST_CHECK( kernel3( test_point ) == 1 );
	BOOST_CHECK( kernel4( test_point ) == 1 );
	BOOST_CHECK( kernel5( test_point ) == std::numeric_limits<double>::max() );
	BOOST_CHECK( kernel6( test_point ) == std::numeric_limits<double>::max() );
	
}



BOOST_AUTO_TEST_CASE(Distance_to_k_th_closest_point_test_2)
{
	std::vector< double > point1(2);
    std::vector< double > point2(2);
    std::vector< double > point3(2);
    std::vector< double > point4(2);
    
    point1[0] = 0;point1[1] = 1;
    point2[0] = 0;point2[1] = 2;
    point3[0] = 0;point3[1] = 3;
    point4[0] = 0;point4[1] = 4;
    
    std::vector< std::vector< double > > point_cloud(4);
    point_cloud[0] = point1;
    point_cloud[1] = point2;
    point_cloud[2] = point3;
    point_cloud[3] = point4;
    
    Gudhi::Topological_inference_with_cubical_complexes::Euclidan_distance_squared eu;
    
	Gudhi::Topological_inference_with_cubical_complexes::Distance_to_k_th_closest_point<Gudhi::Topological_inference_with_cubical_complexes::Euclidan_distance_squared> kernel1( point_cloud,eu,1 );
	Gudhi::Topological_inference_with_cubical_complexes::Distance_to_k_th_closest_point<Gudhi::Topological_inference_with_cubical_complexes::Euclidan_distance_squared> kernel2( point_cloud,eu,2 );
	Gudhi::Topological_inference_with_cubical_complexes::Distance_to_k_th_closest_point<Gudhi::Topological_inference_with_cubical_complexes::Euclidan_distance_squared> kernel3( point_cloud,eu,3 );
	Gudhi::Topological_inference_with_cubical_complexes::Distance_to_k_th_closest_point<Gudhi::Topological_inference_with_cubical_complexes::Euclidan_distance_squared> kernel4( point_cloud,eu,4 );
	Gudhi::Topological_inference_with_cubical_complexes::Distance_to_k_th_closest_point<Gudhi::Topological_inference_with_cubical_complexes::Euclidan_distance_squared> kernel5( point_cloud,eu,5 );	
	
	std::vector< double > test_point(2);
	test_point[0] = test_point[1] = 0;
	BOOST_CHECK( kernel1( test_point ) == 1 );
	BOOST_CHECK( kernel2( test_point ) == 4 );
	BOOST_CHECK( kernel3( test_point ) == 9 );
	BOOST_CHECK( kernel4( test_point ) == 16 );
	BOOST_CHECK( kernel5( test_point ) == std::numeric_limits<double>::max() );
		
	test_point[0] = 10; test_point[1] = 0;
	BOOST_CHECK( kernel1( test_point ) == 101 );
	BOOST_CHECK( kernel2( test_point ) == 104 );
	BOOST_CHECK( kernel3( test_point ) == 109 );
	BOOST_CHECK( kernel4( test_point ) == 116 );
	BOOST_CHECK( kernel5( test_point ) == std::numeric_limits<double>::max() );	
}



BOOST_AUTO_TEST_CASE(test_of_top_inference_2d)
{
	//this is a sampling of a circle from [-1,1]^2
	std::vector< std::vector< double > > point_cloud = 
	{ {0,1}, {0.1986693308,0.9800665778}, {0.3894183423,0.921060994}, {0.5646424734,0.8253356149},
	{0.7173560909,0.6967067093}, {0.8414709848,0.5403023059}, {0.932039086,0.3623577545}, {0.98544973,0.1699671429},
	{0.999573603,-0.0291995223}, {0.9738476309,-0.2272020947}, {0.9092974268,-0.4161468365}, {0.8084964038,-0.5885011173},
	{0.6754631806,-0.7373937155}, {0.5155013718,-0.8568887534}, {0.3349881502,-0.9422223407}, {0.1411200081,-0.9899924966},
	{-0.0583741434,-0.9982947758}, {-0.255541102,-0.9667981926}, {-0.4425204433,-0.8967584163}, {-0.6118578909,-0.7909677119},
	{-0.7568024953,-0.6536436209}, {-0.8715757724,-0.4902608213}, {-0.9516020739,-0.30733287}, {-0.9936910036,-0.1121525269},
	{-0.9961646088,0.0874989834}, {-0.9589242747,0.2836621855}, {-0.8834546557,0.4685166713}, {-0.7727644876,0.6346928759},
	{-0.6312666379,0.7755658785}, {-0.4646021794,0.8855195169}, {-0.2794154982,0.9601702867}, {-0.0830894028,0.996542097},
	{0.1165492049,0.9931849188} };
	
	//std::cout << "point_cloud.size() : " << point_cloud.size() << std::endl;
	//for ( size_t i = 0 ; i != point_cloud.size() ; ++i )
	//{
	//	std::cout << point_cloud[i][0] << " " << point_cloud[i][1] << std::endl;
	//}	
	
	//creating 100 by 100 grid on the domain [-1.1]^2
	std::vector< std::pair< double,double > > coorfinates_of_grid(2);
	coorfinates_of_grid[0] = coorfinates_of_grid[1] = std::pair<double,double>( -2.0 , 2.0 );
	std::vector< unsigned > resolution_of_a_grid(2);	
	resolution_of_a_grid[0] = resolution_of_a_grid[1] = 100;
	unsigned number_of_nearest_neighbors = 5;
	
	//typedefs:
    Gudhi::Topological_inference_with_cubical_complexes::Euclidan_distance_squared eu;
    Gudhi::Topological_inference_with_cubical_complexes::Distance_to_k_th_closest_point<Gudhi::Topological_inference_with_cubical_complexes::Euclidan_distance_squared> 
    f( point_cloud ,eu ,  number_of_nearest_neighbors );
  
    typedef Gudhi::Cubical_complex::Bitmap_cubical_complex_base<double> Bitmap_cubical_complex_base;
    typedef Gudhi::Cubical_complex::Bitmap_cubical_complex<Bitmap_cubical_complex_base> Bitmap_cubical_complex;
    typedef Gudhi::Topological_inference_with_cubical_complexes::Topological_inference< Bitmap_cubical_complex , double ,   
    Gudhi::Topological_inference_with_cubical_complexes::Distance_to_k_th_closest_point<Gudhi::Topological_inference_with_cubical_complexes::Euclidan_distance_squared> > topological_inference;
  
    typedef Gudhi::persistent_cohomology::Field_Zp Field_Zp;
    typedef Gudhi::persistent_cohomology::Persistent_cohomology<topological_inference, Field_Zp> Persistent_cohomology;

    topological_inference b( coorfinates_of_grid , resolution_of_a_grid , f );
    
    b.write_to_file_with_newlines_at_the_ends_of_structure("bla");
    b.write_to_file_Perseus_format("bla_pers");

    // Compute the persistence diagram of the complex
    Persistent_cohomology pcoh(b);
    pcoh.init_coefficients(2);  // initializes the coefficient field for homology
    pcoh.compute_persistent_cohomology(0.2);
    
    //we should get the follwong two intervals:
    //0.226802 0.951097
	//0.0625682 inf
    
    std::vector< std::tuple<size_t, size_t, int> > intervals = pcoh.get_persistent_pairs();   
    for ( size_t i = 0 ; i != intervals.size() ; ++i )
    {
		if ( b.filtration(std::get<1>(intervals[i])) == std::numeric_limits<double>::infinity() )
		{
			//we get this one: 0.0625682 inf			
		    BOOST_CHECK( fabs(   b.filtration(std::get<0>(intervals[i]))-0.0625682) < 0.0000001 );		   
		}
		else
		{			
			//we get this one 0.226802 0.951097
			 BOOST_CHECK( fabs(  b.filtration(std::get<0>(intervals[i])) - 0.226802 ) < 0.000001 );
			 BOOST_CHECK( fabs(  b.filtration(std::get<1>(intervals[i])) - 0.951097 ) < 0.000001 );			 
		}
	}
}


class distnace_from_unit_sphere
{
public:
	distnace_from_unit_sphere(){}
	double operator()( const std::vector< double >& point )
	{
		double norm = 0;
		for ( size_t i = 0 ; i != point.size() ; ++i )
		{
			norm += point[i]*point[i];
		}
		norm = sqrt(norm);
		return fabs( norm-1 );
	}
};

BOOST_AUTO_TEST_CASE(test_of_top_inference_distance_funtion_to_unit_circle)
{
	std::vector< std::pair< double,double > > coorfinates_of_grid(2);
	coorfinates_of_grid[0] = coorfinates_of_grid[1] = std::pair<double,double>( -2.0 , 2.0 );
	std::vector< unsigned > resolution_of_a_grid(2);	
	resolution_of_a_grid[0] = resolution_of_a_grid[1] = 100;	
	
    distnace_from_unit_sphere dist_unit_circle;    
  
    typedef Gudhi::Cubical_complex::Bitmap_cubical_complex_base<double> Bitmap_cubical_complex_base;
    typedef Gudhi::Cubical_complex::Bitmap_cubical_complex<Bitmap_cubical_complex_base> Bitmap_cubical_complex;
    typedef Gudhi::Topological_inference_with_cubical_complexes::Topological_inference< Bitmap_cubical_complex , double ,   
    distnace_from_unit_sphere > topological_inference;
  
    typedef Gudhi::persistent_cohomology::Field_Zp Field_Zp;
    typedef Gudhi::persistent_cohomology::Persistent_cohomology<topological_inference, Field_Zp> Persistent_cohomology;

    topological_inference b( coorfinates_of_grid , resolution_of_a_grid , dist_unit_circle );
    //b.write_to_file_Perseus_format("perse");

   // Compute the persistence diagram of the complex
    Persistent_cohomology pcoh(b);
    pcoh.init_coefficients(2);  // initializes the coefficient field for homology
    pcoh.compute_persistent_cohomology(0.1);

    std::vector< std::tuple<size_t, size_t, int> > intervals = pcoh.get_persistent_pairs();   
    for ( size_t i = 0 ; i != intervals.size() ; ++i )
    {
		//std::cout << b.filtration(std::get<0>(intervals[i])) << " " << b.filtration(std::get<1>(intervals[i])) << " " << b.filtration(std::get<2>(intervals[i])) << std::endl;
		if ( b.filtration(std::get<1>(intervals[i])) == std::numeric_limits<double>::infinity() )
		{
			//we get this one: 0.00040008 inf 			
		    BOOST_CHECK( fabs(   b.filtration(std::get<0>(intervals[i]))-0.00040008) < 0.0000001 );		   
		}
		else
		{			
			//we get this one 0.0197959 0.971716 
			 BOOST_CHECK( fabs(  b.filtration(std::get<0>(intervals[i])) - 0.0197959 ) < 0.000001 );
			 BOOST_CHECK( fabs(  b.filtration(std::get<1>(intervals[i])) - 0.971716 ) < 0.000001 );	
		}
	}	
}


BOOST_AUTO_TEST_CASE(periodic_domain_distance_test_dim2)
{
	//In this test we check the prformance of periodic_domain_distance subroutine.
	//We assume that our periodic domain is [-1,1]^2.
	std::vector< std::pair< double , double > > coordinates_of_grid(2);
	coordinates_of_grid[0] = std::pair< double,double >(-1,1);
	coordinates_of_grid[1] = std::pair< double,double >(-1,1);	
	
	Gudhi::Topological_inference_with_cubical_complexes::Euclidan_distance_squared eu;
	
	Gudhi::Topological_inference_with_cubical_complexes::periodic_domain_distance
	<Gudhi::Topological_inference_with_cubical_complexes::Euclidan_distance_squared> 
	periodic_distance( coordinates_of_grid , eu );
	
	std::vector<double> point1(2);
	point1[0] = 0; point1[1] = 0;
	
	std::vector<double> point2(2);
	point2[0] = -0.999; point2[1] = -0.999;
	
	std::vector<double> point3(2);
	point3[0] = -0.999; point3[1] = 0.999;
	
	std::vector<double> point4(2);
	point4[0] = 0.999; point4[1] = -0.999;
	
	std::vector<double> point5(2);
	point5[0] = 0.999; point5[1] = 0.999;
		
	BOOST_CHECK( fabs( periodic_distance( point1,point2 ) - 1.996 ) < 0.0005 );
	
	BOOST_CHECK( fabs( periodic_distance( point1,point3 ) - 1.996 ) < 0.0005 );
	BOOST_CHECK( fabs( periodic_distance( point1,point4 ) - 1.996 ) < 0.0005 );
	BOOST_CHECK( fabs( periodic_distance( point1,point5 ) - 1.996 ) < 0.0005 );
	
	BOOST_CHECK( fabs( periodic_distance( point2,point3 ) - 4e-06 ) < 0.0005 );
	BOOST_CHECK( fabs( periodic_distance( point2,point4 ) - 4e-06 ) < 0.0005 );
	BOOST_CHECK( fabs( periodic_distance( point2,point5 ) - 4e-06 ) < 0.0005 );
	
	BOOST_CHECK( fabs( periodic_distance( point3,point4 ) - 4e-06 ) < 0.0005 );
	BOOST_CHECK( fabs( periodic_distance( point3,point5 ) - 4e-06 ) < 0.0005 );
	
	BOOST_CHECK( fabs( periodic_distance( point3,point5 ) - 4e-06 ) < 0.0005 );
	
	BOOST_CHECK( fabs( periodic_distance( point1,point1 ) - 0 ) < 0.0005 );
	BOOST_CHECK( fabs( periodic_distance( point2,point2 ) - 0 ) < 0.0005 );
	BOOST_CHECK( fabs( periodic_distance( point3,point3 ) - 0 ) < 0.0005 );
	BOOST_CHECK( fabs( periodic_distance( point4,point4 ) - 0 ) < 0.0005 );
	BOOST_CHECK( fabs( periodic_distance( point5,point5 ) - 0 ) < 0.0005 );
	
}//periodic_domain_distance_test_dim2





BOOST_AUTO_TEST_CASE(periodic_domain_distance_test_dim3)
{
	//In this test we check the prformance of periodic_domain_distance subroutine.
	//We assume that our periodic domain is [-1,1]^2.
	std::vector< std::pair< double , double > > coordinates_of_grid(3);
	coordinates_of_grid[0] = std::pair< double,double >(-1,1);
	coordinates_of_grid[1] = std::pair< double,double >(-1,1);	
	coordinates_of_grid[2] = std::pair< double,double >(-1,1);	
	
	Gudhi::Topological_inference_with_cubical_complexes::Euclidan_distance_squared eu;
	
	Gudhi::Topological_inference_with_cubical_complexes::periodic_domain_distance
	<Gudhi::Topological_inference_with_cubical_complexes::Euclidan_distance_squared> 
	periodic_distance( coordinates_of_grid , eu );
	
	std::vector<double> point1(3);
	point1[0] = 0; point1[1] = 0; point1[2] = 0;
	
	std::vector<double> point2(3);
	point2[0] = -0.999; point2[1] = -0.999; point2[2] = -0.999;
	
	std::vector<double> point3(3);
	point3[0] = -0.999; point3[1] = -0.999; point3[2] = 0.999;
	
	std::vector<double> point4(3);
	point4[0] = -0.999; point4[1] = 0.999; point4[2] = -0.999;
	
	std::vector<double> point5(3);
	point5[0] = -0.999; point5[1] = 0.999; point5[2] = 0.999;
	
	std::vector<double> point6(3);
	point6[0] = 0.999; point6[1] = 0.999; point6[2] = 0.999;
	
	std::vector<double> point7(3);
	point7[0] = 0.999; point7[1] = 0.999; point7[2] = -0.999;
	
	std::vector<double> point8(3);
	point8[0] = 0.999; point8[1] = -0.999; point8[2] = 0.999;
	
	std::vector<double> point9(3);
	point9[0] = 0.999; point9[1] = -0.999; point9[2] = -0.999;
	
	BOOST_CHECK( fabs( periodic_distance( point1,point1 ) - 0 ) < 0.0005 );
	BOOST_CHECK( fabs( periodic_distance( point2,point2 ) - 0 ) < 0.0005 );
	BOOST_CHECK( fabs( periodic_distance( point3,point3 ) - 0 ) < 0.0005 );
	BOOST_CHECK( fabs( periodic_distance( point4,point4 ) - 0 ) < 0.0005 );
	BOOST_CHECK( fabs( periodic_distance( point5,point5 ) - 0 ) < 0.0005 );
	BOOST_CHECK( fabs( periodic_distance( point6,point6 ) - 0 ) < 0.0005 );
	BOOST_CHECK( fabs( periodic_distance( point7,point7 ) - 0 ) < 0.0005 );
	BOOST_CHECK( fabs( periodic_distance( point8,point8 ) - 0 ) < 0.0005 );
	BOOST_CHECK( fabs( periodic_distance( point9,point9 ) - 0 ) < 0.0005 );
	
	BOOST_CHECK( fabs( periodic_distance( point1,point2 ) - 2.994 ) < 0.0005 );
	BOOST_CHECK( fabs( periodic_distance( point1,point3 ) - 2.994 ) < 0.0005 );
	BOOST_CHECK( fabs( periodic_distance( point1,point4 ) - 2.994 ) < 0.0005 );
	BOOST_CHECK( fabs( periodic_distance( point1,point5 ) - 2.994 ) < 0.0005 );
	BOOST_CHECK( fabs( periodic_distance( point1,point6 ) - 2.994 ) < 0.0005 );
	BOOST_CHECK( fabs( periodic_distance( point1,point7 ) - 2.994 ) < 0.0005 );
	BOOST_CHECK( fabs( periodic_distance( point1,point8 ) - 2.994 ) < 0.0005 );
	BOOST_CHECK( fabs( periodic_distance( point1,point9 ) - 2.994 ) < 0.0005 );
	
	BOOST_CHECK( fabs( periodic_distance( point2,point3 ) - 4e-06 ) < 0.0005 );
	BOOST_CHECK( fabs( periodic_distance( point2,point4 ) - 4e-06 ) < 0.0005 );
	BOOST_CHECK( fabs( periodic_distance( point2,point5 ) - 8e-06 ) < 0.0005 );
	BOOST_CHECK( fabs( periodic_distance( point2,point6 ) - 1.2e-05 ) < 0.0005 );
	BOOST_CHECK( fabs( periodic_distance( point2,point7 ) - 8e-06 ) < 0.0005 );
	BOOST_CHECK( fabs( periodic_distance( point2,point8 ) - 8e-06 ) < 0.0005 );
	BOOST_CHECK( fabs( periodic_distance( point2,point9 ) - 4e-06 ) < 0.0005 );

	
	BOOST_CHECK( fabs( periodic_distance( point3,point4 ) - 8e-06 ) < 0.0005 );
	BOOST_CHECK( fabs( periodic_distance( point3,point5 ) - 4e-06 ) < 0.0005 );
	BOOST_CHECK( fabs( periodic_distance( point3,point6 ) - 8e-06 ) < 0.0005 );
	BOOST_CHECK( fabs( periodic_distance( point3,point7 ) - 1.2e-05 ) < 0.0005 );
	BOOST_CHECK( fabs( periodic_distance( point3,point8 ) - 4e-06 ) < 0.0005 );
	BOOST_CHECK( fabs( periodic_distance( point3,point9 ) - 8e-06 ) < 0.0005 );
	
	

	 
	BOOST_CHECK( fabs(  periodic_distance( point4,point5 ) - 4e-06 ) < 0.0005 );
	BOOST_CHECK( fabs(  periodic_distance( point4,point6 ) - 8e-06 ) < 0.0005 );
	BOOST_CHECK( fabs(  periodic_distance( point4,point7 ) - 4e-06 ) < 0.0005 );
	BOOST_CHECK( fabs(  periodic_distance( point4,point8 ) - 1.2e-05 ) < 0.0005 );
	BOOST_CHECK( fabs(  periodic_distance( point4,point9 ) - 8e-06 ) < 0.0005 );
	
	BOOST_CHECK( fabs( periodic_distance( point5,point6 ) - 4e-06 ) < 0.0005 );
	BOOST_CHECK( fabs( periodic_distance( point5,point7 ) - 8e-06 ) < 0.0005 );
	BOOST_CHECK( fabs( periodic_distance( point5,point8 ) - 8e-06 ) < 0.0005 );
	BOOST_CHECK( fabs( periodic_distance( point5,point9 )  - 1.2e-05 ) < 0.0005 );
	
	BOOST_CHECK( fabs( periodic_distance( point6,point7 ) - 4e-06 ) < 0.0005 );
	BOOST_CHECK( fabs( periodic_distance( point6,point8 ) - 4e-06 ) < 0.0005 );
	BOOST_CHECK( fabs( periodic_distance( point6,point9 ) - 8e-06 ) < 0.0005 );
	
	BOOST_CHECK( fabs( periodic_distance( point7,point8 ) - 8e-06 ) < 0.0005 );
	BOOST_CHECK( fabs( periodic_distance( point7,point9 ) - 4e-06 ) < 0.0005 );
	
	BOOST_CHECK( fabs( periodic_distance( point8,point9 ) - 4e-06 ) < 0.0005 );	
}//periodic_domain_distance_test_dim3



 
BOOST_AUTO_TEST_CASE(periodic_domain_nonperiodic_function)
{
	//in this test we will create a topological inference object based on a splitted 
	//circle point cloud. When imposing periodic boundary conditions we should see
	// a circle iff they are imposed in the x direction. 
	std::vector< std::vector< double > > point_cloud = {
	{0,1},{0.0499791693,0.9987502604},{0.0998334166,0.9950041653},{0.1494381325,0.9887710779},
    {0.1986693308,0.9800665778},{0.2474039593,0.9689124217},{0.2955202067,0.9553364891},
    {0.3428978075,0.9393727128},{0.3894183423,0.921060994},{0.4349655341,0.9004471024},
    {0.4794255386,0.8775825619},{0.5226872289,0.8525245221},{0.5646424734,0.8253356149},
    {0.6051864057,0.7960837985},{0.6442176872,0.7648421873},{0.68163876,0.7316888689},
    {0.7173560909,0.6967067093},{0.7512804051,0.6599831459},{0.7833269096,0.6216099683},
    {0.8134155048,0.5816830895},{0.8414709848,0.5403023059},{0.8674232256,0.4975710479},
    {0.8912073601,0.4535961214},{0.9127639403,0.4084874409},{0.932039086,0.3623577545},
    {0.9489846194,0.3153223624},{0.9635581854,0.2674988286},{0.9757233578,0.2190066871},
    {0.98544973,0.1699671429},{0.992712991,0.1205027694},{0.9974949866,0.0707372017},
    {0.9997837642,0.0207948278},{0.999573603,-0.0291995223},{0.9968650285,-0.0791208888},
    {0.9916648105,-0.1288444943},{0.9839859469,-0.1782460556},{0.9738476309,-0.2272020947},
    {0.961275203,-0.2755902468},{0.9463000877,-0.3232895669},{0.928959715,-0.3701808314},
    {0.9092974268,-0.4161468365},{0.8873623686,-0.4610726914},{0.8632093666,-0.5048461046},
    {0.8368987908,-0.5473576655},{0.8084964038,-0.5885011173},{0.7780731969,-0.6281736227},
    {0.7457052122,-0.6662760213},{0.7114733528,-0.7027130768},{0.6754631806,-0.7373937155},
    {0.6377647021,-0.770231254},{0.5984721441,-0.8011436155},{0.5576837174,-0.8300535352},
    {0.5155013718,-0.8568887534},{0.4720305413,-0.8815821959},{0.4273798802,-0.904072142},
    {0.3816609921,-0.9243023786},{0.3349881502,-0.9422223407},{0.2874780123,-0.9577872376},
    {0.2392493292,-0.9709581651},{0.1904226474,-0.981702203},{0.1411200081,-0.9899924966},
    {0.0914646422,-0.9958083245},{0.0415806624,-0.9991351503},{2.9915927526,-0.9999646585},
    {2.9416258566,-0.9982947758},{2.8918048655,-0.9941296761},{2.8422543059,-0.9874797699},
    {2.7930980283,-0.9783616786},{2.744458898,-0.9667981926},{2.6964584873,-0.9528182146},
    {2.6492167723,-0.9364566873},{2.6028518327,-0.917754506},{2.5574795567,-0.8967584163},
    {2.5132133513,-0.8735208977},{2.4701638591,-0.8481000317},{2.4284386813,-0.8205593573},
    {2.3881421091,-0.7909677119},{2.3493748629,-0.7593990591},{2.3122338408,-0.7259323042},
    {2.2768118759,-0.6906510966},{2.2431975047,-0.6536436209},{2.2114747456,-0.6150023765},
    {2.1817228889,-0.5748239465},{2.1540162989,-0.533208756},{2.1284242276,-0.4902608213},
    {2.1050106418,-0.4460874899},{2.0838340633,-0.4007991721},{2.0649474224,-0.354509065},
    {2.0483979261,-0.30733287},{2.0342269394,-0.2593885028},{2.0224698823,-0.2107957994},
    {2.0131561415,-0.1616762164},{2.0063089964,-0.1121525269},{2.0019455612,-0.0623485146},
    {2.0000767424,-0.0123886635},{2.000707211,0.0376021529},{2.0038353912,0.0874989834},
    {2.009453464,0.1371771121},{2.0175473874,0.1865123694},{2.0280969306,0.235381443},
    {2.0410757253,0.2836621855},{2.0564513314,0.3312339202},{2.0741853177,0.3779777427},
    {2.0942333585,0.4237768177},{2.1165453443,0.4685166713},{2.1410655066,0.5120854772},
    {2.1677325578,0.5543743362},{2.1964798441,0.595277548},{2.2272355124,0.6346928759},
    {2.2599226895,0.6725218022},{2.2944596744,0.7086697743},{2.3307601427,0.743046441},
    {2.3687333621,0.7755658785},{2.4082844194,0.8061468053},{2.4493144574,0.8347127848},
    {2.4917209225,0.8611924172},{2.5353978206,0.8855195169},{2.5802359822,0.9076332791},
    {2.6261233352,0.9274784307},{2.6729451851,0.9450053693},{2.7205845018,0.9601702867},
    {2.7689222117,0.9729352783},{2.8178374957,0.9832684384},{2.8672080911,0.9911439396},
    {2.9169105972,0.996542097},{2.9668207835,0.9994494182},{3.0168139005,0.9998586364}
	};
	
	std::vector< std::pair< double,double > > coorfinates_of_grid(2);
	coorfinates_of_grid[0] = std::pair<double,double>( 0 , 3 );
	coorfinates_of_grid[1] = std::pair<double,double>( -1 , 1 );
	std::vector< unsigned > resolution_of_a_grid(2);	
	resolution_of_a_grid[0] = resolution_of_a_grid[1] = 100;
	std::vector<bool> directions_in_which_periodic_b_cond_are_to_be_imposed(2);
	directions_in_which_periodic_b_cond_are_to_be_imposed[0] = directions_in_which_periodic_b_cond_are_to_be_imposed[1] = true;
	unsigned number_of_nearest_neighbors = 5;
	
	
	
	//first we do the check in the periodic case
	{
	//typedefs:
    Gudhi::Topological_inference_with_cubical_complexes::Euclidan_distance_squared eu;
    Gudhi::Topological_inference_with_cubical_complexes::Distance_to_k_th_closest_point<Gudhi::Topological_inference_with_cubical_complexes::Euclidan_distance_squared> 
    f( point_cloud ,eu ,  number_of_nearest_neighbors );
  
    typedef Gudhi::Cubical_complex::Bitmap_cubical_complex_periodic_boundary_conditions_base<double> Periodic_bitmap_cubical_complex_base;
    typedef Gudhi::Cubical_complex::Bitmap_cubical_complex<Periodic_bitmap_cubical_complex_base> Periodic_bitmap_cubical_complex;
    typedef Gudhi::Topological_inference_with_cubical_complexes::Topological_inference< Periodic_bitmap_cubical_complex , double ,   
    Gudhi::Topological_inference_with_cubical_complexes::Distance_to_k_th_closest_point<Gudhi::Topological_inference_with_cubical_complexes::Euclidan_distance_squared> > topological_inference;
  
    typedef Gudhi::persistent_cohomology::Field_Zp Field_Zp;
    typedef Gudhi::persistent_cohomology::Persistent_cohomology<topological_inference, Field_Zp> Persistent_cohomology;

    topological_inference b( coorfinates_of_grid , resolution_of_a_grid , f , directions_in_which_periodic_b_cond_are_to_be_imposed );

    // Compute the persistence diagram of the complex
    Persistent_cohomology pcoh(b,true);
    pcoh.init_coefficients(2);  // initializes the coefficient field for homology
    pcoh.compute_persistent_cohomology(0.05);
    
    std::vector< std::tuple<size_t, size_t, int> > intervals = pcoh.get_persistent_pairs();   
    for ( size_t i = 0 ; i != intervals.size() ; ++i )
    {
		//std::cout << b.filtration(std::get<0>(intervals[i])) << " " << b.filtration(std::get<1>(intervals[i]))  << " " << b.get_dimension_of_a_cell(std::get<0>(intervals[i])) << std::endl;
		//if dimension is 0:
		if ( b.get_dimension_of_a_cell(std::get<0>(intervals[i])) == 0 )
		{
			//intervals in dimension 0 we should get are:
			//0.0100037 inf
			BOOST_CHECK( fabs(  b.filtration(std::get<0>(intervals[i])) - 0.0100037 ) < 0.000001 );
			BOOST_CHECK( b.filtration(std::get<1>(intervals[i])) == std::numeric_limits<double>::infinity() );
		}
		else
		{
			//if dimension is 1:
			if ( b.get_dimension_of_a_cell(std::get<0>(intervals[i])) == 1 )
			{
				//intervals in dimension 1 we should get are:
				//0.0150491 0.643407
				//0.0146383 inf
				//0.252084 inf
				if ( b.filtration(std::get<1>(intervals[i])) == std::numeric_limits<double>::infinity() )
				{
					BOOST_CHECK( 
								fabs(  b.filtration(std::get<0>(intervals[i])) - 0.0146383 ) < 0.000001 
								||
								fabs(  b.filtration(std::get<0>(intervals[i])) - 0.252084 ) < 0.000001 
					           );
				}
				else
				{
					BOOST_CHECK( fabs(  b.filtration(std::get<0>(intervals[i])) - 0.0150491 ) < 0.000001 );
					BOOST_CHECK( fabs(  b.filtration(std::get<1>(intervals[i])) - 0.643407 ) < 0.000001 );
				}
			}
			else
			{
				//in this case, there should bo no other option, but I will check it anyway.
				if ( b.get_dimension_of_a_cell(std::get<0>(intervals[i])) == 2 )
				{
					//intervals in dimension 1 we should get are:
					//0.964548 inf
					BOOST_CHECK( fabs(  b.filtration(std::get<0>(intervals[i])) - 0.964548 ) < 0.000001 );
					BOOST_CHECK( b.filtration(std::get<1>(intervals[i])) == std::numeric_limits<double>::infinity() );
				}
				else
				{
					//we should not get anything here, so we want the test to fail
					BOOST_CHECK(1==0);
				}
			}
		}
	}	
	}
	
	
	
	
	//then we do a similar check in the non periodic case
	{
	//typedefs:
    Gudhi::Topological_inference_with_cubical_complexes::Euclidan_distance_squared eu;
    Gudhi::Topological_inference_with_cubical_complexes::Distance_to_k_th_closest_point<Gudhi::Topological_inference_with_cubical_complexes::Euclidan_distance_squared> 
    f( point_cloud ,eu ,  number_of_nearest_neighbors );
  
    typedef Gudhi::Cubical_complex::Bitmap_cubical_complex_base<double> Periodic_bitmap_cubical_complex_base;
    typedef Gudhi::Cubical_complex::Bitmap_cubical_complex<Periodic_bitmap_cubical_complex_base> Bitmap_cubical_complex;
    typedef Gudhi::Topological_inference_with_cubical_complexes::Topological_inference< Bitmap_cubical_complex , double ,   
    Gudhi::Topological_inference_with_cubical_complexes::Distance_to_k_th_closest_point<Gudhi::Topological_inference_with_cubical_complexes::Euclidan_distance_squared> > topological_inference;
  
    typedef Gudhi::persistent_cohomology::Field_Zp Field_Zp;
    typedef Gudhi::persistent_cohomology::Persistent_cohomology<topological_inference, Field_Zp> Persistent_cohomology;

    topological_inference b( coorfinates_of_grid , resolution_of_a_grid , f , directions_in_which_periodic_b_cond_are_to_be_imposed );

    // Compute the persistence diagram of the complex
    Persistent_cohomology pcoh(b,true);
    pcoh.init_coefficients(2);  // initializes the coefficient field for homology
    pcoh.compute_persistent_cohomology(0.05);
    
    std::vector< std::tuple<size_t, size_t, int> > intervals = pcoh.get_persistent_pairs();   
    for ( size_t i = 0 ; i != intervals.size() ; ++i )
    {
		//std::cout << b.filtration(std::get<0>(intervals[i])) << " " << b.filtration(std::get<1>(intervals[i]))  << " " << b.get_dimension_of_a_cell(std::get<0>(intervals[i])) << std::endl;
		//if dimension is 0:
		if ( b.get_dimension_of_a_cell(std::get<0>(intervals[i])) == 0 )
		{
			//intervals in dimension 0 we should get are:
			//0.0100117 0.252084
			//0.0100037 inf
			if ( b.filtration(std::get<1>(intervals[i])) == std::numeric_limits<double>::infinity() )
			{
				BOOST_CHECK( fabs(  b.filtration(std::get<0>(intervals[i])) - 0.0100037 ) < 0.000001 );
			}
			else
			{
				BOOST_CHECK( fabs(  b.filtration(std::get<0>(intervals[i])) - 0.0100117 ) < 0.000001 );
				BOOST_CHECK( fabs(  b.filtration(std::get<1>(intervals[i])) - 0.252084 ) < 0.000001 );
			}
		}
		else
		{
			//we should not get anything here, so we want the test to fail
			BOOST_CHECK(1==0);	
		}
	}
    }
    
}//periodic_domain_nonperiodic_function

