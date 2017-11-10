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


#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE "cubical_complex"
#include <boost/test/unit_test.hpp>

#include <gudhi/reader_utils.h>
#include <gudhi/Off_reader.h>

#include <gudhi/functions_for_topological_inference/functions_for_topological_inference.h>
#include <gudhi/Topological_inference.h>
#include <gudhi/Morphological_operations_cubical_complex.h>
#include <gudhi/Bitmap_cubical_complex.h>
#include <gudhi/Persistent_cohomology.h>

// standard stuff
#include <iostream>
#include <sstream>
#include <vector>

//first a few tests for functions_for_topological_inference

double epsilon = 0.0000001;

//auxiliary distance functions:
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

class distnace_from_x_axis
{
public:
	distnace_from_x_axis(){}
	double operator()( const std::vector< double >& point )
	{
		double norm = 0;
		for ( size_t i = 1 ; i != point.size() ; ++i )
		{
			norm += point[i]*point[i];
		}
		norm = sqrt(norm);
		return norm;
	}
};


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
    std::vector< std::vector< double > > point_cloud = { {0,1}, {0,-1}, {-1,0}, {1,0} };
    
    Gudhi::Topological_inference_with_cubical_complexes::Euclidan_distance_squared eu;
    
	Gudhi::Topological_inference_with_cubical_complexes::Distance_to_k_th_closest_point<Gudhi::Topological_inference_with_cubical_complexes::Euclidan_distance_squared> kernel1( point_cloud,eu,1 );
	Gudhi::Topological_inference_with_cubical_complexes::Distance_to_k_th_closest_point<Gudhi::Topological_inference_with_cubical_complexes::Euclidan_distance_squared> kernel2( point_cloud,eu,2 );
	Gudhi::Topological_inference_with_cubical_complexes::Distance_to_k_th_closest_point<Gudhi::Topological_inference_with_cubical_complexes::Euclidan_distance_squared> kernel3( point_cloud,eu,3 );
	Gudhi::Topological_inference_with_cubical_complexes::Distance_to_k_th_closest_point<Gudhi::Topological_inference_with_cubical_complexes::Euclidan_distance_squared> kernel4( point_cloud,eu,4 );
	Gudhi::Topological_inference_with_cubical_complexes::Distance_to_k_th_closest_point<Gudhi::Topological_inference_with_cubical_complexes::Euclidan_distance_squared> kernel5( point_cloud,eu,5 );	
	Gudhi::Topological_inference_with_cubical_complexes::Distance_to_k_th_closest_point<Gudhi::Topological_inference_with_cubical_complexes::Euclidan_distance_squared> kernel6( point_cloud,eu,6 );
	
	
	std::vector< double > test_point = {0,0};	
	
	//std::cout << "kernel1( test_point ) : " << kernel1( test_point ) << std::endl;
	//std::cout << "kernel2( test_point ) : " << kernel2( test_point ) << std::endl;
	//std::cout << "kernel3( test_point ) : " << kernel3( test_point ) << std::endl;
	//std::cout << "kernel4( test_point ) : " << kernel4( test_point ) << std::endl;
	//std::cout << "kernel5( test_point ) : " << kernel5( test_point ) << std::endl;
	//std::cout << "kernel6( test_point ) : " << kernel6( test_point ) << std::endl;
	
	BOOST_CHECK( kernel1( test_point ) == 1 );
	BOOST_CHECK( kernel2( test_point ) == 1 );
	BOOST_CHECK( kernel3( test_point ) == 1 );
	BOOST_CHECK( kernel4( test_point ) == 1 );
	BOOST_CHECK( kernel5( test_point ) == std::numeric_limits<double>::infinity() );
	BOOST_CHECK( kernel6( test_point ) == std::numeric_limits<double>::infinity() );	
}



BOOST_AUTO_TEST_CASE(Distance_to_k_th_closest_point_k_d_tree_test)
{	
    std::vector< std::vector< double > > point_cloud = { {0,1}, {0,-1}, {-1,0}, {1,0} };
    
	Gudhi::Topological_inference_with_cubical_complexes::Distance_to_k_th_closest_point_k_d_tree kernel1( point_cloud , 1 );
	Gudhi::Topological_inference_with_cubical_complexes::Distance_to_k_th_closest_point_k_d_tree kernel2( point_cloud , 2 );
	Gudhi::Topological_inference_with_cubical_complexes::Distance_to_k_th_closest_point_k_d_tree kernel3( point_cloud , 3 );
	Gudhi::Topological_inference_with_cubical_complexes::Distance_to_k_th_closest_point_k_d_tree kernel4( point_cloud , 4 );
	Gudhi::Topological_inference_with_cubical_complexes::Distance_to_k_th_closest_point_k_d_tree kernel5( point_cloud , 5 );
	Gudhi::Topological_inference_with_cubical_complexes::Distance_to_k_th_closest_point_k_d_tree kernel6( point_cloud , 6 );
	
	
	std::vector< double > test_point = {0,0};	
	
	//std::cout << "kernel1( test_point ) : " << kernel1( test_point ) << std::endl;
	//std::cout << "kernel2( test_point ) : " << kernel2( test_point ) << std::endl;
	//std::cout << "kernel3( test_point ) : " << kernel3( test_point ) << std::endl;
	//std::cout << "kernel4( test_point ) : " << kernel4( test_point ) << std::endl;
	//std::cout << "kernel5( test_point ) : " << kernel5( test_point ) << std::endl;
	//std::cout << "kernel6( test_point ) : " << kernel6( test_point ) << std::endl;
	
	BOOST_CHECK( kernel1( test_point ) == 1 );
	BOOST_CHECK( kernel2( test_point ) == 1 );
	BOOST_CHECK( kernel3( test_point ) == 1 );
	BOOST_CHECK( kernel4( test_point ) == 1 );
	BOOST_CHECK( kernel5( test_point ) == std::numeric_limits<double>::infinity() );
	BOOST_CHECK( kernel6( test_point ) == std::numeric_limits<double>::infinity() );	
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
	BOOST_CHECK( kernel5( test_point ) == std::numeric_limits<double>::infinity() );
		
	test_point[0] = 10; test_point[1] = 0;
	BOOST_CHECK( kernel1( test_point ) == 101 );
	BOOST_CHECK( kernel2( test_point ) == 104 );
	BOOST_CHECK( kernel3( test_point ) == 109 );	
	BOOST_CHECK( kernel4( test_point ) == 116 );
	BOOST_CHECK( kernel5( test_point ) == std::numeric_limits<double>::infinity() );	
}



BOOST_AUTO_TEST_CASE(Distance_to_k_th_closest_point_k_d_tree_test_2)
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
    
	Gudhi::Topological_inference_with_cubical_complexes::Distance_to_k_th_closest_point_k_d_tree kernel1( point_cloud,1 );
	Gudhi::Topological_inference_with_cubical_complexes::Distance_to_k_th_closest_point_k_d_tree kernel2( point_cloud,2 );
	Gudhi::Topological_inference_with_cubical_complexes::Distance_to_k_th_closest_point_k_d_tree kernel3( point_cloud,3 );
	Gudhi::Topological_inference_with_cubical_complexes::Distance_to_k_th_closest_point_k_d_tree kernel4( point_cloud,4 );
	Gudhi::Topological_inference_with_cubical_complexes::Distance_to_k_th_closest_point_k_d_tree kernel5( point_cloud,5 );	
	
	std::vector< double > test_point(2);
	test_point[0] = test_point[1] = 0;
	BOOST_CHECK( kernel1( test_point ) == 1 );
	BOOST_CHECK( kernel2( test_point ) == 4 );
	BOOST_CHECK( kernel3( test_point ) == 9 );
	BOOST_CHECK( kernel4( test_point ) == 16 );
	BOOST_CHECK( kernel5( test_point ) == std::numeric_limits<double>::infinity() );
		
	test_point[0] = 10; test_point[1] = 0;
	BOOST_CHECK( kernel1( test_point ) == 101 );
	BOOST_CHECK( kernel2( test_point ) == 104 );
	BOOST_CHECK( kernel3( test_point ) == 109 );
	BOOST_CHECK( kernel4( test_point ) == 116 );
	BOOST_CHECK( kernel5( test_point ) == std::numeric_limits<double>::infinity() );	
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



BOOST_AUTO_TEST_CASE(test_of_top_inference_with_k_d_trees_2d)
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
    Gudhi::Topological_inference_with_cubical_complexes::Distance_to_k_th_closest_point_k_d_tree
    f( point_cloud ,  number_of_nearest_neighbors );
  
    typedef Gudhi::Cubical_complex::Bitmap_cubical_complex_base<double> Bitmap_cubical_complex_base;
    typedef Gudhi::Cubical_complex::Bitmap_cubical_complex<Bitmap_cubical_complex_base> Bitmap_cubical_complex;
    typedef Gudhi::Topological_inference_with_cubical_complexes::Topological_inference< Bitmap_cubical_complex , double ,   
    Gudhi::Topological_inference_with_cubical_complexes::Distance_to_k_th_closest_point_k_d_tree > topological_inference;
  
    typedef Gudhi::persistent_cohomology::Field_Zp Field_Zp;
    typedef Gudhi::persistent_cohomology::Persistent_cohomology<topological_inference, Field_Zp> Persistent_cohomology;

    topological_inference b( coorfinates_of_grid , resolution_of_a_grid , f );   

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








BOOST_AUTO_TEST_CASE(periodic_domain_nonperiodic_function_k_d_tree)
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
    
    Gudhi::Topological_inference_with_cubical_complexes::Distance_to_k_th_closest_point_k_d_tree
    f( point_cloud,  number_of_nearest_neighbors );
  
    typedef Gudhi::Cubical_complex::Bitmap_cubical_complex_periodic_boundary_conditions_base<double> Periodic_bitmap_cubical_complex_base;
    typedef Gudhi::Cubical_complex::Bitmap_cubical_complex<Periodic_bitmap_cubical_complex_base> Periodic_bitmap_cubical_complex;
    typedef Gudhi::Topological_inference_with_cubical_complexes::Topological_inference< Periodic_bitmap_cubical_complex , double ,   
    Gudhi::Topological_inference_with_cubical_complexes::Distance_to_k_th_closest_point_k_d_tree > topological_inference;
  
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
    Gudhi::Topological_inference_with_cubical_complexes::Distance_to_k_th_closest_point_k_d_tree
    f( point_cloud,  number_of_nearest_neighbors );
  
    typedef Gudhi::Cubical_complex::Bitmap_cubical_complex_base<double> Periodic_bitmap_cubical_complex_base;
    typedef Gudhi::Cubical_complex::Bitmap_cubical_complex<Periodic_bitmap_cubical_complex_base> Bitmap_cubical_complex;
    typedef Gudhi::Topological_inference_with_cubical_complexes::Topological_inference< Bitmap_cubical_complex , double ,   
    Gudhi::Topological_inference_with_cubical_complexes::Distance_to_k_th_closest_point_k_d_tree > topological_inference;
  
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
    
}//periodic_domain_nonperiodic_function_k_d_tree







BOOST_AUTO_TEST_CASE(periodic_domain_function_having_minimum_on_x_axis)
{
	std::vector< std::pair< double,double > > coorfinates_of_grid(2);
	coorfinates_of_grid[0] = coorfinates_of_grid[1] = std::pair<double,double>( -2.0 , 2.0 );
	std::vector< unsigned > resolution_of_a_grid(2);	
	resolution_of_a_grid[0] = resolution_of_a_grid[1] = 100;
	std::vector<bool> directions_in_which_periodic_b_cond_are_to_be_imposed(2);
	directions_in_which_periodic_b_cond_are_to_be_imposed[0] = true;
	directions_in_which_periodic_b_cond_are_to_be_imposed[1] = true;
		
	
    distnace_from_x_axis dist_x_axis;    
  
    typedef Gudhi::Cubical_complex::Bitmap_cubical_complex_periodic_boundary_conditions_base<double> Periodic_bitmap_cubical_complex_base;
    typedef Gudhi::Cubical_complex::Bitmap_cubical_complex<Periodic_bitmap_cubical_complex_base> Periodic_bitmap_cubical_complex;
    typedef Gudhi::Topological_inference_with_cubical_complexes::Topological_inference< Periodic_bitmap_cubical_complex , double ,   
    distnace_from_x_axis > topological_inference;
  
    typedef Gudhi::persistent_cohomology::Field_Zp Field_Zp;
    typedef Gudhi::persistent_cohomology::Persistent_cohomology<topological_inference, Field_Zp> Persistent_cohomology;

    topological_inference b( coorfinates_of_grid , resolution_of_a_grid , dist_x_axis , directions_in_which_periodic_b_cond_are_to_be_imposed );
    
   // Compute the persistence diagram of the complex
    Persistent_cohomology pcoh(b,true);
    pcoh.init_coefficients(2);  // initializes the coefficient field for homology
    pcoh.compute_persistent_cohomology(0);    

    std::vector< std::tuple<size_t, size_t, int> > intervals = pcoh.get_persistent_pairs();   
    for ( size_t i = 0 ; i != intervals.size() ; ++i )
    {
		//std::cout << b.filtration(std::get<0>(intervals[i])) << " " << b.filtration(std::get<1>(intervals[i])) << " " << b.get_dimension_of_a_cell(std::get<0>(intervals[i])) << std::endl;
		if ( b.get_dimension_of_a_cell(std::get<0>(intervals[i])) == 0 )
		{
			//We expect to see:
			//0.02 inf 0
			BOOST_CHECK( fabs(  b.filtration(std::get<0>(intervals[i])) - 0.02 ) < 0.001 );
			BOOST_CHECK( b.filtration(std::get<1>(intervals[i])) == std::numeric_limits<double>::infinity() );
		}
		else
		{
			if ( b.get_dimension_of_a_cell(std::get<0>(intervals[i])) == 1 )
			{
				//We expect to see:
				//0.02 inf 1
				//1.98 inf 1
				BOOST_CHECK( 
							(fabs(  b.filtration(std::get<0>(intervals[i])) - 0.02 ) < 0.001)
							||
							(fabs(  b.filtration(std::get<0>(intervals[i])) - 1.98 ) < 0.001)
				           );
				BOOST_CHECK( b.filtration(std::get<1>(intervals[i])) == std::numeric_limits<double>::infinity() );
			}
			else
			{
				if ( b.get_dimension_of_a_cell(std::get<0>(intervals[i])) == 2  )
				{
					//We expect to see:
					//1.98 inf 2
					BOOST_CHECK( fabs(  b.filtration(std::get<0>(intervals[i])) - 1.98 ) < 0.001 );
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




BOOST_AUTO_TEST_CASE(Morphological_operations_cubical_complex_dylation_test_1_2d)
{
	//This is the cubical complex we are going to consider:
	// 2  2  2  2 2  
	// 2 -1 -1 -1 2
	// 2 -1  2 -1 2
	// 2 -1 -1 -1 2   
	// 2  2  2  2 2
	std::vector< unsigned > sizes(2);
	sizes[0] = sizes[1] = 5;
	std::vector< double > top_dimensional_cells = {
		                                          2, 2, 2, 2,2,
		                                          2,-1,-1,-1,2,   
		                                          2,-1, 2,-1,2,
		                                          2,-1,-1,-1,2,
		                                          2, 2, 2, 2,2   
		                                         };	
	typedef Gudhi::Cubical_complex::Bitmap_cubical_complex_base<double> Bitmap_cubical_complex_base;
    typedef Gudhi::Cubical_complex::Bitmap_cubical_complex<Bitmap_cubical_complex_base> Bitmap_cubical_complex;
    typedef Gudhi::Topological_inference_with_cubical_complexes::Filtration_below_certain_value<double> Predictor_type;
    typedef Gudhi::Topological_inference_with_cubical_complexes::Morphological_operations_cubical_complex<Bitmap_cubical_complex,Predictor_type> MOCC;
    
    Bitmap_cubical_complex cmplx( sizes,top_dimensional_cells );
    Predictor_type pred(0);//cutoff value at zero.
        
	MOCC mor( &cmplx , pred );
	
	
	std::vector< double > after_runninng_predicate = {std::numeric_limits<double>::infinity(),std::numeric_limits<double>::infinity(),
	std::numeric_limits<double>::infinity(),std::numeric_limits<double>::infinity(),std::numeric_limits<double>::infinity(),std::numeric_limits<double>::infinity(),
	0,0,0,std::numeric_limits<double>::infinity(),std::numeric_limits<double>::infinity(),0,std::numeric_limits<double>::infinity(),0,std::numeric_limits<double>::infinity(),
	std::numeric_limits<double>::infinity(),0,0,0,std::numeric_limits<double>::infinity(),std::numeric_limits<double>::infinity(),std::numeric_limits<double>::infinity(),
	std::numeric_limits<double>::infinity(),std::numeric_limits<double>::infinity(),std::numeric_limits<double>::infinity()};

	
	//now iterate through all the top dimensional cells of the complex:
	size_t count = 0;
	for ( auto it = cmplx.top_dimensional_cells_iterator_begin() ; it != cmplx.top_dimensional_cells_iterator_end() ; ++it )
	{
		BOOST_CHECK( after_runninng_predicate[count] == cmplx.get_cell_data( *it ) );
		++count;
	}
	
	//now run dilation:
	using Gudhi::Topological_inference_with_cubical_complexes::considered_neighberhoods;
	mor.dilation( 1 , Gudhi::Topological_inference_with_cubical_complexes::considered_neighberhoods::full_face );
			
	std::vector<double> result_of_inference = {2,1,1,1,2,		
											   1,0,0,0,1,
											   1,0,1,0,1,
											   1,0,0,0,1, 
											   2,1,1,1,2};
	count = 0;
	for ( auto it = cmplx.top_dimensional_cells_iterator_begin() ; it != cmplx.top_dimensional_cells_iterator_end() ; ++it )
	{
		BOOST_CHECK( result_of_inference[count] == cmplx.get_cell_data( *it ) );		
		++count;
	}
	
	//compute persistence:
	typedef Gudhi::persistent_cohomology::Field_Zp Field_Zp;
    typedef Gudhi::persistent_cohomology::Persistent_cohomology<Bitmap_cubical_complex, Field_Zp> Persistent_cohomology;
    
        
	Persistent_cohomology pcoh(cmplx,true);
    pcoh.init_coefficients(2);  // initializes the coefficient field for homology
    pcoh.compute_persistent_cohomology(0.05);
    
    std::vector< std::tuple<size_t, size_t, int> > intervals = pcoh.get_persistent_pairs();  
    //we expect to have the following persistence intervals:
    //dim1 : 0 1
	//dim 0: 0 inf.
 
    for ( size_t i = 0 ; i != intervals.size() ; ++i )
    {
		//std::cout << b->filtration(std::get<0>(intervals[i])) << " " << b->filtration(std::get<1>(intervals[i]))  << " " << b->get_dimension_of_a_cell(std::get<0>(intervals[i])) << std::endl;		
		if ( cmplx.filtration(std::get<1>(intervals[i])) == std::numeric_limits< double >::infinity() )
		{
			BOOST_CHECK( cmplx.filtration(std::get<0>(intervals[i])) == 0 );
		}
		else
		{
			BOOST_CHECK( cmplx.filtration(std::get<0>(intervals[i])) == 0 );
			BOOST_CHECK( cmplx.filtration(std::get<1>(intervals[i])) == 1 );
		}
	}	
	
}//Morphological_operations_cubical_complex_test_1





BOOST_AUTO_TEST_CASE(Morphological_operations_cubical_complex_erosion_test_1_2d)
{
	//This is the cubical complex we are going to consider:
	// 2 -1 -1 -1 2  
	// 2 -1 -1 -1 2
	// 2 -1  1 -1 2
	// 2 -1 -1 -1 2   
	// 2 -1 -1 -1 2
	std::vector< unsigned > sizes(2);
	sizes[0] = sizes[1] = 5;
	std::vector< double > top_dimensional_cells = {
		                                          2,-1,-1,-1,2,
		                                          2,-1,-1,-1,2,   
		                                          2,-1,-1,-1,2,
		                                          2,-1,-1,-1,2,
		                                          2,-1,-1,-1,2   
		                                         };	
	
		                                         
	typedef Gudhi::Cubical_complex::Bitmap_cubical_complex_base<double> Bitmap_cubical_complex_base;
    typedef Gudhi::Cubical_complex::Bitmap_cubical_complex<Bitmap_cubical_complex_base> Bitmap_cubical_complex;
    typedef Gudhi::Topological_inference_with_cubical_complexes::Filtration_below_certain_value<double> Predictor_type;
    typedef Gudhi::Topological_inference_with_cubical_complexes::Morphological_operations_cubical_complex<Bitmap_cubical_complex,Predictor_type> MOCC;
    
    Bitmap_cubical_complex cmplx( sizes,top_dimensional_cells );
    Predictor_type pred(0);//cutoff value at zero.
        
	MOCC mor( &cmplx , pred );	
	
	//now run erosion:
	using Gudhi::Topological_inference_with_cubical_complexes::considered_neighberhoods;
	mor.erosion( 1 , Gudhi::Topological_inference_with_cubical_complexes::considered_neighberhoods::full_face );
	
	std::vector<double> result_of_inference = {
		std::numeric_limits< double >::infinity(),0,-1,0,std::numeric_limits< double >::infinity(),
        std::numeric_limits< double >::infinity(),0,-1,0,std::numeric_limits< double >::infinity(),
        std::numeric_limits< double >::infinity(),0,-1,0,std::numeric_limits< double >::infinity(),
        std::numeric_limits< double >::infinity(),0,-1,0,std::numeric_limits< double >::infinity(),
        std::numeric_limits< double >::infinity(),0,-1,0,std::numeric_limits< double >::infinity()
		};
	size_t count = 0;
	for ( auto it = cmplx.top_dimensional_cells_iterator_begin() ; it != cmplx.top_dimensional_cells_iterator_end() ; ++it )
	{
		BOOST_CHECK( result_of_inference[count] == cmplx.get_cell_data( *it ) );		
		//std::cout << cmplx.get_cell_data( *it ) << std::endl;
		++count;
	}
}






BOOST_AUTO_TEST_CASE(Morphological_operations_cubical_complex_both_erosion_and_dilation_test_2d)
{
	//This is the cubical complex we are going to consider:
	//2	2	2	2	2	2	2	2	2	2	2
	//2	2	2	2	2	2	2	2	2	2	2
	//2	2	3	3	3	3	3	3	3	2	2
	//2	2	3	3	3	3	3	3	3	2	2
	//2	2	3	3	3	3	3	3	3	2	2
	//2	2	3	3	3	2	3	3	3	2	2
	//2	2	3	3	3	3	3	3	3	2	2
	//2	2	3	3	3	3	3	3	3	2	2
	//2	2	3	3	3	3	3	3	3	2	2
	//2	2	2	2	2	2	2	2	2	2	2
	//2	2	2	2	2	2	2	2	2	2	2
	

	std::vector< unsigned > sizes(2);
	sizes[0] = sizes[1] = 11;
	
	std::vector< double > top_dimensional_cells = {
													2,2,2,2,2,2,2,2,2,2,2,
													2,2,2,2,2,2,2,2,2,2,2,
													2,2,3,3,3,3,3,3,3,2,2,
													2,2,3,3,3,3,3,3,3,2,2,
													2,2,3,3,3,3,3,3,3,2,2,
													2,2,3,3,3,2,3,3,3,2,2,
													2,2,3,3,3,3,3,3,3,2,2,
													2,2,3,3,3,3,3,3,3,2,2,
													2,2,3,3,3,3,3,3,3,2,2,
													2,2,2,2,2,2,2,2,2,2,2,
													2,2,2,2,2,2,2,2,2,2,2
		                                         };	
		                                         
	typedef Gudhi::Cubical_complex::Bitmap_cubical_complex_base<double> Bitmap_cubical_complex_base;
    typedef Gudhi::Cubical_complex::Bitmap_cubical_complex<Bitmap_cubical_complex_base> Bitmap_cubical_complex;
    typedef Gudhi::Topological_inference_with_cubical_complexes::Filtration_above_certain_value<double> Predictor_type;
    typedef Gudhi::Topological_inference_with_cubical_complexes::Morphological_operations_cubical_complex<Bitmap_cubical_complex,Predictor_type> MOCC;
    
    Bitmap_cubical_complex cmplx( sizes,top_dimensional_cells );
    Predictor_type pred(2.5);//cutoff value at 2.5.
        
	MOCC mor( &cmplx , pred );		
	
	double inf = std::numeric_limits< double >::infinity();	
	std::vector<double> set_and_complement =
	{
		inf,inf,inf,inf,inf,inf,inf,inf,inf,inf,inf,
		inf,inf,inf,inf,inf,inf,inf,inf,inf,inf,inf,
		inf,inf,0,0,0,0,0,0,0,inf,inf,
		inf,inf,0,0,0,0,0,0,0,inf,inf,
		inf,inf,0,0,0,0,0,0,0,inf,inf,
		inf,inf,0,0,0,inf,0,0,0,inf,inf,
		inf,inf,0,0,0,0,0,0,0,inf,inf,
		inf,inf,0,0,0,0,0,0,0,inf,inf,
		inf,inf,0,0,0,0,0,0,0,inf,inf,
		inf,inf,inf,inf,inf,inf,inf,inf,inf,inf,inf,
		inf,inf,inf,inf,inf,inf,inf,inf,inf,inf,inf 	
	};
	
	size_t count = 0;
	for ( auto it = cmplx.top_dimensional_cells_iterator_begin() ; it != cmplx.top_dimensional_cells_iterator_end() ; ++it )
	{
		BOOST_CHECK( set_and_complement[count] == cmplx.get_cell_data( *it ) );				
		++count;
	}		
	//now run erosion:
	mor.both_erosion_and_dilation( 1 , Gudhi::Topological_inference_with_cubical_complexes::considered_neighberhoods::full_face );	
	
	
	std::vector<double> result_of_inference = {	
		4,3,2, 2, 2, 2, 2, 2,2,3,4,
		3,2,1, 1, 1, 1, 1, 1,1,2,3,
		2,1,0, 0, 0, 0, 0, 0,0,1,2,
		2,1,0,-1,-1,-1,-1,-1,0,1,2,
		2,1,0,-1,-1, 0,-1,-1,0,1,2,
		2,1,0,-1, 0, 1, 0,-1,0,1,2,
		2,1,0,-1,-1, 0,-1,-1,0,1,2,
		2,1,0,-1,-1,-1,-1,-1,0,1,2,
		2,1,0, 0, 0, 0, 0, 0,0,1,2,
		3,2,1, 1, 1, 1, 1, 1,1,2,3,
		4,3,2, 2, 2, 2, 2, 2,2,3,4	
		};
		
	count = 0;
	for ( auto it = cmplx.top_dimensional_cells_iterator_begin() ; it != cmplx.top_dimensional_cells_iterator_end() ; ++it )
	{
		BOOST_CHECK( result_of_inference[count] == cmplx.get_cell_data( *it ) );		
		//std::cout << cmplx.get_cell_data( *it ) << " ";
		++count;
	}	
}

BOOST_AUTO_TEST_CASE(Morphological_operations_cubical_complex_dilation_test_2d_point_cloud)
{
	std::vector< std::vector<double> > point_cloud = 
	{
	{0.023875911,0.3819170594},{0.8273467752,0.2570094084},{0.2246383489,0.8974934225},{0.1000787409,0.5004473543},
	{0.392248322,0.1793392003},{0.0330905591,0.6366295242},{0.6505259087,0.6786510674},{0.9932294958,0.959717162},
	{0.3341948902,0.8669278021},{0.6569386947,0.8378932301},{0.6144418737,0.3400138197},{0.916686812,0.4007980435},
	{0.1678472976,0.9727493613},{0.4554182182,0.3805637609},{0.2032318886,0.7224295405},{0.7785319716,0.7035288359},
	{0.8183017471,0.3757726145},{0.6548842834,0.9780715294},{0.0597022544,0.0439436375},{0.8900984228,0.0563630168},
	{0.6750005598,0.6869571838},{0.2697740323,0.3053609708},{0.3748816792,0.5403353469},{0.8445334057,0.6003744896},
	{0.9362474468,0.8929047585},{0.0744779278,0.4765794161},{0.186476755,0.2841711254},{0.0092118601,0.6635902137},
	{0.4242709139,0.5454525957},{0.7765266201,0.7759299499},{0.8012906341,0.5699157678},{0.2628584001,0.9410507081},
	{0.2137785519,0.2910071306},{0.4552011553,0.1213308845},{0.1057428913,0.0729731533},{0.5868712259,0.3754960189},
	{0.9556682452,0.1321922166},{0.4998900071,0.6114943502},{0.1164180543,0.4240477013},{0.9319797899,0.8147350496},
	{-0.976124089,0.3819170594},{-0.1726532248,0.2570094084},{-0.7753616511,0.8974934225},{-0.8999212591,0.5004473543},
	{-0.607751678,0.1793392003},{-0.9669094409,0.6366295242},{-0.3494740913,0.6786510674},{-0.0067705042,0.959717162},
	{-0.6658051098,0.8669278021},{-0.3430613053,0.8378932301},{-0.3855581263,0.3400138197},{-0.083313188,0.4007980435},
	{-0.8321527024,0.9727493613},{-0.5445817818,0.3805637609},{-0.7967681114,0.7224295405},{-0.2214680284,0.7035288359},
	{-0.1816982529,0.3757726145},{-0.3451157166,0.9780715294},{-0.9402977456,0.0439436375},{-0.1099015772,0.0563630168},
	{-0.3249994402,0.6869571838},{-0.7302259677,0.3053609708},{-0.6251183208,0.5403353469},{-0.1554665943,0.6003744896},
	{-0.0637525532,0.8929047585},{-0.9255220722,0.4765794161},{-0.813523245,0.2841711254},{-0.9907881399,0.6635902137},
	{-0.5757290861,0.5454525957},{-0.2234733799,0.7759299499},{-0.1987093659,0.5699157678},{-0.7371415999,0.9410507081},
	{-0.7862214481,0.2910071306},{-0.5447988447,0.1213308845},{-0.8942571087,0.0729731533},{-0.4131287741,0.3754960189},
	{-0.0443317548,0.1321922166},{-0.5001099929,0.6114943502},{-0.8835819457,0.4240477013},{-0.0680202101,0.8147350496},
	{-1.976124089,0.3819170594},{-1.1726532248,0.2570094084},{-1.7753616511,0.8974934225},{-1.8999212591,0.5004473543},
	{-1.607751678,0.1793392003},{-1.9669094409,0.6366295242},{-1.3494740913,0.6786510674},{-1.0067705042,0.959717162},
	{-1.6658051098,0.8669278021},{-1.3430613053,0.8378932301},{-1.3855581263,0.3400138197},{-1.083313188,0.4007980435},
	{-1.8321527024,0.9727493613},{-1.5445817818,0.3805637609},{-1.7967681114,0.7224295405},{-1.2214680284,0.7035288359},
	{-1.1816982529,0.3757726145},{-1.3451157166,0.9780715294},{-1.9402977456,0.0439436375},{-1.1099015772,0.0563630168},
	{-1.3249994402,0.6869571838},{-1.7302259677,0.3053609708},{-1.6251183208,0.5403353469},{-1.1554665943,0.6003744896},
	{-1.0637525532,0.8929047585},{-1.9255220722,0.4765794161},{-1.813523245,0.2841711254},{-1.9907881399,0.6635902137},
	{-1.5757290861,0.5454525957},{-1.2234733799,0.7759299499},{-1.1987093659,0.5699157678},{-1.7371415999,0.9410507081},
	{-1.7862214481,0.2910071306},{-1.5447988447,0.1213308845},{-1.8942571087,0.0729731533},{-1.4131287741,0.3754960189},
	{-1.0443317548,0.1321922166},{-1.5001099929,0.6114943502},{-1.8835819457,0.4240477013},{-1.0680202101,0.8147350496},
	{0.023875911,-0.6180829406},{0.8273467752,-0.7429905916},{0.2246383489,-0.1025065775},{0.1000787409,-0.4995526457},
	{0.392248322,-0.8206607997},{0.0330905591,-0.3633704758},{0.6505259087,-0.3213489326},{0.9932294958,-0.040282838},
	{0.3341948902,-0.1330721979},{0.6569386947,-0.1621067699},{0.6144418737,-0.6599861803},{0.916686812,-0.5992019565},
	{0.1678472976,-0.0272506387},{0.4554182182,-0.6194362391},{0.2032318886,-0.2775704595},{0.7785319716,-0.2964711641},
	{0.8183017471,-0.6242273855},{0.6548842834,-0.0219284706},{0.0597022544,-0.9560563625},{0.8900984228,-0.9436369832},
	{0.6750005598,-0.3130428162},{0.2697740323,-0.6946390292},{0.3748816792,-0.4596646531},{0.8445334057,-0.3996255104},
	{0.9362474468,-0.1070952415},{0.0744779278,-0.5234205839},{0.186476755,-0.7158288746},{0.0092118601,-0.3364097863},
	{0.4242709139,-0.4545474043},{0.7765266201,-0.2240700501},{0.8012906341,-0.4300842322},{0.2628584001,-0.0589492919},
	{0.2137785519,-0.7089928694},{0.4552011553,-0.8786691155},{0.1057428913,-0.9270268467},{0.5868712259,-0.6245039811},
	{0.9556682452,-0.8678077834},{0.4998900071,-0.3885056498},{0.1164180543,-0.5759522987},{0.9319797899,-0.1852649504},
	{-1.976124089,0.3819170594},{-1.1726532248,0.2570094084},{-1.7753616511,0.8974934225},{-1.8999212591,0.5004473543},
	{-1.607751678,0.1793392003},{-1.9669094409,0.6366295242},{-1.3494740913,0.6786510674},{-1.0067705042,0.959717162},
	{-1.6658051098,0.8669278021},{-1.3430613053,0.8378932301},{-1.3855581263,0.3400138197},{-1.083313188,0.4007980435},
	{-1.8321527024,0.9727493613},{-1.5445817818,0.3805637609},{-1.7967681114,0.7224295405},{-1.2214680284,0.7035288359},
	{-1.1816982529,0.3757726145},{-1.3451157166,0.9780715294},{-1.9402977456,0.0439436375},{-1.1099015772,0.0563630168},
	{-1.3249994402,0.6869571838},{-1.7302259677,0.3053609708},{-1.6251183208,0.5403353469},{-1.1554665943,0.6003744896},
	{-1.0637525532,0.8929047585},{-1.9255220722,0.4765794161},{-1.813523245,0.2841711254},{-1.9907881399,0.6635902137},
	{-1.5757290861,0.5454525957},{-1.2234733799,0.7759299499},{-1.1987093659,0.5699157678},{-1.7371415999,0.9410507081},
	{-1.7862214481,0.2910071306},{-1.5447988447,0.1213308845},{-1.8942571087,0.0729731533},{-1.4131287741,0.3754960189},
	{-1.0443317548,0.1321922166},{-1.5001099929,0.6114943502},{-1.8835819457,0.4240477013},{-1.0680202101,0.8147350496},
	{-1.976124089,-0.6180829406},{-1.1726532248,-0.7429905916},{-1.7753616511,-0.1025065775},{-1.8999212591,-0.4995526457},
	{-1.607751678,-0.8206607997},{-1.9669094409,-0.3633704758},{-1.3494740913,-0.3213489326},{-1.0067705042,-0.040282838},
	{-1.6658051098,-0.1330721979},{-1.3430613053,-0.1621067699},{-1.3855581263,-0.6599861803},{-1.083313188,-0.5992019565},
	{-1.8321527024,-0.0272506387},{-1.5445817818,-0.6194362391},{-1.7967681114,-0.2775704595},{-1.2214680284,-0.2964711641},
	{-1.1816982529,-0.6242273855},{-1.3451157166,-0.0219284706},{-1.9402977456,-0.9560563625},{-1.1099015772,-0.9436369832},
	{-1.3249994402,-0.3130428162},{-1.7302259677,-0.6946390292},{-1.6251183208,-0.4596646531},{-1.1554665943,-0.3996255104},
	{-1.0637525532,-0.1070952415},{-1.9255220722,-0.5234205839},{-1.813523245,-0.7158288746},{-1.9907881399,-0.3364097863},
	{-1.5757290861,-0.4545474043},{-1.2234733799,-0.2240700501},{-1.1987093659,-0.4300842322},{-1.7371415999,-0.0589492919},
	{-1.7862214481,-0.7089928694},{-1.5447988447,-0.8786691155},{-1.8942571087,-0.9270268467},{-1.4131287741,-0.6245039811},
	{-1.0443317548,-0.8678077834},{-1.5001099929,-0.3885056498},{-1.8835819457,-0.5759522987},{-1.0680202101,-0.1852649504},
	{-1.976124089,-1.6180829406},{-1.1726532248,-1.7429905916},{-1.7753616511,-1.1025065775},{-1.8999212591,-1.4995526457},
	{-1.607751678,-1.8206607997},{-1.9669094409,-1.3633704758},{-1.3494740913,-1.3213489326},{-1.0067705042,-1.040282838},
	{-1.6658051098,-1.1330721979},{-1.3430613053,-1.1621067699},{-1.3855581263,-1.6599861803},{-1.083313188,-1.5992019565},
	{-1.8321527024,-1.0272506387},{-1.5445817818,-1.6194362391},{-1.7967681114,-1.2775704595},{-1.2214680284,-1.2964711641},
	{-1.1816982529,-1.6242273855},{-1.3451157166,-1.0219284706},{-1.9402977456,-1.9560563625},{-1.1099015772,-1.9436369832},
	{-1.3249994402,-1.3130428162},{-1.7302259677,-1.6946390292},{-1.6251183208,-1.4596646531},{-1.1554665943,-1.3996255104},
	{-1.0637525532,-1.1070952415},{-1.9255220722,-1.5234205839},{-1.813523245,-1.7158288746},{-1.9907881399,-1.3364097863},
	{-1.5757290861,-1.4545474043},{-1.2234733799,-1.2240700501},{-1.1987093659,-1.4300842322},{-1.7371415999,-1.0589492919},
	{-1.7862214481,-1.7089928694},{-1.5447988447,-1.8786691155},{-1.8942571087,-1.9270268467},{-1.4131287741,-1.6245039811},
	{-1.0443317548,-1.8678077834},{-1.5001099929,-1.3885056498},{-1.8835819457,-1.5759522987},{-1.0680202101,-1.1852649504},
	{-0.976124089,-1.6180829406},{-0.1726532248,-1.7429905916},{-0.7753616511,-1.1025065775},{-0.8999212591,-1.4995526457},
	{-0.607751678,-1.8206607997},{-0.9669094409,-1.3633704758},{-0.3494740913,-1.3213489326},{-0.0067705042,-1.040282838},
	{-0.6658051098,-1.1330721979},{-0.3430613053,-1.1621067699},{-0.3855581263,-1.6599861803},{-0.083313188,-1.5992019565},
	{-0.8321527024,-1.0272506387},{-0.5445817818,-1.6194362391},{-0.7967681114,-1.2775704595},{-0.2214680284,-1.2964711641},
	{-0.1816982529,-1.6242273855},{-0.3451157166,-1.0219284706},{-0.9402977456,-1.9560563625},{-0.1099015772,-1.9436369832},
	{-0.3249994402,-1.3130428162},{-0.7302259677,-1.6946390292},{-0.6251183208,-1.4596646531},{-0.1554665943,-1.3996255104},
	{-0.0637525532,-1.1070952415},{-0.9255220722,-1.5234205839},{-0.813523245,-1.7158288746},{-0.9907881399,-1.3364097863},
	{-0.5757290861,-1.4545474043},{-0.2234733799,-1.2240700501},{-0.1987093659,-1.4300842322},{-0.7371415999,-1.0589492919},
	{-0.7862214481,-1.7089928694},{-0.5447988447,-1.8786691155},{-0.8942571087,-1.9270268467},{-0.4131287741,-1.6245039811},
	{-0.0443317548,-1.8678077834},{-0.5001099929,-1.3885056498},{-0.8835819457,-1.5759522987},{-0.0680202101,-1.1852649504},
	{0.7756642951,-1.3875174306},{0.0007477794,-1.7535263395},{0.1486754336,-1.6949861362},{0.4112578223,-1.5778728577},
	{0.2404401968,-1.2271274638},{0.0103501924,-1.7481459239},{0.0296794723,-1.3997980794},{0.4371474069,-1.3633766191},
	{0.3145105685,-1.0780697165},{0.9038771724,-1.9043819481},{0.5861932419,-1.7660002813},{0.1688468996,-1.5831418226},
	{0.787674079,-1.423473675},{0.9695329207,-1.8228401768},{0.3954715841,-1.7778132861},{0.0765325117,-1.8046477444},
	{0.7355093854,-1.5665583133},{0.0726364553,-1.8255274778},{0.3235042021,-1.7918204865},{0.5873798011,-1.5158537489},
	{0.4781294588,-1.1906945852},{0.602906262,-1.78618134},{0.795998682,-1.3270501599},{0.8979716203,-1.5026051586},
	{0.5480689027,-1.9169896541},{0.8907558299,-1.7365587598},{0.4805090041,-1.2965910463},{0.5968469442,-1.917501193},
	{0.3574692875,-1.0557353233},{0.6534147975,-1.2346216682},{0.4900979213,-1.3873518452},{0.3194338754,-1.6314795793},
	{0.8625183925,-1.0589964024},{0.4055090565,-1.3075602588},{0.0782075296,-1.2786135017},{0.2948574151,-1.8746728378},
	{0.3118235911,-1.2126379462},{0.9623518405,-1.2105501306},{0.2637937481,-1.4808701964},{0.4598901793,-1.5049176225}
	};
	
	//now creat an object of topoogical inference based on it:
	std::vector< std::pair< double,double > > coorfinates_of_grid(2);
	coorfinates_of_grid[0] = coorfinates_of_grid[1] = std::pair<double,double>( -2.0 , 1 );
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
    
    
    //now build Morphological_operations_cubical_complex based on topological inference object:
    typedef Gudhi::Topological_inference_with_cubical_complexes::Filtration_below_certain_value<double> Predictor_type;
    typedef Gudhi::Topological_inference_with_cubical_complexes::Morphological_operations_cubical_complex<topological_inference,Predictor_type> MOCC;
        
    Predictor_type pred(0.1);
    
        
	MOCC mor( &b , pred );	
	
	//mor.compute_complement();
	
	mor.dilation(1 , Gudhi::Topological_inference_with_cubical_complexes::considered_neighberhoods::all );
	
	//b.write_to_file_with_newlines_at_the_ends_of_structure("annulus");
	
    
     

    // Compute the persistence diagram of the complex
    Persistent_cohomology pcoh(b);
    pcoh.init_coefficients(2);  // initializes the coefficient field for homology
    pcoh.compute_persistent_cohomology(0.2);
    
 
    
    std::vector< std::tuple<size_t, size_t, int> > intervals = pcoh.get_persistent_pairs();   
    //Here are the intervals we expect to see:
    //dimenion : 1, 0 1
	//dimenion : 1, 0 10
	//dimenion : 0, 0 inf

    for ( size_t i = 0 ; i != intervals.size() ; ++i )
    {
		//std::cout << "dimenion : " << b.get_dimension_of_a_cell(std::get<0>(intervals[i])) << " " << b.filtration(std::get<0>(intervals[i])) << " " << b.filtration(std::get<1>(intervals[i])) << std::endl;
		if ( b.get_dimension_of_a_cell(std::get<0>(intervals[i])) == 0 )
		{			
		    BOOST_CHECK( b.filtration(std::get<0>(intervals[i])) == 0) ;
		    BOOST_CHECK( b.filtration(std::get<1>(intervals[i])) == std::numeric_limits<double>::infinity()	);
	   
		}
		else
		{		
			if ( b.get_dimension_of_a_cell(std::get<0>(intervals[i])) == 1 )
		    {				
				if ( b.filtration(std::get<1>(intervals[i])) == 1 )
				{
					 BOOST_CHECK( b.filtration(std::get<0>(intervals[i]))  == 0 );
					 BOOST_CHECK( b.filtration(std::get<1>(intervals[i])) == 1 );			 
				}
				else
				{
					 BOOST_CHECK( b.filtration(std::get<0>(intervals[i]))  == 0 );
					 BOOST_CHECK( b.filtration(std::get<1>(intervals[i])) == 10 );				
				}
			}
		}
	}
}




BOOST_AUTO_TEST_CASE(Morphological_operations_cubical_complex_dilation_test_2d_point_cloud_k_d_tree)
{
	std::vector< std::vector<double> > point_cloud = 
	{
	{0.023875911,0.3819170594},{0.8273467752,0.2570094084},{0.2246383489,0.8974934225},{0.1000787409,0.5004473543},
	{0.392248322,0.1793392003},{0.0330905591,0.6366295242},{0.6505259087,0.6786510674},{0.9932294958,0.959717162},
	{0.3341948902,0.8669278021},{0.6569386947,0.8378932301},{0.6144418737,0.3400138197},{0.916686812,0.4007980435},
	{0.1678472976,0.9727493613},{0.4554182182,0.3805637609},{0.2032318886,0.7224295405},{0.7785319716,0.7035288359},
	{0.8183017471,0.3757726145},{0.6548842834,0.9780715294},{0.0597022544,0.0439436375},{0.8900984228,0.0563630168},
	{0.6750005598,0.6869571838},{0.2697740323,0.3053609708},{0.3748816792,0.5403353469},{0.8445334057,0.6003744896},
	{0.9362474468,0.8929047585},{0.0744779278,0.4765794161},{0.186476755,0.2841711254},{0.0092118601,0.6635902137},
	{0.4242709139,0.5454525957},{0.7765266201,0.7759299499},{0.8012906341,0.5699157678},{0.2628584001,0.9410507081},
	{0.2137785519,0.2910071306},{0.4552011553,0.1213308845},{0.1057428913,0.0729731533},{0.5868712259,0.3754960189},
	{0.9556682452,0.1321922166},{0.4998900071,0.6114943502},{0.1164180543,0.4240477013},{0.9319797899,0.8147350496},
	{-0.976124089,0.3819170594},{-0.1726532248,0.2570094084},{-0.7753616511,0.8974934225},{-0.8999212591,0.5004473543},
	{-0.607751678,0.1793392003},{-0.9669094409,0.6366295242},{-0.3494740913,0.6786510674},{-0.0067705042,0.959717162},
	{-0.6658051098,0.8669278021},{-0.3430613053,0.8378932301},{-0.3855581263,0.3400138197},{-0.083313188,0.4007980435},
	{-0.8321527024,0.9727493613},{-0.5445817818,0.3805637609},{-0.7967681114,0.7224295405},{-0.2214680284,0.7035288359},
	{-0.1816982529,0.3757726145},{-0.3451157166,0.9780715294},{-0.9402977456,0.0439436375},{-0.1099015772,0.0563630168},
	{-0.3249994402,0.6869571838},{-0.7302259677,0.3053609708},{-0.6251183208,0.5403353469},{-0.1554665943,0.6003744896},
	{-0.0637525532,0.8929047585},{-0.9255220722,0.4765794161},{-0.813523245,0.2841711254},{-0.9907881399,0.6635902137},
	{-0.5757290861,0.5454525957},{-0.2234733799,0.7759299499},{-0.1987093659,0.5699157678},{-0.7371415999,0.9410507081},
	{-0.7862214481,0.2910071306},{-0.5447988447,0.1213308845},{-0.8942571087,0.0729731533},{-0.4131287741,0.3754960189},
	{-0.0443317548,0.1321922166},{-0.5001099929,0.6114943502},{-0.8835819457,0.4240477013},{-0.0680202101,0.8147350496},
	{-1.976124089,0.3819170594},{-1.1726532248,0.2570094084},{-1.7753616511,0.8974934225},{-1.8999212591,0.5004473543},
	{-1.607751678,0.1793392003},{-1.9669094409,0.6366295242},{-1.3494740913,0.6786510674},{-1.0067705042,0.959717162},
	{-1.6658051098,0.8669278021},{-1.3430613053,0.8378932301},{-1.3855581263,0.3400138197},{-1.083313188,0.4007980435},
	{-1.8321527024,0.9727493613},{-1.5445817818,0.3805637609},{-1.7967681114,0.7224295405},{-1.2214680284,0.7035288359},
	{-1.1816982529,0.3757726145},{-1.3451157166,0.9780715294},{-1.9402977456,0.0439436375},{-1.1099015772,0.0563630168},
	{-1.3249994402,0.6869571838},{-1.7302259677,0.3053609708},{-1.6251183208,0.5403353469},{-1.1554665943,0.6003744896},
	{-1.0637525532,0.8929047585},{-1.9255220722,0.4765794161},{-1.813523245,0.2841711254},{-1.9907881399,0.6635902137},
	{-1.5757290861,0.5454525957},{-1.2234733799,0.7759299499},{-1.1987093659,0.5699157678},{-1.7371415999,0.9410507081},
	{-1.7862214481,0.2910071306},{-1.5447988447,0.1213308845},{-1.8942571087,0.0729731533},{-1.4131287741,0.3754960189},
	{-1.0443317548,0.1321922166},{-1.5001099929,0.6114943502},{-1.8835819457,0.4240477013},{-1.0680202101,0.8147350496},
	{0.023875911,-0.6180829406},{0.8273467752,-0.7429905916},{0.2246383489,-0.1025065775},{0.1000787409,-0.4995526457},
	{0.392248322,-0.8206607997},{0.0330905591,-0.3633704758},{0.6505259087,-0.3213489326},{0.9932294958,-0.040282838},
	{0.3341948902,-0.1330721979},{0.6569386947,-0.1621067699},{0.6144418737,-0.6599861803},{0.916686812,-0.5992019565},
	{0.1678472976,-0.0272506387},{0.4554182182,-0.6194362391},{0.2032318886,-0.2775704595},{0.7785319716,-0.2964711641},
	{0.8183017471,-0.6242273855},{0.6548842834,-0.0219284706},{0.0597022544,-0.9560563625},{0.8900984228,-0.9436369832},
	{0.6750005598,-0.3130428162},{0.2697740323,-0.6946390292},{0.3748816792,-0.4596646531},{0.8445334057,-0.3996255104},
	{0.9362474468,-0.1070952415},{0.0744779278,-0.5234205839},{0.186476755,-0.7158288746},{0.0092118601,-0.3364097863},
	{0.4242709139,-0.4545474043},{0.7765266201,-0.2240700501},{0.8012906341,-0.4300842322},{0.2628584001,-0.0589492919},
	{0.2137785519,-0.7089928694},{0.4552011553,-0.8786691155},{0.1057428913,-0.9270268467},{0.5868712259,-0.6245039811},
	{0.9556682452,-0.8678077834},{0.4998900071,-0.3885056498},{0.1164180543,-0.5759522987},{0.9319797899,-0.1852649504},
	{-1.976124089,0.3819170594},{-1.1726532248,0.2570094084},{-1.7753616511,0.8974934225},{-1.8999212591,0.5004473543},
	{-1.607751678,0.1793392003},{-1.9669094409,0.6366295242},{-1.3494740913,0.6786510674},{-1.0067705042,0.959717162},
	{-1.6658051098,0.8669278021},{-1.3430613053,0.8378932301},{-1.3855581263,0.3400138197},{-1.083313188,0.4007980435},
	{-1.8321527024,0.9727493613},{-1.5445817818,0.3805637609},{-1.7967681114,0.7224295405},{-1.2214680284,0.7035288359},
	{-1.1816982529,0.3757726145},{-1.3451157166,0.9780715294},{-1.9402977456,0.0439436375},{-1.1099015772,0.0563630168},
	{-1.3249994402,0.6869571838},{-1.7302259677,0.3053609708},{-1.6251183208,0.5403353469},{-1.1554665943,0.6003744896},
	{-1.0637525532,0.8929047585},{-1.9255220722,0.4765794161},{-1.813523245,0.2841711254},{-1.9907881399,0.6635902137},
	{-1.5757290861,0.5454525957},{-1.2234733799,0.7759299499},{-1.1987093659,0.5699157678},{-1.7371415999,0.9410507081},
	{-1.7862214481,0.2910071306},{-1.5447988447,0.1213308845},{-1.8942571087,0.0729731533},{-1.4131287741,0.3754960189},
	{-1.0443317548,0.1321922166},{-1.5001099929,0.6114943502},{-1.8835819457,0.4240477013},{-1.0680202101,0.8147350496},
	{-1.976124089,-0.6180829406},{-1.1726532248,-0.7429905916},{-1.7753616511,-0.1025065775},{-1.8999212591,-0.4995526457},
	{-1.607751678,-0.8206607997},{-1.9669094409,-0.3633704758},{-1.3494740913,-0.3213489326},{-1.0067705042,-0.040282838},
	{-1.6658051098,-0.1330721979},{-1.3430613053,-0.1621067699},{-1.3855581263,-0.6599861803},{-1.083313188,-0.5992019565},
	{-1.8321527024,-0.0272506387},{-1.5445817818,-0.6194362391},{-1.7967681114,-0.2775704595},{-1.2214680284,-0.2964711641},
	{-1.1816982529,-0.6242273855},{-1.3451157166,-0.0219284706},{-1.9402977456,-0.9560563625},{-1.1099015772,-0.9436369832},
	{-1.3249994402,-0.3130428162},{-1.7302259677,-0.6946390292},{-1.6251183208,-0.4596646531},{-1.1554665943,-0.3996255104},
	{-1.0637525532,-0.1070952415},{-1.9255220722,-0.5234205839},{-1.813523245,-0.7158288746},{-1.9907881399,-0.3364097863},
	{-1.5757290861,-0.4545474043},{-1.2234733799,-0.2240700501},{-1.1987093659,-0.4300842322},{-1.7371415999,-0.0589492919},
	{-1.7862214481,-0.7089928694},{-1.5447988447,-0.8786691155},{-1.8942571087,-0.9270268467},{-1.4131287741,-0.6245039811},
	{-1.0443317548,-0.8678077834},{-1.5001099929,-0.3885056498},{-1.8835819457,-0.5759522987},{-1.0680202101,-0.1852649504},
	{-1.976124089,-1.6180829406},{-1.1726532248,-1.7429905916},{-1.7753616511,-1.1025065775},{-1.8999212591,-1.4995526457},
	{-1.607751678,-1.8206607997},{-1.9669094409,-1.3633704758},{-1.3494740913,-1.3213489326},{-1.0067705042,-1.040282838},
	{-1.6658051098,-1.1330721979},{-1.3430613053,-1.1621067699},{-1.3855581263,-1.6599861803},{-1.083313188,-1.5992019565},
	{-1.8321527024,-1.0272506387},{-1.5445817818,-1.6194362391},{-1.7967681114,-1.2775704595},{-1.2214680284,-1.2964711641},
	{-1.1816982529,-1.6242273855},{-1.3451157166,-1.0219284706},{-1.9402977456,-1.9560563625},{-1.1099015772,-1.9436369832},
	{-1.3249994402,-1.3130428162},{-1.7302259677,-1.6946390292},{-1.6251183208,-1.4596646531},{-1.1554665943,-1.3996255104},
	{-1.0637525532,-1.1070952415},{-1.9255220722,-1.5234205839},{-1.813523245,-1.7158288746},{-1.9907881399,-1.3364097863},
	{-1.5757290861,-1.4545474043},{-1.2234733799,-1.2240700501},{-1.1987093659,-1.4300842322},{-1.7371415999,-1.0589492919},
	{-1.7862214481,-1.7089928694},{-1.5447988447,-1.8786691155},{-1.8942571087,-1.9270268467},{-1.4131287741,-1.6245039811},
	{-1.0443317548,-1.8678077834},{-1.5001099929,-1.3885056498},{-1.8835819457,-1.5759522987},{-1.0680202101,-1.1852649504},
	{-0.976124089,-1.6180829406},{-0.1726532248,-1.7429905916},{-0.7753616511,-1.1025065775},{-0.8999212591,-1.4995526457},
	{-0.607751678,-1.8206607997},{-0.9669094409,-1.3633704758},{-0.3494740913,-1.3213489326},{-0.0067705042,-1.040282838},
	{-0.6658051098,-1.1330721979},{-0.3430613053,-1.1621067699},{-0.3855581263,-1.6599861803},{-0.083313188,-1.5992019565},
	{-0.8321527024,-1.0272506387},{-0.5445817818,-1.6194362391},{-0.7967681114,-1.2775704595},{-0.2214680284,-1.2964711641},
	{-0.1816982529,-1.6242273855},{-0.3451157166,-1.0219284706},{-0.9402977456,-1.9560563625},{-0.1099015772,-1.9436369832},
	{-0.3249994402,-1.3130428162},{-0.7302259677,-1.6946390292},{-0.6251183208,-1.4596646531},{-0.1554665943,-1.3996255104},
	{-0.0637525532,-1.1070952415},{-0.9255220722,-1.5234205839},{-0.813523245,-1.7158288746},{-0.9907881399,-1.3364097863},
	{-0.5757290861,-1.4545474043},{-0.2234733799,-1.2240700501},{-0.1987093659,-1.4300842322},{-0.7371415999,-1.0589492919},
	{-0.7862214481,-1.7089928694},{-0.5447988447,-1.8786691155},{-0.8942571087,-1.9270268467},{-0.4131287741,-1.6245039811},
	{-0.0443317548,-1.8678077834},{-0.5001099929,-1.3885056498},{-0.8835819457,-1.5759522987},{-0.0680202101,-1.1852649504},
	{0.7756642951,-1.3875174306},{0.0007477794,-1.7535263395},{0.1486754336,-1.6949861362},{0.4112578223,-1.5778728577},
	{0.2404401968,-1.2271274638},{0.0103501924,-1.7481459239},{0.0296794723,-1.3997980794},{0.4371474069,-1.3633766191},
	{0.3145105685,-1.0780697165},{0.9038771724,-1.9043819481},{0.5861932419,-1.7660002813},{0.1688468996,-1.5831418226},
	{0.787674079,-1.423473675},{0.9695329207,-1.8228401768},{0.3954715841,-1.7778132861},{0.0765325117,-1.8046477444},
	{0.7355093854,-1.5665583133},{0.0726364553,-1.8255274778},{0.3235042021,-1.7918204865},{0.5873798011,-1.5158537489},
	{0.4781294588,-1.1906945852},{0.602906262,-1.78618134},{0.795998682,-1.3270501599},{0.8979716203,-1.5026051586},
	{0.5480689027,-1.9169896541},{0.8907558299,-1.7365587598},{0.4805090041,-1.2965910463},{0.5968469442,-1.917501193},
	{0.3574692875,-1.0557353233},{0.6534147975,-1.2346216682},{0.4900979213,-1.3873518452},{0.3194338754,-1.6314795793},
	{0.8625183925,-1.0589964024},{0.4055090565,-1.3075602588},{0.0782075296,-1.2786135017},{0.2948574151,-1.8746728378},
	{0.3118235911,-1.2126379462},{0.9623518405,-1.2105501306},{0.2637937481,-1.4808701964},{0.4598901793,-1.5049176225}
	};
	
	//now creat an object of topoogical inference based on it:
	std::vector< std::pair< double,double > > coorfinates_of_grid(2);
	coorfinates_of_grid[0] = coorfinates_of_grid[1] = std::pair<double,double>( -2.0 , 1 );
	std::vector< unsigned > resolution_of_a_grid(2);	
	resolution_of_a_grid[0] = resolution_of_a_grid[1] = 100;
	unsigned number_of_nearest_neighbors = 5;
	
	//typedefs:
    Gudhi::Topological_inference_with_cubical_complexes::Distance_to_k_th_closest_point_k_d_tree
    f( point_cloud ,  number_of_nearest_neighbors );
  
    typedef Gudhi::Cubical_complex::Bitmap_cubical_complex_base<double> Bitmap_cubical_complex_base;
    typedef Gudhi::Cubical_complex::Bitmap_cubical_complex<Bitmap_cubical_complex_base> Bitmap_cubical_complex;
    typedef Gudhi::Topological_inference_with_cubical_complexes::Topological_inference< Bitmap_cubical_complex , double ,   
    Gudhi::Topological_inference_with_cubical_complexes::Distance_to_k_th_closest_point_k_d_tree > topological_inference;
  
    typedef Gudhi::persistent_cohomology::Field_Zp Field_Zp;
    typedef Gudhi::persistent_cohomology::Persistent_cohomology<topological_inference, Field_Zp> Persistent_cohomology;

    topological_inference b( coorfinates_of_grid , resolution_of_a_grid , f );  
    
    
    //now build Morphological_operations_cubical_complex based on topological inference object:
    typedef Gudhi::Topological_inference_with_cubical_complexes::Filtration_below_certain_value<double> Predictor_type;
    typedef Gudhi::Topological_inference_with_cubical_complexes::Morphological_operations_cubical_complex<topological_inference,Predictor_type> MOCC;
        
    Predictor_type pred(0.1);
    
        
	MOCC mor( &b , pred );	
	
	//mor.compute_complement();
	
	mor.dilation(1 , Gudhi::Topological_inference_with_cubical_complexes::considered_neighberhoods::all );
	
	//b.write_to_file_with_newlines_at_the_ends_of_structure("annulus");
	
    
     

    // Compute the persistence diagram of the complex
    Persistent_cohomology pcoh(b);
    pcoh.init_coefficients(2);  // initializes the coefficient field for homology
    pcoh.compute_persistent_cohomology(0.2);
    
 
    
    std::vector< std::tuple<size_t, size_t, int> > intervals = pcoh.get_persistent_pairs();   
    //Here are the intervals we expect to see:
    //dimenion : 1, 0 1
	//dimenion : 1, 0 10
	//dimenion : 0, 0 inf

    for ( size_t i = 0 ; i != intervals.size() ; ++i )
    {
		//std::cout << "dimenion : " << b.get_dimension_of_a_cell(std::get<0>(intervals[i])) << " " << b.filtration(std::get<0>(intervals[i])) << " " << b.filtration(std::get<1>(intervals[i])) << std::endl;
		if ( b.get_dimension_of_a_cell(std::get<0>(intervals[i])) == 0 )
		{			
		    BOOST_CHECK( b.filtration(std::get<0>(intervals[i])) == 0) ;
		    BOOST_CHECK( b.filtration(std::get<1>(intervals[i])) == std::numeric_limits<double>::infinity()	);
	   
		}
		else
		{		
			if ( b.get_dimension_of_a_cell(std::get<0>(intervals[i])) == 1 )
		    {				
				if ( b.filtration(std::get<1>(intervals[i])) == 1 )
				{
					 BOOST_CHECK( b.filtration(std::get<0>(intervals[i]))  == 0 );
					 BOOST_CHECK( b.filtration(std::get<1>(intervals[i])) == 1 );			 
				}
				else
				{
					 BOOST_CHECK( b.filtration(std::get<0>(intervals[i]))  == 0 );
					 BOOST_CHECK( b.filtration(std::get<1>(intervals[i])) == 10 );				
				}
			}
		}
	}
}



BOOST_AUTO_TEST_CASE(Morphological_operations_cubical_complex_dilation_test_2d_periodic_point_cloud)
{
	std::vector< std::vector<double> > point_cloud = 
	{
		{0.8969653195,0.8669824321},{0.948787326,0.1550864119},{0.6808391865,0.1174833223},{0.380992071,0.685172674},
		{0.5667654248,0.7725275909},{0.7717394382,0.8585212138},{0.0832100753,0.6969464312},{0.6902645568,0.895427923},
		{0.1910863747,0.1352325361},{0.4057740739,0.7557078151},{0.3033710681,0.1300876634},{0.9716873935,0.5043000265},
		{0.1003751704,0.4719773203},{0.8046362116,0.9671608615},{0.807659992,0.2279810838},{0.2642463297,0.5731411912},
		{0.9149513317,0.0878322234},{0.1725596942,0.8915775996},{0.1494219254,0.5552247122},{0.6168100361,0.0783634349},
		{0.2417251917,0.9171568255},{0.5437938895,0.3794184935},{0.1963619224,0.5923618199},{0.2248295015,0.3775884192},
		{0.3289951666,0.1755200389},{0.4835049524,0.5809588695},{0.8769045991,0.1403569882},{0.5677420071,0.6600178697},
		{0.0566576829,0.2155916393},{0.2042365933,0.6464692687},{0.624759804,0.8349595866},{0.3600770941,0.7845088479},
		{0.0542968004,0.4673010844},{0.0068651608,0.2392006286},{0.695079451,0.6422875135},{0.6452568541,0.5775681087},
		{0.6109693649,0.429406171},{0.0660107879,0.63594877},{0.2174559585,0.2505642476},{0.3623426859,0.5123264501}
	};
	
	//now creat an object of topoogical inference based on it:
	std::vector< std::pair< double,double > > coorfinates_of_grid(2);
	coorfinates_of_grid[0] = std::pair<double,double>( -1 , 2 );
	coorfinates_of_grid[1] = std::pair<double,double>( 0 , 1 );
	std::vector< unsigned > resolution_of_a_grid(2);	
	resolution_of_a_grid[0] = resolution_of_a_grid[1] = 100;
	unsigned number_of_nearest_neighbors = 5;
	
	std::vector<bool> directions_in_which_periodic_boudary_conditions_are_to_be_imposed(2);
	directions_in_which_periodic_boudary_conditions_are_to_be_imposed[0] = true;
	directions_in_which_periodic_boudary_conditions_are_to_be_imposed[1] = true;
	
	
	//typedefs:
    Gudhi::Topological_inference_with_cubical_complexes::Euclidan_distance_squared eu;
    Gudhi::Topological_inference_with_cubical_complexes::Distance_to_k_th_closest_point<Gudhi::Topological_inference_with_cubical_complexes::Euclidan_distance_squared> 
    f( point_cloud ,eu ,  number_of_nearest_neighbors );
  
    typedef Gudhi::Cubical_complex::Bitmap_cubical_complex_periodic_boundary_conditions_base<double> Bitmap_cubical_complex_periodic_base;
    typedef Gudhi::Cubical_complex::Bitmap_cubical_complex<Bitmap_cubical_complex_periodic_base> Bitmap_cubical_complex_periodic;
    typedef Gudhi::Topological_inference_with_cubical_complexes::Topological_inference< Bitmap_cubical_complex_periodic , double ,   
    Gudhi::Topological_inference_with_cubical_complexes::Distance_to_k_th_closest_point<Gudhi::Topological_inference_with_cubical_complexes::Euclidan_distance_squared> > topological_inference;
  
    typedef Gudhi::persistent_cohomology::Field_Zp Field_Zp;
    typedef Gudhi::persistent_cohomology::Persistent_cohomology<topological_inference, Field_Zp> Persistent_cohomology;

    topological_inference b( coorfinates_of_grid , resolution_of_a_grid , f , directions_in_which_periodic_boudary_conditions_are_to_be_imposed );  
    
    
    //now build Morphological_operations_cubical_complex based on topological inference object:    
    typedef Gudhi::Topological_inference_with_cubical_complexes::Filtration_below_certain_value<double> Predictor_type;
    typedef Gudhi::Topological_inference_with_cubical_complexes::Morphological_operations_cubical_complex<topological_inference,Predictor_type> MOCC;
        
    Predictor_type pred(0.3);        
	MOCC mor( &b , pred );			
	
	mor.dilation(1 , Gudhi::Topological_inference_with_cubical_complexes::considered_neighberhoods::all );
	    

    // Compute the persistence diagram of the complex
    Persistent_cohomology pcoh(b);
    pcoh.init_coefficients(2);  // initializes the coefficient field for homology
    pcoh.compute_persistent_cohomology(2);
    
    
    std::vector< std::tuple<size_t, size_t, int> > intervals = pcoh.get_persistent_pairs();   
    //Here are the intervals we expect to see:
	//dimenion : 0 0 inf
	//dimenion : 1 0 inf
	//dimenion : 1 21 inf

    for ( size_t i = 0 ; i != intervals.size() ; ++i )
    {
		//std::cout << "dimenion : " << b.get_dimension_of_a_cell(std::get<0>(intervals[i])) << " " << b.filtration(std::get<0>(intervals[i])) << " " << b.filtration(std::get<1>(intervals[i])) << std::endl;
		if ( b.get_dimension_of_a_cell(std::get<0>(intervals[i])) == 0 )
		{			
		    BOOST_CHECK( b.filtration(std::get<0>(intervals[i])) == 0) ;
		    BOOST_CHECK( b.filtration(std::get<1>(intervals[i])) == std::numeric_limits<double>::infinity()	);
	   
		}
		else
		{		
			if ( b.get_dimension_of_a_cell(std::get<0>(intervals[i])) == 0 )
		    {				
				if ( b.filtration(std::get<0>(intervals[i])) == 1 )
				{
					 BOOST_CHECK( b.filtration(std::get<0>(intervals[i]))  == 1 );
					 BOOST_CHECK( b.filtration(std::get<1>(intervals[i])) == std::numeric_limits<double>::infinity() );			 
				}
				else
				{								
					 BOOST_CHECK( fabs(b.filtration(std::get<0>(intervals[i]))  - 21) < 0.0001 );
					 BOOST_CHECK( b.filtration(std::get<1>(intervals[i])) == std::numeric_limits<double>::infinity() );				
				}
			}
		}
	}
}




BOOST_AUTO_TEST_CASE(Distance_to_k_th_nearest_neighbor_periodic_domain_2d_brute_force)
{
	//first we test the brute force algorithm:

	typedef Gudhi::Topological_inference_with_cubical_complexes::Euclidan_distance_squared Euclidean_distance_squared;		
	typedef Gudhi::Topological_inference_with_cubical_complexes::periodic_domain_distance<Euclidean_distance_squared>
	periodic_Euclidean_distance_squared;
	
	
	std::vector< std::pair< double , double > > coordinates_of_grid(2);
	coordinates_of_grid[0] = std::pair< double,double >(-1,1);
	coordinates_of_grid[1] = std::pair< double,double >(-1,1);
	
	Euclidean_distance_squared eu;
	periodic_Euclidean_distance_squared period_eu( coordinates_of_grid , eu );
	
	std::vector<double> point1 = {-0.999,-0.999};	
	std::vector<double> point2 = {-0.999,0};	
	std::vector<double> point3 = {-0.999,0.999};		
	std::vector<double> point4 = {0,0.999};	
	std::vector<double> point5 = {0.999,0.999};	
	std::vector< std::vector<double> > point_cloud = {point1,point2,point3,point4,point5};	
	
	Gudhi::Topological_inference_with_cubical_complexes::Distance_to_k_th_closest_point<periodic_Euclidean_distance_squared> 
	periodic_1_n_n_brute_force( point_cloud, period_eu , 1 );    
	Gudhi::Topological_inference_with_cubical_complexes::Distance_to_k_th_closest_point<periodic_Euclidean_distance_squared> 
	periodic_2_n_n_brute_force( point_cloud, period_eu , 2 );    
	Gudhi::Topological_inference_with_cubical_complexes::Distance_to_k_th_closest_point<periodic_Euclidean_distance_squared> 
	periodic_3_n_n_brute_force( point_cloud, period_eu , 3 );    
	Gudhi::Topological_inference_with_cubical_complexes::Distance_to_k_th_closest_point<periodic_Euclidean_distance_squared> 
	periodic_4_n_n_brute_force( point_cloud, period_eu , 4 );    
	Gudhi::Topological_inference_with_cubical_complexes::Distance_to_k_th_closest_point<periodic_Euclidean_distance_squared> 
	periodic_5_n_n_brute_force( point_cloud, period_eu , 5 );
	Gudhi::Topological_inference_with_cubical_complexes::Distance_to_k_th_closest_point<periodic_Euclidean_distance_squared> 
	periodic_6_n_n_brute_force( point_cloud, period_eu , 6 );
	
	{
		std::vector<double> test1 = {1,1};	
		BOOST_CHECK( fabs(periodic_1_n_n_brute_force( test1 ) - 2e-06) <= 2e-07 );
		BOOST_CHECK( fabs(periodic_2_n_n_brute_force( test1 ) - 2e-06) <= 2e-07 );
		BOOST_CHECK( fabs(periodic_3_n_n_brute_force( test1 ) - 2e-06) <= 2e-07 );
		BOOST_CHECK( fabs(periodic_4_n_n_brute_force( test1 ) - 1) <= 2e-05 );
		BOOST_CHECK( fabs(periodic_5_n_n_brute_force( test1 ) - 1) <= 2e-05 );
		BOOST_CHECK( periodic_6_n_n_brute_force( test1 ) == std::numeric_limits<double>::infinity() );
	}
 
	
	{
		std::vector<double> test1 = {0,0};	
		BOOST_CHECK( fabs(periodic_1_n_n_brute_force( test1 ) - 0.998001) <= 5e-06 );
		BOOST_CHECK( fabs(periodic_2_n_n_brute_force( test1 ) - 0.998001) <= 5e-06 );
		BOOST_CHECK( fabs(periodic_3_n_n_brute_force( test1 ) - 1.996) <= 5e-06 );
		BOOST_CHECK( fabs(periodic_4_n_n_brute_force( test1 ) - 1.996) <= 5e-06 );
		BOOST_CHECK( fabs(periodic_5_n_n_brute_force( test1 ) - 1.996) <= 5e-06 );
		BOOST_CHECK( periodic_6_n_n_brute_force( test1 ) == std::numeric_limits<double>::infinity() );		
	}
	
	{
		std::vector<double> test1 = {1,0};	
		BOOST_CHECK( fabs(periodic_1_n_n_brute_force( test1 ) - 1e-06) <= 5e-06 );
		BOOST_CHECK( fabs(periodic_2_n_n_brute_force( test1 ) - 0.998002) <= 5e-06 );
		BOOST_CHECK( fabs(periodic_3_n_n_brute_force( test1 ) - 0.998002) <= 5e-06 );
		BOOST_CHECK( fabs(periodic_4_n_n_brute_force( test1 ) - 0.998002) <= 5e-06 );
		BOOST_CHECK( fabs(periodic_5_n_n_brute_force( test1 ) - 1.998) <= 5e-06 );
		BOOST_CHECK( periodic_6_n_n_brute_force( test1 ) == std::numeric_limits<double>::infinity() );	
	}
	
	{
		std::vector<double> test1 = {0,-1};	
		BOOST_CHECK( fabs(periodic_1_n_n_brute_force( test1 ) - 1e-06) <= 5e-06 );
		BOOST_CHECK( fabs(periodic_2_n_n_brute_force( test1 ) - 0.998002) <= 5e-06 );
		BOOST_CHECK( fabs(periodic_3_n_n_brute_force( test1 ) - 0.998002) <= 5e-06 );
		BOOST_CHECK( fabs(periodic_4_n_n_brute_force( test1 ) - 0.998002) <= 5e-06 );
		BOOST_CHECK( fabs(periodic_5_n_n_brute_force( test1 ) - 1.998) <= 5e-06 );
		BOOST_CHECK( periodic_6_n_n_brute_force( test1 ) == std::numeric_limits<double>::infinity() );	
	}
}	
	
	
BOOST_AUTO_TEST_CASE(Distance_to_k_th_nearest_neighbor_periodic_domain_2d_kd_trees)	
{			
	std::vector< std::pair< double , double > > coordinates_of_grid(2);
	coordinates_of_grid[0] = std::pair< double,double >(-1,1);
	coordinates_of_grid[1] = std::pair< double,double >(-1,1);
	
	std::vector<double> point1 = {-0.999,-0.999};	
	std::vector<double> point2 = {-0.999,0};	
	std::vector<double> point3 = {-0.999,0.999};		
	std::vector<double> point4 = {0,0.999};	
	std::vector<double> point5 = {0.999,0.999};	
	std::vector< std::vector<double> > point_cloud = {point1,point2,point3,point4,point5};	
	
	Gudhi::Topological_inference_with_cubical_complexes::Distance_to_k_th_closest_point_periodic_k_d_tree
	periodic_1_n_n_kd_tree( point_cloud, coordinates_of_grid , 1 );    
	Gudhi::Topological_inference_with_cubical_complexes::Distance_to_k_th_closest_point_periodic_k_d_tree
	periodic_2_n_n_kd_tree( point_cloud, coordinates_of_grid , 2 );    
	Gudhi::Topological_inference_with_cubical_complexes::Distance_to_k_th_closest_point_periodic_k_d_tree
	periodic_3_n_n_kd_tree( point_cloud, coordinates_of_grid , 3 );    
	Gudhi::Topological_inference_with_cubical_complexes::Distance_to_k_th_closest_point_periodic_k_d_tree
	periodic_4_n_n_kd_tree( point_cloud, coordinates_of_grid , 4 );    
	Gudhi::Topological_inference_with_cubical_complexes::Distance_to_k_th_closest_point_periodic_k_d_tree
	periodic_5_n_n_kd_tree( point_cloud, coordinates_of_grid , 5 );
	Gudhi::Topological_inference_with_cubical_complexes::Distance_to_k_th_closest_point_periodic_k_d_tree
	periodic_6_n_n_kd_tree( point_cloud, coordinates_of_grid , 6 );
	
	{
		std::vector<double> test1 = {1,1};	
		BOOST_CHECK( fabs(periodic_1_n_n_kd_tree( test1 ) - 2e-06) <= 2e-07 );
		BOOST_CHECK( fabs(periodic_2_n_n_kd_tree( test1 ) - 2e-06) <= 2e-07 );
		BOOST_CHECK( fabs(periodic_3_n_n_kd_tree( test1 ) - 2e-06) <= 2e-07 );
		BOOST_CHECK( fabs(periodic_4_n_n_kd_tree( test1 ) - 1) <= 2e-05 );
		BOOST_CHECK( fabs(periodic_5_n_n_kd_tree( test1 ) - 1) <= 2e-05 );
		BOOST_CHECK( periodic_6_n_n_kd_tree( test1 ) == std::numeric_limits<double>::infinity() );
	}
 
	
	{
		std::vector<double> test1 = {0,0};	
		BOOST_CHECK( fabs(periodic_1_n_n_kd_tree( test1 ) - 0.998001) <= 5e-06 );
		BOOST_CHECK( fabs(periodic_2_n_n_kd_tree( test1 ) - 0.998001) <= 5e-06 );
		BOOST_CHECK( fabs(periodic_3_n_n_kd_tree( test1 ) - 1.996) <= 5e-06 );
		BOOST_CHECK( fabs(periodic_4_n_n_kd_tree( test1 ) - 1.996) <= 5e-06 );
		BOOST_CHECK( fabs(periodic_5_n_n_kd_tree( test1 ) - 1.996) <= 5e-06 );
		BOOST_CHECK( periodic_6_n_n_kd_tree( test1 ) == std::numeric_limits<double>::infinity() );		
	}
	
	{
		std::vector<double> test1 = {1,0};	
		BOOST_CHECK( fabs(periodic_1_n_n_kd_tree( test1 ) - 1e-06) <= 5e-06 );
		BOOST_CHECK( fabs(periodic_2_n_n_kd_tree( test1 ) - 0.998002) <= 5e-06 );
		BOOST_CHECK( fabs(periodic_3_n_n_kd_tree( test1 ) - 0.998002) <= 5e-06 );
		BOOST_CHECK( fabs(periodic_4_n_n_kd_tree( test1 ) - 0.998002) <= 5e-06 );
		BOOST_CHECK( fabs(periodic_5_n_n_kd_tree( test1 ) - 1.998) <= 5e-06 );
		BOOST_CHECK( periodic_6_n_n_kd_tree( test1 ) == std::numeric_limits<double>::infinity() );	
	}
	
	{
		std::vector<double> test1 = {0,-1};	
		BOOST_CHECK( fabs(periodic_1_n_n_kd_tree( test1 ) - 1e-06) <= 5e-06 );
		BOOST_CHECK( fabs(periodic_2_n_n_kd_tree( test1 ) - 0.998002) <= 5e-06 );
		BOOST_CHECK( fabs(periodic_3_n_n_kd_tree( test1 ) - 0.998002) <= 5e-06 );
		BOOST_CHECK( fabs(periodic_4_n_n_kd_tree( test1 ) - 0.998002) <= 5e-06 );
		BOOST_CHECK( fabs(periodic_5_n_n_kd_tree( test1 ) - 1.998) <= 5e-06 );
		BOOST_CHECK( periodic_6_n_n_kd_tree( test1 ) == std::numeric_limits<double>::infinity() );	
	}
}	


/*
BOOST_AUTO_TEST_CASE(Distance_to_k_th_nearest_neighbor_periodic_domain_3d_brute_force)
{
	//first we test the brute force algorithm:

	typedef Gudhi::Topological_inference_with_cubical_complexes::Euclidan_distance_squared Euclidean_distance_squared;		
	typedef Gudhi::Topological_inference_with_cubical_complexes::periodic_domain_distance<Euclidean_distance_squared>
	periodic_Euclidean_distance_squared;
	
	
	std::vector< std::pair< double , double > > coordinates_of_grid(3);
	coordinates_of_grid[0] = std::pair< double,double >(0,1);
	coordinates_of_grid[1] = std::pair< double,double >(0,1);
	coordinates_of_grid[2] = std::pair< double,double >(0,1);
	
	Euclidean_distance_squared eu;
	periodic_Euclidean_distance_squared period_eu( coordinates_of_grid , eu );
	
	
	//100 random points from R^3.
	std::vector< std::vector<double> > point_cloud = 
	{
	{0.5972671069,0.5497298294,0.4969925818},{0.8563019917,0.6795220347,0.2069536853},{0.0856518536,0.8227915512,0.2959526246},
	{0.0867100377,0.8575005436,0.217976413},{0.9527630396,0.8458755189,0.2894110898},{0.2544106271,0.4017220777,0.7145846318},
	{0.4149809659,0.5972550039,0.7904258594},{0.2131361254,0.3645701057,0.4913443741},{0.31539013,0.202036649,0.3210915248},
	{0.7574333525,0.0939333802,0.4100728314},{0.7516370339,0.5429121791,0.5922422118},{0.2195665021,0.1388090414,0.7959989444},
	{0.899406763,0.3000081666,0.5997520362},{0.5769034552,0.2301354171,0.6781405173},{0.3654143855,0.921912255,0.7039054288},
	{0.5972661381,0.1523829449,0.9388168044},{0.0848337784,0.6187418809,0.7124021121},{0.8300760316,0.7855182383,0.9261578321},
	{0.1363582865,0.6787659361,0.4440253698},{0.9333335564,0.2108276933,0.3914245311},{0.1965385748,0.5836126527,0.0487841049},
	{0.5397062108,0.5839183424,0.6932072763},{0.139552644,0.2591635338,0.3764203726},{0.6378182741,0.1488169709,0.3233298115},
	{0.7192198473,0.4706952481,0.3648997273},{0.947802094,0.4032824829,0.2739019038},{0.0970854505,0.9110400244,0.3080223324},
	{0.1985638938,0.3845967632,0.2895803584},{0.1684866478,0.5466994711,0.3734291738},{0.4505842216,0.6338260781,0.4444476345},
	{0.1916122308,0.5137707235,0.9476647123},{0.3408366609,0.0379224224,0.7878450612},{0.3489282646,0.9416533255,0.3851815299},
	{0.8091952498,0.1183833291,0.3632246931},{0.8706514349,0.1368547638,0.1946790947},{0.4391879758,0.0198513791,0.9993616869},
	{0.8185067559,0.8835033791,0.4361149771},{0.3187520679,0.0242436903,0.692015507},{0.1223861212,0.3934785489,0.051226496},
	{0.2450209786,0.1274023014,0.8058765808},{0.7720149281,0.7432372535,0.2363865508},{0.3001267086,0.7566496953,0.1571175677},
	{0.7458658114,0.0800891486,0.5941177257},{0.0080535775,0.4755061602,0.6350037972},{0.0827106193,0.6615890143,0.028390774},
	{0.7883587447,0.8537377736,0.1755873731},{0.9213447317,0.2877281052,0.0517333397},{0.858527428,0.6600245954,0.9033006672},
	{0.426437065,0.6410801688,0.1879744353},{0.7098625924,0.6176833413,0.6763019743},{0.6594055216,0.3670684479,0.0904959645},
	{0.2070961983,0.5454259284,0.2464962576},{0.109380943,0.4127504176,0.9255148144},{0.8236308896,0.2717485328,0.0889611356},
	{0.6204018274,0.2531314562,0.0879562788},{0.933535517,0.8034739352,0.1976074078},{0.1333329342,0.9984488715,0.4651318933},
	{0.2929060075,0.1734483752,0.3304254888},{0.4427554065,0.7031817061,0.7188385194},{0.7861618018,0.4392747791,0.4609192256},
	{0.4297306708,0.0910678173,0.0856644441},{0.3263101599,0.1900679297,0.8656203696},{0.9598891845,0.8739098033,0.5946957348},
	{0.8601532679,0.3356839109,0.445534257},{0.5811327847,0.0941536282,0.8464888963},{0.3362390934,0.3802505163,0.7506839186},
	{0.5980010226,0.8490654575,0.9566466322},{0.6723008321,0.4367131393,0.853445705},{0.1793271508,0.1859787316,0.1367681916},
	{0.4056227694,0.8398745777,0.9611199587},{0.6466252799,0.3384586412,0.1059214673},{0.1947013228,0.7682137974,0.615737533},
	{0.9453567625,0.7602070384,0.893087266},{0.3547970066,0.9525313561,0.9873204262},{0.0786859554,0.6566422957,0.5613909846},
	{0.6394291024,0.7515391838,0.9339341335},{0.0306413302,0.2682827814,0.3345994204},{0.3132445228,0.0521786655,0.1866439714},
	{0.3550304624,0.814520665,0.0284737975},{0.0295199333,0.3506491578,0.0804069594},{0.0772441155,0.4198115987,0.5326562512},
	{0.6100857158,0.0040999444,0.0831175786},{0.2230040201,0.0384767053,0.9185578863},{0.4115153016,0.9873048151,0.7218966454},
	{0.0724200262,0.6555718367,0.0178539874},{0.581852132,0.2605712679,0.1275672165},{0.2533805799,0.2765149043,0.3912544115},
	{0.3027908381,0.9986928527,0.4299810582},{0.3472439814,0.7639185048,0.4064361725},{0.9034222241,0.4131227473,0.6427890444},
	{0.7460419929,0.6069772805,0.7167975691},{0.3979357823,0.1354945437,0.9507125753},{0.0570875325,0.028278504,0.9812576636},
	{0.1229622571,0.0364863402,0.7952201411},{0.5176588467,0.7044958638,0.110351444},{0.5363182717,0.3783060559,0.0327049491},
	{0.980905558,0.2990631682,0.435821288},{0.6473647165,0.2050705738,0.3516822953},{0.9423140211,0.3945645185,0.4104162175},
	{0.1950914182,0.8122796561,0.3180316554}
	};	
	
	std::vector< Gudhi::Topological_inference_with_cubical_complexes::Distance_to_k_th_closest_point<periodic_Euclidean_distance_squared>* >
	classes_to_compute_knn;
	classes_to_compute_knn.reserve(91);
	for ( size_t i = 0 ; i != 90 ; i=i+5 )
	{
		classes_to_compute_knn.push_back(
		new Gudhi::Topological_inference_with_cubical_complexes::Distance_to_k_th_closest_point<periodic_Euclidean_distance_squared>
		( point_cloud, period_eu , i+1 ) );
	}  
	
	
	
	//100 random points to test the distance:
	std::vector< std::vector<double> > test_points =
	{
	{0.0278817846,0.2916815139,0.3887573609},{0.3553810648,0.9343284944,0.5285676813},{0.6866485006,0.8520970619,0.8668888684},
	{0.8956627958,0.1156968276,0.0739414582},{0.6355211725,0.3484526002,0.743152576},{0.0726810361,0.7873784064,0.385060332},
	{0.7665069464,0.2014935308,0.6368220472},{0.6570753611,0.5955838757,0.8398550118},{0.9891614937,0.7330921197,0.7102547155},
	{0.7883918153,0.014290286,0.7783693536},{0.0336034244,0.6878910807,0.2098890906},{0.8529858275,0.2412898415,0.3730617869},
	{0.6427641907,0.200828871,0.252669279},{0.1852649881,0.2579907831,0.5425100487},{0.663123504,0.5577539282,0.5905491377},
	{0.3354653688,0.8227474587,0.9516165522},{0.6034799691,0.6799063885,0.8784505285},{0.719459174,0.5953775912,0.4358448796},
	{0.42901467,0.0499089628,0.338926865},{0.6868318615,0.4042496786,0.219024647},{0.4858837686,0.3048673074,0.2064413528},
	{0.3251100981,0.6204368386,0.1215259079},{0.384246537,0.4619099458,0.577011209},{0.1957502051,0.2889291102,0.4228005044},
	{0.6831766628,0.2632829093,0.5414486954},{0.4916156209,0.6092843614,0.7656351738},{0.4322729874,0.1720699086,0.8144672341},
	{0.0017052968,0.2426031898,0.6799070754},{0.993835246,0.3934641157,0.8028928195},{0.8344388744,0.1533741863,0.052077933},
	{0.3089801904,0.0926956008,0.3230055785},{0.0271170982,0.1870652179,0.3214361172},{0.832586685,0.5277374319,0.2314551908},
	{0.9314843167,0.6746024499,0.2297013158},{0.0021536089,0.0794244623,0.2783408458},{0.2407027101,0.1303185185,0.7768267517},
	{0.8024219871,0.8166264719,0.9408042894},{0.4125754705,0.1070645854,0.3917626143},{0.7881310335,0.6119206273,0.3612281107},
	{0.2093140029,0.3191626819,0.4877916307},{0.8854867723,0.9801039936,0.1372925977},{0.8612136859,0.7057745687,0.8465774728},
	{0.8764370135,0.2966763664,0.8315236054},{0.3119052621,0.5994925359,0.3814946879},{0.7070578339,0.9780796578,0.5635364696},
	{0.5851224931,0.4823992441,0.8954798079},{0.3127624721,0.8012911966,0.4995430515},{0.5117247337,0.9834836028,0.6449389677},
	{0.2507975928,0.6793389737,0.3178833791},{0.0095620689,0.4389932367,0.9797510416},{0.9370902732,0.0666056497,0.6635932433},
	{0.8764814257,0.5038128984,0.8082919149},{0.5889022176,0.2642431161,0.7928514399},{0.9185145744,0.046059252,0.5190000213},
	{0.4452630868,0.7133746753,0.3376530202},{0.9111526795,0.5850539894,0.4890879511},{0.2747024964,0.1768817641,0.6366421564},
	{0.5592753026,0.6720825401,0.5405601789},{0.6614501821,0.5494211896,0.4482777803},{0.331644475,0.4239660688,0.2038797489},
	{0.0581143782,0.0254596837,0.7859259995},{0.1199549222,0.0386884976,0.4347598413},{0.292417774,0.3479369427,0.1531838665},
	{0.1128242332,0.7438441454,0.745681555},{0.2401035188,0.1771410096,0.807042615},{0.1234689017,0.6644280744,0.8244124809},
	{0.6452131516,0.879331308,0.2716083024},{0.9479403694,0.332957929,0.1844342141},{0.8303372548,0.1892285589,0.4289489985},
	{0.7557857418,0.0460203474,0.2242759808},{0.8151567031,0.5541405319,0.8722450782},{0.2119767461,0.9610106237,0.6250581758},
	{0.2892467363,0.073150466,0.657205262},{0.8341697995,0.5837534382,0.5041492598},{0.7957765597,0.0565284048,0.8094913445},
	{0.8098659823,0.7344695956,0.4325275582},{0.4223770299,0.7119579143,0.4925132662},{0.7018576597,0.5692464605,0.441388723},
	{0.878948973,0.4988749905,0.3432697712},{0.7926970301,0.7522489438,0.9827959056},{0.0475694353,0.0330607148,0.4641666713},
	{0.9334698259,0.3049679426,0.8354724362},{0.6555848967,0.1480603861,0.7719476013},{0.1087768308,0.9330630642,0.1765642043},
	{0.9145491763,0.3303631963,0.9045855913},{0.9143750975,0.8296238736,0.09732326},{0.0618927069,0.897331967,0.6368782204},
	{0.3150910612,0.9116620645,0.6570473495},{0.3580378599,0.2623571989,0.814647858},{0.9664658206,0.8377381989,0.7721899245},
	{0.9837020198,0.2720242415,0.4563574172},{0.91944414,0.4706337838,0.6330762212},{0.5168302623,0.1900973888,0.9053026016},
	{0.1698169338,0.8422267565,0.4352734836},{0.1704680121,0.7366362549,0.1204422733},{0.1558771764,0.3680324578,0.922587656},
	{0.377208197,0.3812016125,0.2776982288},{0.6930598826,0.2130688538,0.7388268828},{0.4868098402,0.0608007649,0.5235421255},
	{0.103365178,0.8140293758,0.7529569559},{0.5972671069,0.5497298294,0.4969925818},{0.8563019917,0.6795220347,0.2069536853},
	{0.0856518536,0.8227915512,0.2959526246},{0.0867100377,0.8575005436,0.217976413},{0.9527630396,0.8458755189,0.2894110898},
	{0.2544106271,0.4017220777,0.7145846318},{0.4149809659,0.5972550039,0.7904258594},{0.2131361254,0.3645701057,0.4913443741},
	{0.31539013,0.202036649,0.3210915248},{0.7574333525,0.0939333802,0.4100728314},{0.7516370339,0.5429121791,0.5922422118},
	{0.2195665021,0.1388090414,0.7959989444},{0.899406763,0.3000081666,0.5997520362},{0.5769034552,0.2301354171,0.6781405173},
	{0.3654143855,0.921912255,0.7039054288},{0.5972661381,0.1523829449,0.9388168044},{0.0848337784,0.6187418809,0.7124021121},
	{0.8300760316,0.7855182383,0.9261578321},{0.1363582865,0.6787659361,0.4440253698},{0.9333335564,0.2108276933,0.3914245311},
	{0.1965385748,0.5836126527,0.0487841049},{0.5397062108,0.5839183424,0.6932072763},{0.139552644,0.2591635338,0.3764203726},
	{0.6378182741,0.1488169709,0.3233298115},{0.7192198473,0.4706952481,0.3648997273},{0.947802094,0.4032824829,0.2739019038},
	{0.0970854505,0.9110400244,0.3080223324},{0.1985638938,0.3845967632,0.2895803584},{0.1684866478,0.5466994711,0.3734291738},
	{0.4505842216,0.6338260781,0.4444476345},{0.1916122308,0.5137707235,0.9476647123},{0.3408366609,0.0379224224,0.7878450612},
	{0.3489282646,0.9416533255,0.3851815299},{0.8091952498,0.1183833291,0.3632246931},{0.8706514349,0.1368547638,0.1946790947},
	{0.4391879758,0.0198513791,0.9993616869},{0.8185067559,0.8835033791,0.4361149771},{0.3187520679,0.0242436903,0.692015507},
	{0.1223861212,0.3934785489,0.051226496},{0.2450209786,0.1274023014,0.8058765808},{0.7720149281,0.7432372535,0.2363865508},
	{0.3001267086,0.7566496953,0.1571175677},{0.7458658114,0.0800891486,0.5941177257},{0.0080535775,0.4755061602,0.6350037972},
	{0.0827106193,0.6615890143,0.028390774},{0.7883587447,0.8537377736,0.1755873731},{0.9213447317,0.2877281052,0.0517333397},
	{0.858527428,0.6600245954,0.9033006672},{0.426437065,0.6410801688,0.1879744353},{0.7098625924,0.6176833413,0.6763019743},
	{0.6594055216,0.3670684479,0.0904959645},{0.2070961983,0.5454259284,0.2464962576},{0.109380943,0.4127504176,0.9255148144},
	{0.8236308896,0.2717485328,0.0889611356},{0.6204018274,0.2531314562,0.0879562788},{0.933535517,0.8034739352,0.1976074078},
	{0.1333329342,0.9984488715,0.4651318933},{0.2929060075,0.1734483752,0.3304254888},{0.4427554065,0.7031817061,0.7188385194},
	{0.7861618018,0.4392747791,0.4609192256},{0.4297306708,0.0910678173,0.0856644441},{0.3263101599,0.1900679297,0.8656203696},
	{0.9598891845,0.8739098033,0.5946957348},{0.8601532679,0.3356839109,0.445534257},{0.5811327847,0.0941536282,0.8464888963},
	{0.3362390934,0.3802505163,0.7506839186},{0.5980010226,0.8490654575,0.9566466322},{0.6723008321,0.4367131393,0.853445705},
	{0.1793271508,0.1859787316,0.1367681916},{0.4056227694,0.8398745777,0.9611199587},{0.6466252799,0.3384586412,0.1059214673},
	{0.1947013228,0.7682137974,0.615737533},{0.9453567625,0.7602070384,0.893087266},{0.3547970066,0.9525313561,0.9873204262},
	{0.0786859554,0.6566422957,0.5613909846},{0.6394291024,0.7515391838,0.9339341335},{0.0306413302,0.2682827814,0.3345994204},
	{0.3132445228,0.0521786655,0.1866439714},{0.3550304624,0.814520665,0.0284737975},{0.0295199333,0.3506491578,0.0804069594},
	{0.0772441155,0.4198115987,0.5326562512},{0.6100857158,0.0040999444,0.0831175786},{0.2230040201,0.0384767053,0.9185578863},
	{0.4115153016,0.9873048151,0.7218966454},{0.0724200262,0.6555718367,0.0178539874},{0.581852132,0.2605712679,0.1275672165},
	{0.2533805799,0.2765149043,0.3912544115},{0.3027908381,0.9986928527,0.4299810582},{0.3472439814,0.7639185048,0.4064361725},
	{0.9034222241,0.4131227473,0.6427890444},{0.7460419929,0.6069772805,0.7167975691},{0.3979357823,0.1354945437,0.9507125753},
	{0.0570875325,0.028278504,0.9812576636},{0.1229622571,0.0364863402,0.7952201411},{0.5176588467,0.7044958638,0.110351444},
	{0.5363182717,0.3783060559,0.0327049491},{0.980905558,0.2990631682,0.435821288},{0.6473647165,0.2050705738,0.3516822953},
	{0.9423140211,0.3945645185,0.4104162175},{0.1950914182,0.8122796561,0.3180316554}
	};	
	
	
	std::vector<double> result = 
	{
		0.0034882,0.0320593,0.051086,0.0863591,0.0976073,0.124951,0.156195,0.169588,0.212912,0.225603,
		0.231808,0.248094,0.268131,0.286455,0.294607,0.323273,0.350672,0.383048,0.0166278,0.0440219,
		0.0972661,0.120431,0.139304,0.167899,0.178021,0.189844,0.207047,0.225616,0.24783,0.254605,
		0.273384,0.295205,0.31428,0.346701,0.387461,0.413231,0.015924,0.0757216,0.103334,0.115002,
		0.137492,0.170696,0.17525,0.210232,0.229506,0.244188,0.27455,0.283341,0.298765,0.332275,
		0.354117,0.366531,0.411705,0.435295,0.0156508,0.0893522,0.108645,0.119279,0.129082,0.137249,
		0.161551,0.180403,0.193159,0.210556,0.226921,0.243276,0.254512,0.271362,0.303793,0.328946,
		0.336952,0.385805,0.0213072,0.0783045,0.0925462,0.112775,0.139199,0.163583,0.180605,0.204864,
		0.216111,0.231345,0.261933,0.275493,0.297022,0.324413,0.336924,0.352445,0.377303,0.389444,
		0.00936252,0.033031,0.0672403,0.0958084,0.129122,0.163913,0.175481,0.198768,0.207612,0.225629,
		0.23514,0.249451,0.268238,0.28034,0.295024,0.313949,0.352132,0.38366,0.0169887,0.0635683,
		0.0955123,0.119444,0.148647,0.177981,0.205266,0.227535,0.234948,0.245536,0.260792,0.281141,
		0.292393,0.303792,0.318963,0.340777,0.366175,0.396293,0.023188,0.0487609,0.081383,0.123234,
		0.166132,0.181498,0.199758,0.209506,0.227659,0.245042,0.254383,0.271964,0.301809,0.312184,
		0.367229,0.390766,0.413438,0.443326,0.0222338,0.0596709,0.0954963,0.11717,0.163408,0.180182,
		0.191332,0.206201,0.221479,0.230088,0.244653,0.249348,0.267331,0.284904,0.301267,0.346635,
		0.380036,0.387296,0.0400866,0.0953307,0.115428,0.145953,0.171488,0.190681,0.200954,0.217353,
		0.228118,0.238348,0.25203,0.260735,0.281391,0.294462,0.310757,0.331518,0.357426,0.385487,
		0.0235238,0.037818,0.0634556,0.0888267,0.125425,0.134585,0.162313,0.186054,0.197746,0.2175,
		0.228353,0.244493,0.259194,0.281214,0.317533,0.327824,0.374619,0.395035,0.00772089,0.032868,
		0.0513824,0.0825026,0.12339,0.135521,0.158897,0.173305,0.199316,0.225973,0.240542,0.264418,
		0.293049,0.311375,0.335888,0.359915,0.38462,0.396642,0.00772261,0.0467191,0.0685179,0.10292,
		0.129196,0.147612,0.165233,0.186222,0.214987,0.232447,0.261199,0.275102,0.301867,0.314187,
		0.324426,0.351452,0.359149,0.38721,0.0147539,0.0550495,0.0796375,0.0885206,0.112223,0.132822,
		0.169847,0.176402,0.193062,0.227384,0.245333,0.256929,0.271687,0.288861,0.319903,0.337957,
		0.378005,0.398558,0.00805779,0.0459796,0.0861695,0.127718,0.164004,0.188042,0.199291,0.22503,
		0.243783,0.264347,0.283221,0.287153,0.309275,0.335562,0.351789,0.358103,0.376937,0.408084,
		0.00530571,0.0602795,0.0799927,0.0971439,0.108273,0.122235,0.142403,0.165529,0.191646,0.236178,
		0.246674,0.261101,0.277991,0.291826,0.324938,0.349762,0.397726,0.430102,0.00950203,0.0518501,
		0.066062,0.123542,0.14802,0.162495,0.177536,0.216675,0.234396,0.243916,0.261265,0.283349,
		0.310663,0.326889,0.35146,0.377535,0.406929,0.420193,0.020579,0.0644082,0.0906368,0.13732,
		0.148559,0.164621,0.180162,0.207891,0.218741,0.243182,0.259767,0.282366,0.288796,0.300069,
		0.308741,0.334616,0.377559,0.43654,0.0202726,0.053625,0.100324,0.116311,0.137488,0.153638,
		0.168822,0.184725,0.194291,0.213131,0.22341,0.261955,0.270897,0.294042,0.302052,0.316873,
		0.354914,0.379517,0.0186543,0.053187,0.0785278,0.105882,0.127897,0.144544,0.181969,0.204134,
		0.226566,0.24926,0.252311,0.26457,0.277559,0.294151,0.300744,0.334466,0.355782,0.4181,
		0.0173932,0.0527869,0.0890175,0.10728,0.128972,0.159125,0.169946,0.182629,0.197202,0.226249,
		0.241046,0.266963,0.276132,0.286615,0.304884,0.321501,0.34116,0.372843,0.0151087,0.0472225,
		0.092323,0.111125,0.128712,0.16185,0.183127,0.195456,0.21498,0.223237,0.229884,0.236976,
		0.256125,0.270716,0.285884,0.320539,0.34758,0.360937,0.0391352,0.0594931,0.0979901,0.131249,
		0.141218,0.157276,0.177182,0.199423,0.226101,0.237899,0.259706,0.272615,0.290271,0.298388,
		0.321404,0.334539,0.350136,0.381014,0.00447053,0.0322087,0.0755421,0.0967258,0.119254,0.149275,
		0.165709,0.184363,0.212991,0.22495,0.234737,0.245496,0.26694,0.276325,0.290454,0.307212,
		0.320255,0.368186,0.0310774,0.0514529,0.08123,0.101557,0.152825,0.173642,0.182525,0.207458,
		0.227126,0.23553,0.265189,0.28278,0.309248,0.315114,0.33359,0.351497,0.393254,0.40791,
		0.00663215,0.0701387,0.101954,0.126676,0.151229,0.171381,0.206958,0.212671,0.225942,0.244903,
		0.261883,0.27071,0.294639,0.315719,0.319282,0.359004,0.379506,0.414451,0.0141687,0.0428745,
		0.05663,0.0840779,0.128238,0.155355,0.174518,0.192419,0.229111,0.253721,0.266839,0.284006,
		0.304243,0.320144,0.340307,0.359512,0.376748,0.407586,0.0201852,0.0704845,0.0903813,0.11138,
		0.139814,0.147811,0.160405,0.180615,0.199885,0.221949,0.234697,0.243478,0.26641,0.291355,
		0.303552,0.321732,0.359467,0.376999,0.0287589,0.0745483,0.080679,0.120141,0.133892,0.144986,
		0.156694,0.195256,0.215812,0.230496,0.240407,0.254537,0.270993,0.289579,0.327544,0.359125,
		0.370231,0.390361,0.0154897,0.0702369,0.0809899,0.124511,0.137623,0.153429,0.168799,0.189937,
		0.206079,0.214112,0.231868,0.257603,0.283953,0.305237,0.315463,0.330377,0.370045,0.386091,
		0.00683444,0.0415388,0.0781226,0.11144,0.127955,0.148347,0.17391,0.188939,0.20466,0.21753,
		0.2284,0.235986,0.245361,0.261888,0.287317,0.308637,0.326301,0.344025,0.00678199,0.0539534,
		0.0653648,0.0833256,0.101298,0.122965,0.141833,0.164901,0.171233,0.214141,0.226901,0.260293,
		0.2735,0.284798,0.297294,0.331034,0.34253,0.372375,0.0242013,0.0626347,0.087367,0.117181,
		0.126369,0.140796,0.152876,0.160273,0.175916,0.185481,0.191968,0.22104,0.245994,0.281904,
		0.336629,0.352809,0.374228,0.428884,0.00619404,0.0555027,0.0879255,0.0990101,0.114725,0.147327,
		0.173751,0.185931,0.193731,0.209772,0.229725,0.239801,0.259637,0.279044,0.294502,0.316063,
		0.376551,0.42459,0.0275903,0.0571071,0.0731426,0.0960903,0.107201,0.114709,0.13958,0.164667,
		0.186178,0.20948,0.2322,0.243551,0.253546,0.264719,0.284124,0.320801,0.357015,0.36667,
		0.000871044,0.0245365,0.0722767,0.119121,0.136423,0.158547,0.176633,0.187158,0.211702,0.237728,
		0.256056,0.281387,0.291379,0.302148,0.311798,0.329751,0.357408,0.381351,0.00194699,0.0566981,
		0.0973105,0.12241,0.158403,0.172279,0.198994,0.221283,0.233274,0.245871,0.258182,0.271382,
		0.288197,0.29878,0.310681,0.355217,0.375591,0.411152,0.0224898,0.0549528,0.0979109,0.122233,
		0.133939,0.156826,0.190028,0.198672,0.216053,0.225144,0.238872,0.251635,0.255484,0.264564,
		0.290123,0.30415,0.345118,0.389567,0.0247068,0.0594615,0.0869963,0.126487,0.140606,0.14988,
		0.16998,0.187756,0.205392,0.219529,0.225913,0.248166,0.258742,0.275206,0.302196,0.343465,
		0.390069,0.419116,0.00208906,0.0527596,0.0665186,0.108945,0.123972,0.139225,0.172534,0.186426,
		0.20582,0.218231,0.233843,0.260965,0.271891,0.283246,0.310065,0.333346,0.351295,0.380993,
		0.0268688,0.0620323,0.0855119,0.103239,0.128727,0.152209,0.15793,0.176587,0.199959,0.222092,
		0.23818,0.253278,0.267853,0.287702,0.316401,0.328232,0.347138,0.380521,0.0053178,0.059661,
		0.101928,0.130584,0.153681,0.1684,0.207621,0.214945,0.225387,0.238666,0.264474,0.291668,
		0.305208,0.310701,0.326735,0.343591,0.355204,0.413212,0.0499087,0.0765704,0.118131,0.128444,
		0.137495,0.157458,0.165832,0.194351,0.214467,0.224938,0.243452,0.271519,0.287924,0.295988,
		0.322332,0.357655,0.372783,0.395833,0.0234211,0.0522972,0.09002,0.110201,0.144019,0.161432,
		0.177854,0.184483,0.198049,0.218587,0.243281,0.252218,0.268991,0.281849,0.301301,0.3307,
		0.358178,0.398729,0.0128472,0.0916446,0.134995,0.161394,0.169543,0.183708,0.193089,0.200184,
		0.230065,0.240883,0.263443,0.281163,0.28817,0.301072,0.317521,0.325939,0.343139,0.380137,
		0.0114541,0.0687872,0.0933432,0.110936,0.138695,0.158432,0.175357,0.213457,0.232083,0.246766,
		0.264601,0.280999,0.29711,0.310755,0.322257,0.370587,0.398234,0.422456,0.0112546,0.0492131,
		0.0795403,0.101537,0.13626,0.144233,0.192449,0.209303,0.223868,0.24109,0.264076,0.278162,
		0.291552,0.297395,0.310602,0.32488,0.356655,0.389306,0.015979,0.0661877,0.117744,0.129553,
		0.142809,0.153625,0.175735,0.19119,0.212536,0.236436,0.257973,0.286614,0.294478,0.314625,
		0.333331,0.35671,0.378806,0.402957,0.0207764,0.034256,0.07741,0.0997626,0.11737,0.143603,
		0.16892,0.186861,0.201038,0.232,0.245021,0.25699,0.279849,0.286997,0.303664,0.326128,
		0.346817,0.371127,0.0135941,0.0523092,0.0916143,0.121501,0.134675,0.157749,0.175319,0.191969,
		0.209716,0.227174,0.233643,0.24995,0.259296,0.275169,0.289119,0.320225,0.329517,0.347061,
		0.0415754,0.0948899,0.109259,0.146616,0.157691,0.171683,0.192924,0.200402,0.215595,0.22894,
		0.243137,0.253325,0.265278,0.28684,0.299974,0.315167,0.353928,0.369181,0.0337511,0.0581497,
		0.0855508,0.118831,0.132991,0.160127,0.176976,0.191561,0.222033,0.233181,0.251086,0.265119,
		0.286042,0.304068,0.348774,0.358002,0.383313,0.428396,0.0144659,0.0779645,0.104136,0.114541,
		0.140749,0.144383,0.169105,0.199284,0.219504,0.242748,0.259915,0.277109,0.28814,0.305113,
		0.324523,0.354784,0.38559,0.432446,0.0366083,0.0436437,0.0939576,0.118188,0.150287,0.165072,
		0.175963,0.189016,0.205899,0.233201,0.244872,0.251785,0.27675,0.290751,0.317549,0.336147,
		0.360345,0.406905,0.0168935,0.0636505,0.107935,0.143031,0.159351,0.171266,0.187181,0.194819,
		0.208406,0.226893,0.238217,0.256271,0.275629,0.293452,0.322614,0.334431,0.362544,0.374904,
		0.0376678,0.0532441,0.0766315,0.0894936,0.100461,0.133169,0.161513,0.179638,0.211812,0.232852,
		0.253041,0.273043,0.283411,0.29249,0.328758,0.342371,0.386377,0.431053,0.028305,0.05704,
		0.0705965,0.0941118,0.11554,0.15915,0.177817,0.193697,0.217328,0.230805,0.23703,0.25939,
		0.276671,0.29129,0.312627,0.338401,0.357766,0.390151,0.0183117,0.0563591,0.11202,0.142924,
		0.179063,0.196322,0.204156,0.216798,0.228634,0.242704,0.261535,0.280634,0.296057,0.320058,
		0.32859,0.361758,0.381249,0.392275,0.00649269,0.0589985,0.104299,0.133761,0.158918,0.164886,
		0.191798,0.201409,0.219304,0.242286,0.253821,0.258416,0.29107,0.310315,0.329077,0.342713,
		0.387134,0.415628,0.026605,0.0677953,0.0843425,0.1119,0.123521,0.138862,0.15684,0.172831,
		0.192415,0.214571,0.239013,0.248144,0.269685,0.288495,0.296917,0.308316,0.328874,0.358868,
		0.00441322,0.0691846,0.109296,0.133899,0.163549,0.186865,0.191484,0.203393,0.21463,0.227825,
		0.236363,0.251866,0.263414,0.281907,0.285074,0.299578,0.337665,0.369689,0.00272066,0.0643026,
		0.0762565,0.0871308,0.114018,0.129899,0.146326,0.158797,0.173459,0.220572,0.24468,0.251139,
		0.262776,0.273323,0.308727,0.345468,0.375094,0.405731,0.0287565,0.0618611,0.0799007,0.0920612,
		0.108764,0.138976,0.161352,0.177127,0.187747,0.196707,0.222083,0.237173,0.244761,0.272471,
		0.307613,0.329287,0.338284,0.358068,0.0175416,0.0835023,0.0965364,0.114789,0.139506,0.143952,
		0.163813,0.18214,0.201959,0.216928,0.239937,0.263473,0.28056,0.300003,0.329686,0.352192,
		0.365791,0.389427,0.00201307,0.0336456,0.0726696,0.101361,0.132858,0.160234,0.17371,0.193733,
		0.21556,0.229314,0.256325,0.268188,0.28838,0.300186,0.329723,0.354538,0.384338,0.390922,
		0.0161262,0.0593908,0.0848878,0.11779,0.13928,0.150405,0.172817,0.195192,0.208889,0.233478,
		0.24075,0.256611,0.277036,0.296597,0.304819,0.320957,0.365195,0.412102,0.0303657,0.0753523,
		0.0960232,0.120284,0.13631,0.155363,0.170519,0.198272,0.214085,0.22974,0.249355,0.276129,
		0.291382,0.30407,0.318188,0.333733,0.351823,0.40639,0.01295,0.044535,0.0758872,0.0969853,
		0.116156,0.139254,0.146023,0.165522,0.194511,0.215908,0.226947,0.248923,0.267432,0.273816,
		0.300388,0.323227,0.341152,0.379683,0.00978572,0.0397001,0.0552731,0.0956729,0.127998,0.157666,
		0.164432,0.188005,0.212975,0.243114,0.259517,0.279207,0.294964,0.306579,0.333523,0.366232,
		0.410199,0.425374,0.022321,0.0429115,0.0831006,0.111444,0.135313,0.144725,0.162255,0.186322,
		0.201461,0.210899,0.227312,0.252535,0.27967,0.297126,0.328001,0.345807,0.381969,0.410253,
		0.0140569,0.0598499,0.099676,0.109397,0.130692,0.152237,0.176585,0.18443,0.202055,0.254919,
		0.267019,0.285299,0.302074,0.333335,0.349915,0.371294,0.379864,0.408537,0.0198827,0.0477222,
		0.0720566,0.116209,0.138488,0.155921,0.180741,0.19519,0.203484,0.210774,0.237189,0.268502,
		0.277398,0.292938,0.306758,0.323862,0.36839,0.387796,0.00447419,0.0308554,0.0738956,0.107829,
		0.121467,0.152081,0.180459,0.203634,0.216269,0.234634,0.241453,0.26465,0.285705,0.300148,
		0.31892,0.341985,0.365906,0.398475,0.01624,0.0535252,0.0683776,0.0984959,0.108195,0.153209,
		0.16983,0.201543,0.217128,0.237053,0.254663,0.269499,0.282532,0.299911,0.33092,0.37309,
		0.404173,0.43114,0.0488564,0.0985851,0.114012,0.146952,0.164551,0.184456,0.194745,0.213162,
		0.217706,0.235901,0.25105,0.270979,0.284637,0.291793,0.302454,0.346793,0.375057,0.390408,
		0.0222986,0.0682491,0.0834828,0.109832,0.139355,0.156226,0.167562,0.213208,0.222181,0.2347,
		0.254408,0.284482,0.302912,0.319533,0.351087,0.361539,0.391979,0.408131,0.00921053,0.0701852,
		0.097784,0.125928,0.148062,0.166083,0.184102,0.206078,0.230387,0.250565,0.265062,0.281009,
		0.29849,0.313613,0.334311,0.338279,0.361671,0.375241,0.0144119,0.0673182,0.0899209,0.136594,
		0.159068,0.16487,0.182903,0.200792,0.222569,0.231478,0.266405,0.281115,0.290814,0.302909,
		0.32406,0.337756,0.391653,0.429248,0.0186906,0.0517283,0.0825718,0.105562,0.116978,0.128759,
		0.149887,0.155659,0.168873,0.191946,0.210417,0.228383,0.249106,0.275591,0.288285,0.362984,
		0.412682,0.436979,0.00571191,0.0479638,0.0940357,0.1284,0.176201,0.188866,0.197516,0.207719,
		0.223639,0.23189,0.243924,0.271217,0.289136,0.303743,0.333183,0.361387,0.376637,0.394107,
		0.00855429,0.0674878,0.075624,0.0977854,0.114697,0.127053,0.147934,0.164653,0.183218,0.202404,
		0.233354,0.2479,0.268869,0.309191,0.329127,0.337429,0.367731,0.381953,0.047213,0.0748337,
		0.111018,0.12947,0.144018,0.158371,0.165383,0.175493,0.20019,0.210925,0.249843,0.25773,
		0.284244,0.301632,0.331122,0.342337,0.354624,0.392284,0.0140055,0.0902418,0.114982,0.135232,
		0.14787,0.168415,0.18373,0.197428,0.2096,0.225882,0.241616,0.263715,0.276017,0.308906,
		0.33029,0.376238,0.393141,0.427909,0.00791161,0.0479457,0.0881496,0.0985629,0.115362,0.136963,
		0.147231,0.156321,0.176501,0.19374,0.210058,0.224717,0.256494,0.289892,0.296262,0.327672,
		0.357572,0.397703,0.0235164,0.0726099,0.111817,0.126112,0.142818,0.156698,0.169276,0.180537,
		0.203997,0.225185,0.249478,0.266385,0.284886,0.310056,0.325127,0.339837,0.351452,0.384279,
		0.0111079,0.0450336,0.068836,0.108389,0.142494,0.177567,0.199021,0.214064,0.222707,0.238512,
		0.259981,0.27228,0.28185,0.290112,0.313104,0.337933,0.349577,0.365672,0.0127327,0.0838426,
		0.0997338,0.120366,0.13529,0.152328,0.182788,0.203683,0.208371,0.223864,0.240778,0.260535,
		0.274544,0.289484,0.312824,0.321428,0.347675,0.374247,0.00483318,0.0592848,0.0773995,0.112327,
		0.131564,0.143203,0.170416,0.197184,0.217895,0.236236,0.247875,0.271763,0.28499,0.305642,
		0.314344,0.332412,0.391748,0.418216,0.00883058,0.0401738,0.0791539,0.0995136,0.125818,0.141533,
		0.173914,0.186347,0.210683,0.240719,0.248556,0.269177,0.283602,0.292113,0.311454,0.351564,
		0.364964,0.400202,0.0210728,0.0655449,0.104906,0.139108,0.161828,0.170908,0.194708,0.219751,
		0.233329,0.247804,0.256128,0.263074,0.282671,0.29499,0.308238,0.333534,0.363609,0.380973,
		0.00116066,0.0284508,0.0624288,0.0850561,0.116469,0.145699,0.160603,0.177223,0.189699,0.216427,
		0.236505,0.244487,0.269717,0.291593,0.316127,0.344719,0.350681,0.396414,0.00365856,0.0483856,
		0.0650959,0.112601,0.126111,0.15213,0.174902,0.2033,0.217881,0.226618,0.253411,0.262359,
		0.287579,0.29856,0.315483,0.335854,0.362185,0.407558,0.00901551,0.0480627,0.0679272,0.0876796,
		0.109498,0.138139,0.182408,0.200603,0.215998,0.251806,0.272533,0.286793,0.300313,0.320199,
		0.350453,0.37442,0.39582,0.413671,0.0152813,0.0384441,0.0586521,0.113816,0.125123,0.163692,
		0.178216,0.200544,0.211029,0.220559,0.234191,0.251704,0.268003,0.280431,0.310219,0.339497,
		0.349148,0.375177,0.018557,0.0453698,0.070988,0.0983093,0.113787,0.134625,0.156003,0.181997,
		0.194014,0.222254,0.245834,0.255576,0.269945,0.283934,0.291925,0.334109,0.354522,0.376537,
		0.00417016,0.0622306,0.079471,0.113129,0.13678,0.156967,0.174627,0.197882,0.216401,0.226836,
		0.241592,0.255809,0.267635,0.280454,0.303121,0.331618,0.358913,0.386116,0.0320665,0.0728409,
		0.0853461,0.111547,0.124002,0.152251,0.163944,0.188805,0.206926,0.223053,0.238986,0.246895,
		0.267609,0.28034,0.302495,0.31141,0.324083,0.353644,0.0174664,0.0694793,0.130523,0.143086,
		0.155411,0.16353,0.177953,0.182063,0.206877,0.225436,0.236666,0.253953,0.292534,0.314811,
		0.318756,0.357,0.404251,0.436845,0.0464741,0.0665578,0.0875822,0.121332,0.149101,0.176313,
		0.193041,0.206785,0.230349,0.24424,0.251957,0.261625,0.27509,0.291261,0.307285,0.318123,
		0.346808,0.389856,0.0292705,0.0620773,0.0995285,0.114823,0.128637,0.139245,0.182065,0.199272,
		0.209243,0.219747,0.23857,0.26619,0.288849,0.31548,0.326345,0.341234,0.35349,0.386658,
		0,0.0491828,0.117573,0.135897,0.164058,0.186803,0.200008,0.220325,0.231893,0.240959,
		0.248529,0.25717,0.288885,0.300884,0.336598,0.353499,0.388684,0.412127,0,0.0810483,
		0.0891626,0.121792,0.134635,0.149977,0.17516,0.181309,0.207137,0.214424,0.222769,0.247878,
		0.263141,0.289485,0.307403,0.346491,0.362594,0.377572,0,0.0331843,0.0841031,0.0960542,
		0.107677,0.142074,0.167103,0.190934,0.202982,0.215067,0.234584,0.251784,0.263993,0.296326,
		0.310498,0.321332,0.353359,0.367754,0,0.0267962,0.0848857,0.103798,0.112697,0.12745,
		0.156201,0.167991,0.18923,0.202223,0.233724,0.253629,0.272193,0.291968,0.30636,0.326018,
		0.358289,0.404441,0,0.0400464,0.0855384,0.114246,0.143102,0.164467,0.183069,0.190509,
		0.213268,0.226167,0.238223,0.256194,0.264266,0.284133,0.309525,0.327853,0.344547,0.395931,
		0,0.0697679,0.0769657,0.12637,0.145187,0.148123,0.16533,0.18701,0.196477,0.23,
		0.239516,0.249511,0.264943,0.285013,0.306531,0.323899,0.342622,0.372123,0,0.0815873,
		0.107465,0.115547,0.124534,0.168611,0.173604,0.199692,0.211998,0.237202,0.249887,0.271308,
		0.289257,0.317743,0.346098,0.359381,0.375028,0.403084,0,0.0490687,0.0687856,0.106853,
		0.121216,0.143824,0.178742,0.185292,0.202775,0.213988,0.234838,0.255133,0.273001,0.288946,
		0.315563,0.333694,0.345013,0.377458,0,0.0479695,0.0808138,0.11188,0.146728,0.152027,
		0.174626,0.184232,0.198756,0.215995,0.22549,0.235169,0.243321,0.259353,0.274724,0.307506,
		0.332871,0.361587,0,0.0449529,0.0986015,0.123002,0.139113,0.153449,0.177351,0.201225,
		0.222769,0.240959,0.261162,0.277178,0.287514,0.307037,0.331336,0.351689,0.367655,0.386876,
		0,0.042439,0.0808946,0.121898,0.137121,0.167179,0.182879,0.221657,0.24058,0.255993,
		0.277686,0.295097,0.313074,0.325443,0.343983,0.359623,0.380138,0.422698,0,0.0250991,
		0.0739599,0.103959,0.139442,0.161882,0.168884,0.183924,0.21372,0.245538,0.255477,0.265346,
		0.290763,0.30439,0.311178,0.345289,0.36452,0.376304,0,0.04663,0.0808946,0.114347,
		0.134077,0.157629,0.167,0.19153,0.206798,0.228047,0.252122,0.274751,0.29994,0.306671,
		0.331168,0.350741,0.376528,0.398769,0,0.0857164,0.112168,0.134713,0.137017,0.149921,
		0.172112,0.190287,0.204774,0.22102,0.232238,0.25727,0.295806,0.314814,0.330955,0.363203,
		0.394321,0.431779,0,0.0605397,0.0802486,0.102246,0.116732,0.157024,0.17849,0.189928,
		0.227519,0.236975,0.25298,0.264317,0.283736,0.316913,0.345613,0.355586,0.370968,0.414207,
		0,0.046219,0.0729577,0.0937621,0.136095,0.152127,0.162477,0.18895,0.228639,0.249349,
		0.274046,0.295235,0.298108,0.321211,0.343241,0.357891,0.392811,0.414025,0,0.0721134,
		0.0884521,0.114938,0.129931,0.148664,0.183551,0.205801,0.213666,0.237392,0.246537,0.251444,
		0.265097,0.276349,0.300502,0.343361,0.368802,0.388353,0,0.0686093,0.0907693,0.134532,
		0.156201,0.184338,0.203456,0.213064,0.229612,0.252335,0.262966,0.279079,0.293459,0.303897,
		0.326702,0.347068,0.372123,0.408023,0,0.0452401,0.0784073,0.100758,0.117411,0.140387,
		0.168407,0.179124,0.188435,0.223103,0.231905,0.239582,0.25739,0.288475,0.331242,0.342449,
		0.36208,0.441564,0,0.0342006,0.0525045,0.0933213,0.110915,0.123643,0.143974,0.180003,
		0.190006,0.207906,0.232744,0.248767,0.273716,0.29788,0.325648,0.354242,0.371275,0.407649,
		0,0.0416555,0.0831672,0.117071,0.130595,0.149736,0.16615,0.18611,0.21586,0.224561,
		0.240127,0.251626,0.259029,0.265943,0.279816,0.313708,0.332034,0.34943,0,0.0436627,
		0.0959908,0.144732,0.157547,0.189034,0.218723,0.233777,0.242015,0.25657,0.267655,0.277919,
		0.297455,0.311662,0.333536,0.342047,0.365943,0.3889,0,0.0302893,0.0583922,0.0886961,
		0.109215,0.130962,0.164465,0.18923,0.199801,0.216262,0.234449,0.251366,0.263675,0.278969,
		0.289225,0.301346,0.323937,0.372485,0,0.0665855,0.0958115,0.106797,0.1302,0.14743,
		0.163286,0.190063,0.225893,0.24058,0.253178,0.267699,0.28496,0.305326,0.309461,0.335749,
		0.364043,0.393067,0,0.0579508,0.0898262,0.114079,0.126624,0.14545,0.158913,0.18466,
		0.202401,0.226071,0.255755,0.270971,0.27937,0.303684,0.324529,0.34533,0.371747,0.399095,
		0,0.0468884,0.0650749,0.0832122,0.118122,0.12322,0.14474,0.15048,0.166178,0.194091,
		0.228737,0.24777,0.264509,0.28089,0.307177,0.329462,0.365166,0.383497,0,0.0336376,
		0.0813776,0.107706,0.12212,0.129252,0.14212,0.164928,0.192506,0.220601,0.234894,0.252876,
		0.267699,0.283718,0.299331,0.32612,0.37163,0.392683,0,0.0413221,0.0631707,0.0803647,
		0.133641,0.149613,0.166285,0.181179,0.192857,0.212679,0.23159,0.236618,0.246537,0.258054,
		0.302277,0.316412,0.341215,0.371288,0,0.049779,0.0802278,0.0980249,0.108469,0.139609,
		0.143609,0.156646,0.184483,0.199538,0.219018,0.231634,0.277911,0.300394,0.326101,0.348491,
		0.368678,0.381093,0,0.0801619,0.108604,0.121113,0.150734,0.162191,0.192966,0.219936,
		0.231159,0.240048,0.255291,0.278083,0.29326,0.308743,0.320996,0.329462,0.368814,0.390582,
		0,0.0402265,0.0815873,0.129099,0.13497,0.16438,0.182372,0.19414,0.21395,0.230952,
		0.251947,0.256205,0.262816,0.277655,0.296616,0.307983,0.323563,0.355405,0,0.024951,
		0.0475257,0.10165,0.123719,0.14876,0.192515,0.213021,0.221859,0.24544,0.265379,0.280309,
		0.292881,0.305744,0.321666,0.353104,0.375728,0.388759,0,0.0560991,0.0914044,0.108604,
		0.1302,0.158446,0.171471,0.190685,0.210552,0.224376,0.228186,0.244758,0.25925,0.287196,
		0.30532,0.326121,0.362187,0.376528,0,0.033837,0.0723276,0.105679,0.129129,0.136055,
		0.158662,0.182476,0.221779,0.238084,0.248867,0.27901,0.302263,0.3251,0.351606,0.360391,
		0.375289,0.401951,0,0.061055,0.0840054,0.0979643,0.103215,0.120297,0.155117,0.177909,
		0.198065,0.208392,0.225826,0.238402,0.2782,0.291409,0.31137,0.32745,0.342619,0.411063,
		0,0.0364692,0.0536117,0.0867045,0.108986,0.139036,0.155896,0.178986,0.202327,0.222828,
		0.24347,0.256271,0.270263,0.289225,0.314549,0.340022,0.353227,0.415296,0,0.0617275,
		0.0947721,0.122324,0.14301,0.18259,0.199063,0.22648,0.239852,0.247837,0.269127,0.274061,
		0.291275,0.301253,0.325889,0.353227,0.375728,0.38928,0,0.0337756,0.085571,0.101879,
		0.130489,0.154901,0.169206,0.200085,0.219565,0.242709,0.251784,0.270309,0.300151,0.306254,
		0.321158,0.336376,0.356754,0.37755,0,0.0516011,0.0739788,0.124086,0.142531,0.152586,
		0.166956,0.176844,0.192983,0.216591,0.234506,0.243563,0.254782,0.271109,0.296274,0.311104,
		0.335307,0.351689,0,0.0232776,0.0752992,0.113719,0.137017,0.156196,0.172237,0.194216,
		0.224798,0.232237,0.249887,0.275879,0.292839,0.302277,0.31323,0.344018,0.370328,0.385966,
		0,0.0617275,0.109125,0.132204,0.14646,0.158578,0.180758,0.207174,0.219354,0.230214,
		0.24347,0.250544,0.279909,0.28629,0.308389,0.322908,0.347665,0.40298,0,0.0524085,
		0.0696494,0.0883813,0.115202,0.138221,0.164963,0.195586,0.20829,0.222474,0.243949,0.263141,
		0.275123,0.289363,0.302447,0.323296,0.33309,0.386876,0,0.0719709,0.0933213,0.146126,
		0.181817,0.199843,0.214242,0.229,0.234749,0.24949,0.263331,0.277144,0.287651,0.298506,
		0.309833,0.325202,0.341732,0.38115,0,0.0438464,0.0750037,0.0986022,0.130161,0.139241,
		0.16267,0.193302,0.212554,0.231612,0.241946,0.250475,0.267782,0.301438,0.322309,0.340787,
		0.367288,0.392724,0,0.0659081,0.07434,0.0975843,0.137339,0.145198,0.166356,0.194954,
		0.213512,0.227185,0.248257,0.26491,0.281804,0.298498,0.30678,0.30908,0.324446,0.348072,
		0,0.0629407,0.0908245,0.113203,0.140444,0.180569,0.192675,0.209048,0.223785,0.237998,
		0.255833,0.265235,0.277911,0.291154,0.314139,0.350728,0.366161,0.391405,0,0.0634117,
		0.090707,0.123537,0.142805,0.161123,0.17086,0.178312,0.192831,0.216443,0.231634,0.244833,
		0.260519,0.273038,0.29952,0.306976,0.353034,0.372114,0,0.0588923,0.0925903,0.124565,
		0.141256,0.161039,0.186446,0.20741,0.223331,0.245508,0.2681,0.276302,0.288475,0.31475,
		0.338903,0.369103,0.382139,0.424057,0,0.0664142,0.105233,0.129099,0.142396,0.154479,
		0.17259,0.186622,0.212967,0.224211,0.256802,0.267391,0.274061,0.293023,0.316913,0.325318,
		0.343361,0.376902,0,0.065541,0.0892526,0.121242,0.142704,0.174511,0.214442,0.232814,
		0.249541,0.264509,0.288676,0.302902,0.316518,0.325223,0.34174,0.360189,0.37163,0.405106,
		0,0.0360582,0.0946066,0.134233,0.141532,0.167337,0.186342,0.202466,0.21919,0.231236,
		0.254007,0.262107,0.275273,0.287684,0.304689,0.332289,0.354753,0.378125,0,0.0612587,
		0.0825471,0.0941263,0.112697,0.13018,0.142602,0.156797,0.196346,0.228527,0.241996,0.254722,
		0.260731,0.281874,0.291951,0.320742,0.340022,0.370328,0,0.0656468,0.0884521,0.114135,
		0.14129,0.153661,0.193012,0.20556,0.216967,0.231752,0.243321,0.256209,0.268596,0.285971,
		0.303324,0.319679,0.344457,0.360061,0,0.0416497,0.0970678,0.105497,0.126624,0.139113,
		0.182598,0.21395,0.225423,0.242952,0.247077,0.262966,0.281347,0.288187,0.302902,0.322242,
		0.342622,0.385484,0,0.0329281,0.0725883,0.0930757,0.144467,0.171034,0.192158,0.21066,
		0.229359,0.248297,0.255817,0.269066,0.281393,0.290412,0.300853,0.316905,0.353854,0.384115,
		0,0.0312212,0.0765222,0.11282,0.142058,0.158719,0.174967,0.206386,0.216887,0.227999,
		0.244284,0.257621,0.27707,0.28568,0.304943,0.321204,0.341451,0.363568,0,0.0617506,
		0.083126,0.09719,0.112811,0.129835,0.14521,0.163782,0.188739,0.223785,0.240425,0.256375,
		0.27524,0.291728,0.312113,0.345292,0.359984,0.382041,0,0.0405485,0.0742342,0.124232,
		0.148102,0.155442,0.170356,0.189854,0.201169,0.208015,0.224683,0.229897,0.249811,0.260363,
		0.271288,0.288187,0.330994,0.373688,0,0.0787641,0.096637,0.116644,0.137121,0.159502,
		0.188302,0.21281,0.224211,0.238204,0.24783,0.254097,0.28102,0.309034,0.331463,0.357414,
		0.393067,0.423281,0,0.0475106,0.0786771,0.0952075,0.127137,0.14838,0.165367,0.184657,
		0.202422,0.234431,0.249669,0.271308,0.299074,0.321023,0.326773,0.346375,0.390398,0.429252,
		0,0.0400978,0.0734263,0.0853292,0.103526,0.130354,0.153706,0.181179,0.187814,0.204579,
		0.216725,0.234584,0.264943,0.286293,0.321926,0.34168,0.367316,0.389107,0,0.0364544,
		0.0690293,0.080195,0.107454,0.138046,0.168382,0.192851,0.217436,0.234349,0.265019,0.278969,
		0.303881,0.313479,0.33066,0.364618,0.379723,0.396483,0,0.0883163,0.102177,0.134532,
		0.152934,0.16267,0.182707,0.199538,0.209328,0.232752,0.245432,0.262589,0.27331,0.30762,
		0.320986,0.324545,0.346472,0.387059,0,0.0265967,0.0565911,0.0773235,0.117573,0.145477,
		0.161123,0.184196,0.203977,0.225649,0.245432,0.280868,0.303828,0.321351,0.33177,0.352691,
		0.383971,0.39902,0,0.05571,0.0801371,0.0975931,0.125707,0.140149,0.188302,0.209048,
		0.233579,0.255364,0.269127,0.289363,0.308812,0.329958,0.341961,0.357406,0.39278,0.428406,
		0,0.0752992,0.086183,0.119571,0.136234,0.157813,0.178558,0.201502,0.211424,0.220042,
		0.223103,0.239285,0.265829,0.271164,0.289791,0.317435,0.327143,0.378735,0,0.056214,
		0.0841929,0.108995,0.126205,0.143102,0.17795,0.198812,0.225893,0.247077,0.262589,0.277655,
		0.286847,0.306018,0.330727,0.353622,0.421568,0.442714,0,0.065541,0.091391,0.106674,
		0.134626,0.161952,0.186205,0.211887,0.235025,0.245698,0.253802,0.274034,0.291154,0.321187,
		0.33309,0.360391,0.387809,0.412543,0,0.0536158,0.0712803,0.0951419,0.114006,0.12598,
		0.14218,0.16602,0.174566,0.192849,0.211424,0.230422,0.247286,0.265501,0.295404,0.326702,
		0.352856,0.371733,0,0.0531495,0.0746039,0.0875586,0.113956,0.140583,0.155269,0.184338,
		0.207552,0.22282,0.235747,0.260363,0.278083,0.300245,0.327636,0.344346,0.389875,0.398848,
		0,0.0360689,0.0833079,0.123127,0.147408,0.164579,0.183297,0.200008,0.222129,0.239893,
		0.252209,0.265831,0.286011,0.300595,0.312411,0.324064,0.338678,0.371176,0,0.0667511,
		0.0905707,0.108469,0.123719,0.158995,0.171054,0.187078,0.198201,0.22161,0.244284,0.255262,
		0.285407,0.298181,0.306671,0.327974,0.355885,0.394499,0,0.0721134,0.102177,0.122766,
		0.13916,0.164467,0.189112,0.21351,0.232237,0.239536,0.246351,0.263174,0.28102,0.306764,
		0.322579,0.348989,0.379991,0.398356,0,0.034479,0.0707943,0.0813743,0.0976983,0.134104,
		0.160365,0.170356,0.197677,0.220713,0.2439,0.257524,0.275851,0.286235,0.316288,0.332036,
		0.348491,0.407649,0,0.0554812,0.0981116,0.120847,0.134853,0.152505,0.162261,0.172892,
		0.18449,0.202401,0.209825,0.245875,0.265379,0.300245,0.317324,0.332712,0.398424,0.433933,
		0,0.0632076,0.0910209,0.106674,0.12676,0.15969,0.182275,0.213512,0.235726,0.260643,
		0.270971,0.291968,0.301438,0.318725,0.334338,0.345289,0.375212,0.446205,0,0.029497,
		0.0643593,0.0777937,0.100395,0.115707,0.148456,0.183864,0.202982,0.208645,0.228186,0.244334,
		0.262905,0.279053,0.309903,0.325379,0.347601,0.368432,0,0.0513861,0.0802012,0.0901992,
		0.109027,0.132558,0.152127,0.168884,0.181816,0.200452,0.228227,0.239761,0.252728,0.265235,
		0.298221,0.310575,0.348482,0.371747,0,0.0500907,0.0832452,0.105246,0.110912,0.136692,
		0.153999,0.168382,0.197856,0.211237,0.225521,0.237365,0.27365,0.286053,0.317661,0.347827,
		0.366147,0.404841,0,0.0486888,0.0831672,0.108101,0.137356,0.156318,0.162227,0.182476,
		0.198842,0.210745,0.227579,0.246201,0.251293,0.272529,0.283823,0.314139,0.334338,0.364384,
		0,0.042388,0.0617961,0.0750445,0.114528,0.148899,0.185252,0.194757,0.217968,0.223826,
		0.238755,0.256194,0.273038,0.285545,0.292159,0.310575,0.354948,0.399245,0,0.0621465,
		0.079414,0.099047,0.117271,0.136681,0.162782,0.178312,0.190443,0.214774,0.238242,0.263331,
		0.296029,0.309023,0.324772,0.348911,0.372485,0.399475,0,0.0309707,0.0606917,0.0796688,
		0.142285,0.157813,0.174539,0.192739,0.210074,0.219778,0.232697,0.246656,0.262905,0.281027,
		0.315582,0.334654,0.350896,0.375351,0,0.05571,0.0768307,0.091058,0.119378,0.15249,
		0.179977,0.215995,0.233005,0.242144,0.255026,0.265501,0.287282,0.313749,0.329978,0.350741,
		0.374013,0.39902,0,0.0588923,0.0814614,0.0987311,0.140485,0.14821,0.168854,0.189339,
		0.208399,0.228825,0.243783,0.259092,0.289618,0.293229,0.304329,0.321333,0.328886,0.340794,
		0,0.047569,0.0685504,0.11188,0.127437,0.167675,0.184801,0.202935,0.221783,0.234674,
		0.248762,0.257145,0.277144,0.287213,0.313188,0.326336,0.347556,0.383828,0,0.0250242,
		0.0805244,0.106745,0.146459,0.159197,0.169321,0.189122,0.204402,0.221859,0.233377,0.239893,
		0.25398,0.267655,0.281027,0.295604,0.325219,0.354053,0,0.0576494,0.0811262,0.0993094,
		0.136293,0.155842,0.167166,0.190287,0.201742,0.226959,0.238755,0.249555,0.262419,0.287041,
		0.31407,0.328892,0.350964,0.372114,0,0.0576494,0.0841031,0.110407,0.132558,0.176027,
		0.183551,0.207484,0.220592,0.23476,0.246603,0.27244,0.287686,0.300406,0.308784,0.340912,
		0.354753,0.378035,0,0.0467784,0.0800354,0.114468,0.132144,0.150413,0.177351,0.198668,
		0.212679,0.224923,0.247105,0.26677,0.279816,0.293112,0.314585,0.349865,0.375898,0.416438,
		0,0.0531003,0.0926249,0.115117,0.143124,0.177511,0.213302,0.22982,0.244228,0.257673,
		0.278355,0.283836,0.304575,0.325318,0.341786,0.357302,0.379677,0.415124,0,0.0393065,
		0.0557629,0.0821654,0.103723,0.122103,0.148422,0.185672,0.216632,0.23486,0.251285,0.280537,
		0.291951,0.314495,0.341978,0.375411,0.388353,0.405106,0,0.075903,0.0944053,0.114513,
		0.136692,0.142074,0.153706,0.182707,0.210771,0.216876,0.238212,0.244247,0.253263,0.26825,
		0.285756,0.322309,0.331515,0.344183,0,0.0475257,0.093234,0.117452,0.157532,0.167166,
		0.188761,0.200918,0.217968,0.228186,0.237811,0.253768,0.274424,0.28318,0.291777,0.332697,
		0.366667,0.38743,0,0.0522273,0.103182,0.121522,0.140288,0.160084,0.169149,0.201486,
		0.22071,0.233233,0.23628,0.256468,0.267409,0.313074,0.320241,0.333626,0.353055,0.365672,
		0,0.0540365,0.105233,0.121361,0.148013,0.15786,0.179818,0.194952,0.21281,0.216294,
		0.228119,0.246288,0.26506,0.286028,0.307763,0.325202,0.347601,0.365742,0,0.0302893,
		0.061305,0.0767373,0.114462,0.134457,0.1512,0.176844,0.206438,0.222129,0.231992,0.253629,
		0.26923,0.296708,0.31017,0.342047,0.363807,0.376902,0,0.071147,0.0833908,0.111152,
		0.127163,0.142415,0.16225,0.182879,0.231159,0.242121,0.265612,0.279533,0.286689,0.30252,
		0.324783,0.342452,0.359968,0.395432,0,0.029497,0.0576386,0.0807936,0.118092,0.126353,
		0.160694,0.174566,0.183297,0.217026,0.232238,0.25505,0.271164,0.290499,0.307763,0.331515,
		0.369963,0.408905,0,0.0371499,0.0606709,0.0905707,0.119229,0.144741,0.172471,0.187441,
		0.205142,0.225826,0.233886,0.246917,0.279345,0.294288,0.311665,0.335673,0.349749,0.371903
	};
	
	
	//for every test point:
	unsigned counter = 0;
	for ( size_t pts = 0 ; pts != test_points.size() ; ++pts )
	{
		//for every distance class:
		for ( size_t dist = 0 ; dist != classes_to_compute_knn.size() ; ++dist )
		{
			std::cout << (*classes_to_compute_knn[dist])( test_points[pts]) << std::endl;
			//if ( counter % 10 == 9 )  std::cout << std::endl;
			BOOST_CHECK( fabs( (*classes_to_compute_knn[dist])( test_points[pts]) - result[counter] ) <= 5e-07 );
			++counter;
		}
	}		
}	*/






BOOST_AUTO_TEST_CASE(Distance_to_k_th_nearest_neighbor_periodic_domain_3d_brute_forcee_AAA)
{
	//first we test the brute force algorithm:

	typedef Gudhi::Topological_inference_with_cubical_complexes::Euclidan_distance_squared Euclidean_distance_squared;		
	typedef Gudhi::Topological_inference_with_cubical_complexes::periodic_domain_distance<Euclidean_distance_squared>
	periodic_Euclidean_distance_squared;
	
	
	std::vector< std::pair< double , double > > coordinates_of_grid(3);
	coordinates_of_grid[0] = std::pair< double,double >(0,1);
	coordinates_of_grid[1] = std::pair< double,double >(0,1);
	coordinates_of_grid[2] = std::pair< double,double >(0,1);
	
	Euclidean_distance_squared eu;
	periodic_Euclidean_distance_squared period_eu( coordinates_of_grid , eu );
	
	std::vector< std::vector<double> > point_cloud = 
{
{	0.2240710377	,	0.1435912591	,	0.0614581585	},

{	0.4329306507	,	0.8065080105	,	0.2185589848	},
{	0.0574821825	,	0.0626896962	,	0.8454428848	},

{	0.3465973805	,	0.9354056774	,	0.0196996354	},
{	0.1324670461	,	0.8192300708	,	0.9213577402	},
{	0.1694966145	,	0.2638772354	,	0.6784897626	},

{	0.8074921621	,	0.4387856545	,	0.7051981008	},
{	0.5657059196	,	0.3733877705	,	0.9475903716	},
{	0.7062168564	,	0.2639570518	,	0.8835326936	},

{	0.1807706819	,	0.5460551055	,	0.9522723649	},

{	0.9422669522	,	0.7820943766	,	0.2449800323	},
{	0.9724435415	,	0.4551719481	,	0.8986831212	},

{	0.5292712143	,	0.8408309834	,	0.0873390147	},
{	0.3762390846	,	0.8355073435	,	0.9812029754	},
{	0.1998295102	,	0.2246011954	,	0.1189967997	},
{	0.0178428628	,	0.7698484559	,	0.3050332498	},

{	0.1359464219	,	0.1277220107	,	0.999495666	},
{	0.1344095669	,	0.9178393232	,	0.1190430215	},

{	0.0956069152	,	0.9877411523	,	0.9099150554	},
{	0.3834235435	,	0.4715404287	,	0.7691760184	},
{	0.5806148415	,	0.2498710137	,	0.0978167353	},
{	0.7077589587	,	0.0270938836	,	0.9588638127	},

{	0.4059306134	,	0.6081185401	,	0.3019679056	},
{	0.3681147429	,	0.7604678955	,	0.5159981849	}
};	
	
	Gudhi::Topological_inference_with_cubical_complexes::Distance_to_k_th_closest_point_periodic_k_d_tree
	kd_tree( point_cloud, coordinates_of_grid , 9 );
	
	Gudhi::Topological_inference_with_cubical_complexes::Distance_to_k_th_closest_point<periodic_Euclidean_distance_squared>
	brute_force( point_cloud, period_eu , 9 );
	    
	
	std::vector< std::vector<double> > test_points =
	
	{
		{0.6570753611,0.5955838757,0.8398550118}/*
	{0.0278817846,0.2916815139,0.3887573609},{0.3553810648,0.9343284944,0.5285676813},{0.6866485006,0.8520970619,0.8668888684},
	{0.8956627958,0.1156968276,0.0739414582},{0.6355211725,0.3484526002,0.743152576},{0.0726810361,0.7873784064,0.385060332},
	{0.7665069464,0.2014935308,0.6368220472},{0.9891614937,0.7330921197,0.7102547155},
	{0.7883918153,0.014290286,0.7783693536},{0.0336034244,0.6878910807,0.2098890906},{0.8529858275,0.2412898415,0.3730617869},
	{0.6427641907,0.200828871,0.252669279},{0.1852649881,0.2579907831,0.5425100487},{0.663123504,0.5577539282,0.5905491377},
	{0.3354653688,0.8227474587,0.9516165522},{0.6034799691,0.6799063885,0.8784505285},{0.719459174,0.5953775912,0.4358448796},
	{0.42901467,0.0499089628,0.338926865},{0.6868318615,0.4042496786,0.219024647},{0.4858837686,0.3048673074,0.2064413528},
	{0.3251100981,0.6204368386,0.1215259079},{0.384246537,0.4619099458,0.577011209},{0.1957502051,0.2889291102,0.4228005044},
	{0.6831766628,0.2632829093,0.5414486954},{0.4916156209,0.6092843614,0.7656351738},{0.4322729874,0.1720699086,0.8144672341},
	{0.0017052968,0.2426031898,0.6799070754},{0.993835246,0.3934641157,0.8028928195},{0.8344388744,0.1533741863,0.052077933},
	{0.3089801904,0.0926956008,0.3230055785},{0.0271170982,0.1870652179,0.3214361172},{0.832586685,0.5277374319,0.2314551908},
	{0.9314843167,0.6746024499,0.2297013158},{0.0021536089,0.0794244623,0.2783408458},{0.2407027101,0.1303185185,0.7768267517},
	{0.8024219871,0.8166264719,0.9408042894},{0.4125754705,0.1070645854,0.3917626143},{0.7881310335,0.6119206273,0.3612281107},
	{0.2093140029,0.3191626819,0.4877916307},{0.8854867723,0.9801039936,0.1372925977},{0.8612136859,0.7057745687,0.8465774728},
	{0.8764370135,0.2966763664,0.8315236054},{0.3119052621,0.5994925359,0.3814946879},{0.7070578339,0.9780796578,0.5635364696},
	{0.5851224931,0.4823992441,0.8954798079},{0.3127624721,0.8012911966,0.4995430515},{0.5117247337,0.9834836028,0.6449389677},
	{0.2507975928,0.6793389737,0.3178833791},{0.0095620689,0.4389932367,0.9797510416},{0.9370902732,0.0666056497,0.6635932433},
	{0.8764814257,0.5038128984,0.8082919149},{0.5889022176,0.2642431161,0.7928514399},{0.9185145744,0.046059252,0.5190000213},
	{0.4452630868,0.7133746753,0.3376530202},{0.9111526795,0.5850539894,0.4890879511},{0.2747024964,0.1768817641,0.6366421564},
	{0.5592753026,0.6720825401,0.5405601789},{0.6614501821,0.5494211896,0.4482777803},{0.331644475,0.4239660688,0.2038797489},
	{0.0581143782,0.0254596837,0.7859259995},{0.1199549222,0.0386884976,0.4347598413},{0.292417774,0.3479369427,0.1531838665},
	{0.1128242332,0.7438441454,0.745681555},{0.2401035188,0.1771410096,0.807042615},{0.1234689017,0.6644280744,0.8244124809},
	{0.6452131516,0.879331308,0.2716083024},{0.9479403694,0.332957929,0.1844342141},{0.8303372548,0.1892285589,0.4289489985},
	{0.7557857418,0.0460203474,0.2242759808},{0.8151567031,0.5541405319,0.8722450782},{0.2119767461,0.9610106237,0.6250581758},
	{0.2892467363,0.073150466,0.657205262},{0.8341697995,0.5837534382,0.5041492598},{0.7957765597,0.0565284048,0.8094913445},
	{0.8098659823,0.7344695956,0.4325275582},{0.4223770299,0.7119579143,0.4925132662},{0.7018576597,0.5692464605,0.441388723},
	{0.878948973,0.4988749905,0.3432697712},{0.7926970301,0.7522489438,0.9827959056},{0.0475694353,0.0330607148,0.4641666713},
	{0.9334698259,0.3049679426,0.8354724362},{0.6555848967,0.1480603861,0.7719476013},{0.1087768308,0.9330630642,0.1765642043},
	{0.9145491763,0.3303631963,0.9045855913},{0.9143750975,0.8296238736,0.09732326},{0.0618927069,0.897331967,0.6368782204},
	{0.3150910612,0.9116620645,0.6570473495},{0.3580378599,0.2623571989,0.814647858},{0.9664658206,0.8377381989,0.7721899245},
	{0.9837020198,0.2720242415,0.4563574172},{0.91944414,0.4706337838,0.6330762212},{0.5168302623,0.1900973888,0.9053026016},
	{0.1698169338,0.8422267565,0.4352734836},{0.1704680121,0.7366362549,0.1204422733},{0.1558771764,0.3680324578,0.922587656},
	{0.377208197,0.3812016125,0.2776982288},{0.6930598826,0.2130688538,0.7388268828},{0.4868098402,0.0608007649,0.5235421255},
	{0.103365178,0.8140293758,0.7529569559},{0.5972671069,0.5497298294,0.4969925818},{0.8563019917,0.6795220347,0.2069536853},
	{0.0856518536,0.8227915512,0.2959526246},{0.0867100377,0.8575005436,0.217976413},{0.9527630396,0.8458755189,0.2894110898},
	{0.2544106271,0.4017220777,0.7145846318},{0.4149809659,0.5972550039,0.7904258594},{0.2131361254,0.3645701057,0.4913443741},
	{0.31539013,0.202036649,0.3210915248},{0.7574333525,0.0939333802,0.4100728314},{0.7516370339,0.5429121791,0.5922422118},
	{0.2195665021,0.1388090414,0.7959989444},{0.899406763,0.3000081666,0.5997520362},{0.5769034552,0.2301354171,0.6781405173},
	{0.3654143855,0.921912255,0.7039054288},{0.5972661381,0.1523829449,0.9388168044},{0.0848337784,0.6187418809,0.7124021121},
	{0.8300760316,0.7855182383,0.9261578321},{0.1363582865,0.6787659361,0.4440253698},{0.9333335564,0.2108276933,0.3914245311},
	{0.1965385748,0.5836126527,0.0487841049},{0.5397062108,0.5839183424,0.6932072763},{0.139552644,0.2591635338,0.3764203726},
	{0.6378182741,0.1488169709,0.3233298115},{0.7192198473,0.4706952481,0.3648997273},{0.947802094,0.4032824829,0.2739019038},
	{0.0970854505,0.9110400244,0.3080223324},{0.1985638938,0.3845967632,0.2895803584},{0.1684866478,0.5466994711,0.3734291738},
	{0.4505842216,0.6338260781,0.4444476345},{0.1916122308,0.5137707235,0.9476647123},{0.3408366609,0.0379224224,0.7878450612},
	{0.3489282646,0.9416533255,0.3851815299},{0.8091952498,0.1183833291,0.3632246931},{0.8706514349,0.1368547638,0.1946790947},
	{0.4391879758,0.0198513791,0.9993616869},{0.8185067559,0.8835033791,0.4361149771},{0.3187520679,0.0242436903,0.692015507},
	{0.1223861212,0.3934785489,0.051226496},{0.2450209786,0.1274023014,0.8058765808},{0.7720149281,0.7432372535,0.2363865508},
	{0.3001267086,0.7566496953,0.1571175677},{0.7458658114,0.0800891486,0.5941177257},{0.0080535775,0.4755061602,0.6350037972},
	{0.0827106193,0.6615890143,0.028390774},{0.7883587447,0.8537377736,0.1755873731},{0.9213447317,0.2877281052,0.0517333397},
	{0.858527428,0.6600245954,0.9033006672},{0.426437065,0.6410801688,0.1879744353},{0.7098625924,0.6176833413,0.6763019743},
	{0.6594055216,0.3670684479,0.0904959645},{0.2070961983,0.5454259284,0.2464962576},{0.109380943,0.4127504176,0.9255148144},
	{0.8236308896,0.2717485328,0.0889611356},{0.6204018274,0.2531314562,0.0879562788},{0.933535517,0.8034739352,0.1976074078},
	{0.1333329342,0.9984488715,0.4651318933},{0.2929060075,0.1734483752,0.3304254888},{0.4427554065,0.7031817061,0.7188385194},
	{0.7861618018,0.4392747791,0.4609192256},{0.4297306708,0.0910678173,0.0856644441},{0.3263101599,0.1900679297,0.8656203696},
	{0.9598891845,0.8739098033,0.5946957348},{0.8601532679,0.3356839109,0.445534257},{0.5811327847,0.0941536282,0.8464888963},
	{0.3362390934,0.3802505163,0.7506839186},{0.5980010226,0.8490654575,0.9566466322},{0.6723008321,0.4367131393,0.853445705},
	{0.1793271508,0.1859787316,0.1367681916},{0.4056227694,0.8398745777,0.9611199587},{0.6466252799,0.3384586412,0.1059214673},
	{0.1947013228,0.7682137974,0.615737533},{0.9453567625,0.7602070384,0.893087266},{0.3547970066,0.9525313561,0.9873204262},
	{0.0786859554,0.6566422957,0.5613909846},{0.6394291024,0.7515391838,0.9339341335},{0.0306413302,0.2682827814,0.3345994204},
	{0.3132445228,0.0521786655,0.1866439714},{0.3550304624,0.814520665,0.0284737975},{0.0295199333,0.3506491578,0.0804069594},
	{0.0772441155,0.4198115987,0.5326562512},{0.6100857158,0.0040999444,0.0831175786},{0.2230040201,0.0384767053,0.9185578863},
	{0.4115153016,0.9873048151,0.7218966454},{0.0724200262,0.6555718367,0.0178539874},{0.581852132,0.2605712679,0.1275672165},
	{0.2533805799,0.2765149043,0.3912544115},{0.3027908381,0.9986928527,0.4299810582},{0.3472439814,0.7639185048,0.4064361725},
	{0.9034222241,0.4131227473,0.6427890444},{0.7460419929,0.6069772805,0.7167975691},{0.3979357823,0.1354945437,0.9507125753},
	{0.0570875325,0.028278504,0.9812576636},{0.1229622571,0.0364863402,0.7952201411},{0.5176588467,0.7044958638,0.110351444},
	{0.5363182717,0.3783060559,0.0327049491},{0.980905558,0.2990631682,0.435821288},{0.6473647165,0.2050705738,0.3516822953},
	{0.9423140211,0.3945645185,0.4104162175},{0.1950914182,0.8122796561,0.3180316554}*/
	};	
		
	//for every test point:
	unsigned counter = 0;
	for ( size_t pts = 0 ; pts != test_points.size() ; ++pts )
	{		
		if ( brute_force( test_points[pts]) != kd_tree( test_points[pts]) )
		{
			std::cout << "test point : " << test_points[pts][0] << " " << test_points[pts][1] << " " << test_points[pts][2] << std::endl;
			std::cerr << "brute_force( test_points[pts]) : " << brute_force( test_points[pts]) << std::endl;
			std::cerr << "k_d_tree( test_points[pts]) : " << kd_tree( test_points[pts]) << std::endl;
			getchar();
		}	
	}		
}	
