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


/*

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
*/



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
	for ( size_t i = 5 ; i != 90 ; i=i+5 )
	{
		classes_to_compute_knn.push_back(
		new Gudhi::Topological_inference_with_cubical_complexes::Distance_to_k_th_closest_point<periodic_Euclidean_distance_squared>
		( point_cloud, period_eu , i ) );
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
		0.0183759,0.050156,0.0854197,0.0952759,0.116742,0.153669,0.166477,0.21269,0.221567,0.230418,
		0.239377,0.263733,0.283786,0.293594,0.321605,0.349976,0.375912,0.0433336,0.0849141,0.116313,
		0.135558,0.164437,0.177582,0.18574,0.206727,0.220485,0.247747,0.254062,0.268168,0.284448,
		0.313071,0.342483,0.379331,0.405558,0.0701411,0.0918121,0.109618,0.134716,0.160365,0.174209,
		0.203939,0.22239,0.24407,0.272176,0.281825,0.295286,0.317253,0.345377,0.36605,0.411096,
		0.434179,0.0731621,0.0948533,0.114211,0.125408,0.135625,0.155947,0.173453,0.191111,0.20949,
		0.223264,0.241813,0.253938,0.266445,0.299424,0.313299,0.335877,0.374819,0.0781912,0.0906376,
		0.110597,0.135696,0.160852,0.173948,0.198892,0.212169,0.22744,0.256435,0.275039,0.296375,
		0.320722,0.333974,0.343934,0.375764,0.38845,0.026951,0.0641567,0.0901749,0.114452,0.156911,
		0.171941,0.193983,0.205757,0.221333,0.234941,0.249287,0.265062,0.280321,0.294173,0.313646,
		0.337134,0.368897,0.0633677,0.0898456,0.118775,0.146482,0.171232,0.203563,0.212629,0.23453,
		0.244706,0.254888,0.278418,0.289729,0.303192,0.314339,0.335239,0.361461,0.379228,0.0354172,
		0.0734525,0.115046,0.157256,0.180173,0.196595,0.207793,0.221068,0.24053,0.25401,0.264277,
		0.298375,0.309899,0.357612,0.390137,0.401981,0.443137,0.0524137,0.0924799,0.115075,0.151285,
		0.179952,0.191154,0.205104,0.208745,0.227165,0.243142,0.247776,0.267231,0.278602,0.300889,
		0.342265,0.368338,0.386844,0.082854,0.113557,0.142944,0.170812,0.183613,0.197555,0.210331,
		0.225252,0.235496,0.248688,0.259262,0.279481,0.290262,0.308258,0.329787,0.353001,0.384047,
		0.036045,0.0633767,0.0785472,0.124046,0.131458,0.15699,0.17795,0.197117,0.215986,0.226483,
		0.243657,0.257206,0.279858,0.301454,0.327278,0.37453,0.393458,0.032214,0.0450644,0.0824513,
		0.110081,0.134153,0.149288,0.17288,0.195518,0.21965,0.23522,0.263319,0.284309,0.304114,
		0.32376,0.353991,0.38109,0.396045,0.0404918,0.0645427,0.102641,0.125534,0.141795,0.161711,
		0.183496,0.213371,0.228427,0.250066,0.271525,0.293449,0.314178,0.323826,0.340935,0.358331,
		0.379473,0.0548322,0.0760463,0.0872719,0.0951263,0.128908,0.169091,0.172552,0.191023,0.226863,
		0.234642,0.254567,0.264763,0.287955,0.3114,0.336675,0.373286,0.398503,0.0264551,0.0838497,
		0.12244,0.156027,0.183321,0.196407,0.223131,0.239276,0.260064,0.275484,0.286703,0.291929,
		0.332374,0.350049,0.355976,0.373356,0.406378,0.0518879,0.0731502,0.0957513,0.101714,0.120616,
		0.137544,0.161658,0.190597,0.222431,0.245803,0.260857,0.272118,0.285008,0.315308,0.335909,
		0.380657,0.409299,0.0517743,0.0647756,0.122633,0.147034,0.160533,0.175568,0.213859,0.226917,
		0.243124,0.259374,0.274943,0.309568,0.318433,0.343344,0.376622,0.401656,0.418974,0.0584093,
		0.0873295,0.115267,0.14649,0.16062,0.174997,0.205384,0.218438,0.242034,0.258583,0.279537,
		0.288089,0.295363,0.303829,0.331468,0.371085,0.435779,0.036598,0.0930345,0.114859,0.135668,
		0.150891,0.167559,0.179415,0.189934,0.210879,0.222211,0.260687,0.267593,0.287277,0.298102,
		0.316472,0.351848,0.367152,0.0444286,0.0711179,0.104641,0.124858,0.139523,0.162718,0.201008,
		0.220857,0.246115,0.251656,0.259285,0.276124,0.292382,0.30012,0.330874,0.342834,0.398748,
		0.0474221,0.0698834,0.107054,0.126297,0.156784,0.166107,0.177236,0.189338,0.225009,0.23637,
		0.263587,0.274936,0.285105,0.297695,0.318414,0.340936,0.368789,0.0442658,0.0803653,0.102251,
		0.128093,0.153926,0.182844,0.195337,0.214772,0.220044,0.227997,0.233725,0.255517,0.26965,
		0.28309,0.320436,0.346266,0.360576,0.0525553,0.0951872,0.12616,0.14015,0.152358,0.175526,
		0.185971,0.225609,0.237772,0.246985,0.272088,0.285796,0.298134,0.311944,0.329083,0.345535,
		0.376814,0.0313082,0.0696264,0.0957463,0.115328,0.148675,0.16238,0.18206,0.207722,0.222909,
		0.23237,0.242507,0.263473,0.276149,0.287288,0.30686,0.314378,0.365969,0.0480641,0.0754885,
		0.10108,0.146422,0.17138,0.181572,0.204912,0.224544,0.233708,0.245456,0.282586,0.307301,
		0.313004,0.331311,0.35122,0.392331,0.404508,0.0671232,0.0987809,0.117474,0.135886,0.168395,
		0.192357,0.211643,0.223155,0.241482,0.260738,0.268793,0.290637,0.313142,0.318527,0.338364,
		0.375414,0.410864,0.0371323,0.049734,0.0801157,0.11676,0.152013,0.168053,0.189486,0.22202,
		0.244324,0.264825,0.280685,0.297135,0.313464,0.337977,0.35854,0.375598,0.404634,0.0631982,
		0.0889066,0.100867,0.135869,0.147595,0.154582,0.176302,0.199493,0.214219,0.231054,0.241806,
		0.26401,0.290318,0.303318,0.316758,0.35353,0.376627,0.0672194,0.0801206,0.115851,0.130639,
		0.144166,0.155844,0.191913,0.215407,0.224065,0.238423,0.254188,0.265357,0.278857,0.304519,
		0.358729,0.365763,0.387898,0.06908,0.0777779,0.112258,0.131518,0.151193,0.164883,0.188706,
		0.20113,0.213458,0.230321,0.251221,0.27805,0.304053,0.307862,0.326791,0.369559,0.385289,
		0.0282755,0.0709141,0.111284,0.127006,0.140949,0.166251,0.184391,0.204219,0.216793,0.221529,
		0.234272,0.243547,0.260749,0.277108,0.301137,0.323849,0.340155,0.0430699,0.064071,0.0812654,
		0.0949799,0.122868,0.141125,0.15302,0.16967,0.189385,0.224295,0.256332,0.269294,0.28074,
		0.290134,0.317131,0.338201,0.360822,0.0618022,0.0861668,0.116056,0.125858,0.13695,0.149606,
		0.159664,0.173504,0.183437,0.189956,0.219869,0.238984,0.276521,0.321957,0.351506,0.373734,
		0.427229,0.0501169,0.0758345,0.096246,0.112073,0.136846,0.173078,0.181447,0.192743,0.209459,
		0.22707,0.239461,0.251808,0.278287,0.290159,0.315187,0.356398,0.414609,0.045956,0.0627871,
		0.0938921,0.105928,0.113787,0.131833,0.162519,0.179895,0.207667,0.225152,0.240637,0.252442,
		0.260173,0.283449,0.319892,0.337874,0.363972,0.0230056,0.0643037,0.101121,0.132733,0.153865,
		0.176125,0.185305,0.203859,0.226634,0.249492,0.280769,0.288487,0.297031,0.311581,0.327139,
		0.355252,0.378235,0.0430912,0.0936794,0.118107,0.154825,0.171654,0.192803,0.221002,0.23059,
		0.245322,0.257573,0.269542,0.285429,0.298099,0.310291,0.343752,0.372948,0.403576,0.0540567,
		0.0951568,0.122018,0.133266,0.152758,0.175133,0.19687,0.212454,0.220872,0.235354,0.249563,
		0.254713,0.262467,0.287776,0.303007,0.343136,0.389277,0.0587287,0.0846069,0.121343,0.137246,
		0.149819,0.166565,0.185433,0.200631,0.218634,0.225265,0.244027,0.256073,0.265865,0.300756,
		0.339348,0.389875,0.412834,0.0436849,0.0602848,0.0971883,0.121213,0.136553,0.159103,0.181875,
		0.203753,0.215924,0.22888,0.255348,0.270442,0.278257,0.308987,0.333222,0.349902,0.374102,
		0.0561145,0.0793565,0.103113,0.124414,0.152144,0.156687,0.17095,0.199606,0.219787,0.234407,
		0.249603,0.261792,0.284823,0.312216,0.325437,0.335484,0.374705,0.058914,0.101451,0.128956,
		0.141738,0.168254,0.201891,0.214053,0.224111,0.238482,0.248541,0.288767,0.300093,0.309768,
		0.324322,0.340811,0.351974,0.41013,0.069684,0.117674,0.127092,0.133463,0.154886,0.165166,
		0.191253,0.212156,0.223401,0.241783,0.269149,0.284528,0.293538,0.321312,0.343414,0.371537,
		0.386436,0.0410111,0.0770108,0.10837,0.141674,0.153368,0.172721,0.182769,0.192733,0.210766,
		0.242893,0.248462,0.259458,0.277776,0.300959,0.327088,0.343045,0.397832,0.075746,0.112508,
		0.152994,0.169419,0.183271,0.192182,0.200012,0.229489,0.239542,0.26234,0.278087,0.285547,
		0.299532,0.316097,0.325052,0.342937,0.3798,0.0568505,0.0908555,0.106362,0.138351,0.153333,
		0.169967,0.209412,0.231119,0.246289,0.261582,0.278511,0.290448,0.305389,0.315762,0.370019,
		0.385581,0.416217,0.0469136,0.0746141,0.0952416,0.135664,0.140665,0.178766,0.204778,0.223165,
		0.240898,0.262567,0.274251,0.291377,0.296859,0.308792,0.323376,0.348942,0.387432,0.0576877,
		0.0957264,0.127734,0.132302,0.153498,0.174052,0.18605,0.2021,0.219622,0.254438,0.28521,
		0.293709,0.313651,0.329152,0.355085,0.377538,0.393388,0.0290085,0.0686477,0.0904024,0.117256,
		0.137311,0.167978,0.184821,0.199111,0.231361,0.242537,0.255463,0.272599,0.28407,0.297104,
		0.325695,0.343046,0.370224,0.0397635,0.077511,0.120186,0.132964,0.156778,0.170948,0.19184,
		0.208923,0.225255,0.230447,0.243789,0.257981,0.27377,0.288198,0.317917,0.325527,0.340442,
		0.0825434,0.107836,0.125872,0.157662,0.161542,0.18898,0.199796,0.212102,0.224792,0.238337,
		0.250659,0.264685,0.286394,0.299562,0.314323,0.349633,0.36823,0.048231,0.077672,0.115867,
		0.13166,0.14546,0.174574,0.190212,0.21876,0.231986,0.24791,0.262851,0.283987,0.303624,
		0.340752,0.354301,0.374029,0.424759,0.0733051,0.0980453,0.113196,0.137149,0.143607,0.153135,
		0.189788,0.203928,0.242172,0.257076,0.276031,0.284298,0.301236,0.321448,0.330381,0.383087,
		0.430633,0.0432959,0.0926857,0.115719,0.133813,0.162835,0.175908,0.188539,0.195521,0.223389,
		0.243781,0.2507,0.276083,0.287669,0.312577,0.33307,0.3538,0.394453,0.056986,0.105666,
		0.134688,0.155021,0.163682,0.186907,0.193853,0.206537,0.226872,0.234933,0.255932,0.267824,
		0.290616,0.318909,0.330731,0.362296,0.37305,0.0434465,0.0666833,0.0811692,0.0998344,0.129987,
		0.150075,0.178549,0.202285,0.228999,0.247532,0.271072,0.279786,0.291075,0.327971,0.340329,
		0.365974,0.421268,0.0552682,0.067883,0.0927511,0.108818,0.157368,0.173847,0.192659,0.217133,
		0.228796,0.236,0.259053,0.274255,0.291243,0.309449,0.334397,0.347851,0.385243,0.0463272,
		0.0969957,0.142843,0.167481,0.195948,0.204154,0.214556,0.226546,0.240672,0.257632,0.275752,
		0.294313,0.307846,0.327886,0.345863,0.376465,0.39074,0.0516033,0.0946871,0.13139,0.146485,
		0.161933,0.183074,0.198377,0.218142,0.23701,0.248668,0.258273,0.278289,0.303549,0.328798,
		0.340406,0.376907,0.402387,0.0632555,0.0802736,0.100208,0.122046,0.134417,0.155444,0.172287,
		0.183847,0.210055,0.237385,0.245941,0.268401,0.281511,0.295032,0.306887,0.325572,0.346826,
		0.0457244,0.105376,0.130448,0.158144,0.185897,0.190812,0.199931,0.209953,0.22638,0.235898,
		0.251294,0.256923,0.277295,0.284508,0.296012,0.334649,0.368428,0.058958,0.0707226,0.0862562,
		0.111497,0.129428,0.146319,0.156663,0.170095,0.211333,0.243429,0.250409,0.261702,0.27298,
		0.304688,0.322848,0.365077,0.403764,0.0549889,0.0756352,0.0895365,0.107783,0.136396,0.157928,
		0.176623,0.18712,0.195465,0.218391,0.231926,0.242053,0.253855,0.303037,0.328417,0.337456,
		0.356178,0.0631029,0.0957855,0.114255,0.138061,0.14363,0.160491,0.180536,0.200702,0.216205,
		0.236662,0.259601,0.276716,0.298423,0.320526,0.344438,0.363771,0.378743,0.0319558,0.0591898,
		0.0961042,0.124744,0.149179,0.172466,0.182931,0.207354,0.228077,0.245264,0.268113,0.280612,
		0.299027,0.316322,0.342717,0.374282,0.388452,0.0456138,0.0764367,0.114592,0.131464,0.147096,
		0.163055,0.194683,0.206928,0.229464,0.238752,0.25625,0.276758,0.295696,0.301097,0.32061,
		0.35498,0.405465,0.0728413,0.0943603,0.112523,0.130388,0.154447,0.170106,0.197776,0.209261,
		0.229249,0.247555,0.272033,0.289901,0.303806,0.317909,0.329264,0.349262,0.39963,0.0335719,
		0.065431,0.0932403,0.11392,0.134673,0.145522,0.150511,0.194003,0.215533,0.225058,0.248293,
		0.265224,0.272617,0.297945,0.319411,0.334534,0.379182,0.0347817,0.0550451,0.0936592,0.124887,
		0.155997,0.163167,0.187212,0.207174,0.23939,0.259271,0.277748,0.290344,0.302113,0.329596,
		0.363212,0.392761,0.418568,0.0404042,0.0798069,0.0920877,0.130258,0.143753,0.159461,0.184556,
		0.197297,0.208491,0.226119,0.249759,0.273765,0.295871,0.318998,0.338154,0.369596,0.407124,
		0.0566648,0.0976747,0.108814,0.129526,0.149037,0.175067,0.183224,0.201681,0.241032,0.258876,
		0.283003,0.298701,0.324869,0.34192,0.365862,0.374398,0.403066,0.0425753,0.0614734,0.11446,
		0.128529,0.155394,0.171076,0.194635,0.203035,0.208417,0.236862,0.267307,0.27649,0.29215,
		0.304869,0.322211,0.350851,0.380552,0.0284301,0.0667816,0.105257,0.118228,0.142394,0.16424,
		0.195122,0.212529,0.226995,0.24037,0.264302,0.284574,0.29785,0.317419,0.337048,0.351453,
		0.39342,0.0531317,0.0656493,0.0979867,0.107426,0.152119,0.16687,0.194771,0.214932,0.221505,
		0.254021,0.262993,0.279104,0.298404,0.329369,0.36744,0.397281,0.420932,0.0952979,0.112101,
		0.132949,0.162404,0.175205,0.193996,0.211454,0.21751,0.233648,0.242395,0.269821,0.280955,
		0.291114,0.3015,0.338511,0.370296,0.384848,0.0655936,0.0830657,0.102511,0.137811,0.153873,
		0.167359,0.190081,0.218844,0.231135,0.247114,0.277936,0.302853,0.316698,0.337757,0.360282,
		0.387897,0.408031,0.0696748,0.092167,0.125314,0.146934,0.164461,0.179485,0.200123,0.228424,
		0.250489,0.264964,0.279722,0.296576,0.309304,0.33332,0.337939,0.359677,0.374982,0.0575944,
		0.0892923,0.116085,0.150898,0.164038,0.171576,0.1948,0.220909,0.230884,0.26506,0.279729,
		0.289337,0.299789,0.317562,0.33621,0.370949,0.424944,0.0374426,0.0814392,0.102323,0.113739,
		0.126475,0.140023,0.153616,0.16788,0.191583,0.198805,0.227061,0.246482,0.275369,0.286495,
		0.351043,0.389172,0.434464,0.0474874,0.0888205,0.118908,0.157975,0.188603,0.196854,0.203699,
		0.223392,0.230675,0.243109,0.268619,0.28404,0.300557,0.326597,0.361312,0.367757,0.394013,
		0.0672828,0.0745661,0.092963,0.111597,0.125533,0.144653,0.155148,0.182625,0.201131,0.233315,
		0.241404,0.263513,0.30615,0.328256,0.336288,0.367466,0.373472,0.0713053,0.10961,0.126978,
		0.144005,0.157493,0.164447,0.17508,0.19235,0.210703,0.247828,0.256489,0.28283,0.300482,
		0.327151,0.340344,0.352466,0.39147,0.0879174,0.112189,0.135176,0.144561,0.156346,0.183596,
		0.193246,0.209317,0.225081,0.24034,0.260008,0.273831,0.299575,0.327122,0.361353,0.38525,
		0.412472,0.0446764,0.0705274,0.0967181,0.108961,0.136944,0.146145,0.154098,0.1756,0.191598,
		0.203848,0.223812,0.235063,0.289393,0.296176,0.317556,0.355631,0.38292,0.0686833,0.102484,
		0.123527,0.14176,0.153317,0.165279,0.176563,0.202243,0.225123,0.249096,0.261093,0.280871,
		0.306437,0.325101,0.334561,0.351153,0.382699,0.0386355,0.0615877,0.105781,0.127815,0.173211,
		0.198744,0.211987,0.221692,0.235616,0.251412,0.271586,0.279534,0.288084,0.313007,0.321035,
		0.348882,0.364058,0.0639119,0.0980269,0.115029,0.135102,0.149924,0.180839,0.190455,0.207934,
		0.223521,0.239124,0.259606,0.269464,0.288757,0.306717,0.320385,0.346998,0.373416,0.0367776,
		0.0759555,0.10581,0.131481,0.141885,0.169345,0.19093,0.216963,0.23513,0.247151,0.270766,
		0.282122,0.2962,0.313961,0.330749,0.373204,0.414288,0.0361996,0.0790776,0.0967399,0.115986,
		0.135969,0.173123,0.184563,0.202301,0.239789,0.247554,0.265935,0.283096,0.289404,0.310422,
		0.339714,0.362147,0.399068,0.0645223,0.104762,0.136933,0.15527,0.169919,0.194018,0.217203,
		0.228558,0.244311,0.253957,0.262353,0.277148,0.288172,0.306851,0.327866,0.361736,0.377056,
		0.019434,0.0611104,0.0769851,0.112383,0.133212,0.156786,0.168695,0.188895,0.211367,0.229676,
		0.242844,0.268063,0.290517,0.31524,0.343051,0.350508,0.390929,0.0375679,0.0568989,0.112008,
		0.126088,0.142399,0.16574,0.195719,0.209055,0.223734,0.247832,0.260022,0.286349,0.298334,
		0.315013,0.331174,0.360408,0.402616,0.0438591,0.058596,0.085943,0.102944,0.135411,0.167502,
		0.194606,0.213601,0.243318,0.268107,0.278028,0.298248,0.313023,0.338607,0.36906,0.393389,
		0.410219,0.0279155,0.0543581,0.101675,0.125117,0.161619,0.176894,0.197527,0.207601,0.218565,
		0.234052,0.247047,0.266703,0.279215,0.309739,0.330682,0.346121,0.372843,0.0311365,0.0665587,
		0.0913395,0.109447,0.12486,0.153329,0.172252,0.188736,0.209259,0.245096,0.254932,0.269417,
		0.279777,0.291212,0.330218,0.35204,0.371515,0.0541091,0.0781329,0.11208,0.128634,0.147338,
		0.174048,0.195388,0.211804,0.222972,0.237531,0.250581,0.267409,0.279349,0.294047,0.325958,
		0.355269,0.381029,0.0568813,0.0811195,0.10948,0.123813,0.136097,0.161443,0.184892,0.20637,
		0.220622,0.235694,0.245465,0.265285,0.272836,0.30044,0.310084,0.323531,0.338605,0.0635852,
		0.128975,0.138236,0.152637,0.163128,0.169662,0.181001,0.204523,0.219178,0.232979,0.248811,
		0.283738,0.314527,0.31815,0.332819,0.395547,0.436076,0.0606917,0.0872101,0.118522,0.1437,
		0.16821,0.189876,0.206734,0.226129,0.242668,0.247195,0.261269,0.269441,0.287357,0.298631,
		0.311916,0.344182,0.382126,0.0516573,0.0962368,0.107738,0.125946,0.137595,0.178185,0.197577,
		0.207417,0.21841,0.229219,0.251894,0.287382,0.311368,0.323527,0.338282,0.350829,0.382475,
		0.0429823,0.116589,0.13537,0.161218,0.181873,0.19913,0.215,0.231286,0.23894,0.247342,
		0.255817,0.276557,0.295004,0.336326,0.353093,0.377909,0.40763,0.0437774,0.0873471,0.112953,
		0.129996,0.149958,0.170496,0.181269,0.205413,0.212885,0.222347,0.245556,0.259638,0.289259,
		0.300502,0.345816,0.354053,0.375028,0.0182351,0.0810483,0.0946979,0.105476,0.130595,0.16602,
		0.1859,0.199801,0.214786,0.232254,0.24555,0.262998,0.277686,0.307597,0.320662,0.341339,
		0.358361,0.0238025,0.083126,0.0908245,0.112428,0.125263,0.154985,0.166171,0.185487,0.201225,
		0.233377,0.245572,0.267736,0.290685,0.306254,0.324263,0.351315,0.399721,0.0254217,0.0606709,
		0.100385,0.139113,0.156797,0.176027,0.190301,0.201001,0.224683,0.236598,0.252712,0.260757,
		0.282465,0.305395,0.316905,0.343241,0.370027,0.0656468,0.0758587,0.12022,0.144788,0.14787,
		0.163892,0.186803,0.195486,0.218744,0.237597,0.245772,0.261066,0.283836,0.297696,0.318725,
		0.329331,0.364043,0.0697679,0.100397,0.115345,0.124396,0.160097,0.173061,0.195154,0.204569,
		0.22617,0.249143,0.25764,0.288327,0.303659,0.344634,0.358289,0.374295,0.400307,0.0413221,
		0.0671446,0.0926959,0.119182,0.1411,0.167246,0.184673,0.20159,0.210919,0.2335,0.254968,
		0.271204,0.287514,0.314985,0.328676,0.341006,0.359842,0.0405382,0.0730318,0.111152,0.143112,
		0.150991,0.169206,0.182656,0.19843,0.215006,0.222811,0.233877,0.240588,0.252617,0.274436,
		0.300571,0.320439,0.356682,0.034198,0.092681,0.122671,0.138406,0.153367,0.174449,0.194269,
		0.221541,0.234801,0.259002,0.27679,0.28525,0.299955,0.321325,0.347045,0.36452,0.385617,
		0.032949,0.0762426,0.120847,0.136494,0.159318,0.178064,0.219645,0.235987,0.254968,0.277571,
		0.289337,0.305278,0.324466,0.337103,0.352489,0.38003,0.421366,0.024951,0.0729372,0.103741,
		0.13661,0.160528,0.168567,0.178127,0.212368,0.238717,0.253089,0.262281,0.282942,0.300474,
		0.308604,0.342449,0.356115,0.372193,0.0438464,0.0719709,0.109215,0.13145,0.14956,0.164209,
		0.185244,0.204919,0.223854,0.250523,0.265019,0.298928,0.305712,0.325412,0.348072,0.370968,
		0.391405,0.0825069,0.109226,0.126773,0.136214,0.149318,0.17037,0.189122,0.202935,0.219778,
		0.23138,0.247958,0.283153,0.31312,0.328724,0.359866,0.381204,0.422482,0.0540477,0.0799437,
		0.0995883,0.115345,0.144732,0.172892,0.183303,0.213113,0.236045,0.252728,0.261907,0.27848,
		0.305278,0.341025,0.35196,0.368432,0.413581,0.0429749,0.0649845,0.0923199,0.134886,0.149507,
		0.162055,0.181766,0.228054,0.246436,0.259016,0.292552,0.298096,0.311936,0.338366,0.353093,
		0.368549,0.412401,0.0719394,0.0800354,0.101689,0.127091,0.141908,0.182823,0.205142,0.211237,
		0.223211,0.245174,0.249111,0.263141,0.275812,0.292441,0.341732,0.364242,0.387486,0.0588266,
		0.0896343,0.12081,0.151845,0.179223,0.199778,0.212554,0.22617,0.249877,0.257634,0.275123,
		0.288804,0.303532,0.326695,0.33826,0.371608,0.399948,0.0408899,0.0782837,0.0942494,0.115202,
		0.134635,0.164333,0.178285,0.186248,0.214121,0.231893,0.23893,0.251481,0.282123,0.317101,
		0.339093,0.353034,0.419429,0.0247515,0.0510598,0.0905373,0.10723,0.122324,0.134969,0.178449,
		0.184837,0.203502,0.229405,0.245213,0.268701,0.296626,0.321211,0.353104,0.36781,0.395715,
		0.0406598,0.0788507,0.115703,0.126865,0.148123,0.164925,0.184304,0.214826,0.22361,0.237934,
		0.250989,0.256843,0.264579,0.277264,0.298339,0.328892,0.34901,0.0429823,0.086183,0.143102,
		0.155269,0.179977,0.20849,0.233577,0.241447,0.254625,0.266579,0.271759,0.294998,0.311178,
		0.328886,0.338601,0.361918,0.385212,0.0297326,0.0540998,0.083523,0.109027,0.129129,0.153852,
		0.1869,0.196998,0.21586,0.228344,0.247754,0.263283,0.277571,0.286302,0.291752,0.322435,
		0.365943,0.0539442,0.0897239,0.104567,0.125302,0.146155,0.161041,0.173448,0.222517,0.239487,
		0.252943,0.263264,0.281516,0.304864,0.309293,0.327143,0.349615,0.3876,0.0576386,0.0896137,
		0.11196,0.119349,0.143124,0.156404,0.176243,0.201001,0.223559,0.255407,0.268326,0.276349,
		0.300582,0.309023,0.344806,0.364032,0.397675,0.0417095,0.0634765,0.08016,0.108358,0.123127,
		0.140387,0.147878,0.164312,0.178082,0.228119,0.238212,0.264436,0.278355,0.302501,0.324545,
		0.355423,0.378782,0.0254217,0.0739904,0.102383,0.121792,0.128919,0.140626,0.159228,0.184304,
		0.220421,0.233783,0.252764,0.265752,0.275974,0.296616,0.316556,0.342613,0.389735,0.0342126,
		0.0626945,0.0760776,0.128034,0.142805,0.163876,0.180103,0.18854,0.210951,0.229253,0.23628,
		0.242857,0.254157,0.2936,0.313397,0.340697,0.362404,0.0490687,0.0791758,0.0922137,0.108199,
		0.129416,0.143112,0.152048,0.182894,0.19913,0.216031,0.225534,0.275446,0.299995,0.313188,
		0.347549,0.359866,0.375305,0.0723149,0.106286,0.120741,0.149613,0.158578,0.190934,0.214815,
		0.231095,0.234739,0.249811,0.260344,0.293112,0.307037,0.320104,0.328273,0.347247,0.39001,
		0.0392409,0.0777691,0.124351,0.134304,0.152325,0.178082,0.190839,0.213003,0.228026,0.251366,
		0.255753,0.260753,0.27665,0.295235,0.307923,0.321158,0.338041,0.0211083,0.0472769,0.099423,
		0.118591,0.145187,0.171471,0.210919,0.220592,0.244606,0.258733,0.274242,0.289017,0.301421,
		0.319402,0.325585,0.354953,0.388692,0.0529063,0.0886211,0.107023,0.121299,0.150169,0.163345,
		0.188686,0.208569,0.222327,0.227948,0.243561,0.259059,0.281091,0.303144,0.324196,0.359329,
		0.3711,0.0325257,0.0673999,0.100321,0.128919,0.132222,0.158205,0.173906,0.219239,0.235987,
		0.244658,0.271109,0.291777,0.320241,0.347686,0.3575,0.369995,0.398236,0.0481098,0.0832122,
		0.0966171,0.102572,0.118092,0.149805,0.16948,0.197187,0.206881,0.218021,0.237999,0.26992,
		0.28863,0.306098,0.32708,0.337103,0.40836,0.0349807,0.0519815,0.0788119,0.103741,0.117926,
		0.154594,0.177361,0.201774,0.219277,0.242651,0.25398,0.268596,0.283823,0.310947,0.339371,
		0.351005,0.414781,0.0605683,0.0946979,0.120194,0.140154,0.164058,0.199,0.217363,0.23727,
		0.247754,0.258088,0.273642,0.285289,0.297307,0.321723,0.345586,0.370307,0.38574,0.0290423,
		0.0695697,0.0975931,0.119177,0.153677,0.167246,0.195129,0.214121,0.240195,0.249569,0.260409,
		0.288765,0.303927,0.318506,0.335048,0.352125,0.368061,0.0416555,0.0723032,0.117411,0.139442,
		0.149805,0.16255,0.174698,0.190837,0.215152,0.232367,0.243511,0.252533,0.269066,0.290167,
		0.309034,0.331857,0.351417,0.0210896,0.0671188,0.0867045,0.132566,0.155297,0.171647,0.191729,
		0.222284,0.230467,0.245539,0.275812,0.291409,0.302123,0.310163,0.339788,0.368678,0.383497,
		0.0460161,0.108245,0.125355,0.145685,0.158205,0.180017,0.196589,0.217611,0.229142,0.242125,
		0.250036,0.279375,0.285545,0.302115,0.319679,0.344715,0.391573,0.0522273,0.0644326,0.0878344,
		0.114636,0.137458,0.163777,0.18116,0.205413,0.221581,0.240666,0.26264,0.273508,0.288273,
		0.300582,0.32137,0.330416,0.374469,0.0688875,0.091026,0.138104,0.178345,0.195129,0.206663,
		0.228208,0.234449,0.249421,0.261903,0.274724,0.284741,0.296236,0.308743,0.323814,0.34168,
		0.359488,0.0432181,0.072469,0.0942494,0.128388,0.134969,0.161952,0.193042,0.206527,0.226769,
		0.241584,0.248871,0.265148,0.28629,0.317601,0.331733,0.365457,0.390148,0.0468986,0.0739788,
		0.0975462,0.119229,0.144035,0.160365,0.194193,0.213419,0.223368,0.247286,0.256271,0.276607,
		0.292839,0.304239,0.308167,0.323349,0.343683,0.0400464,0.0872918,0.105679,0.131562,0.176364,
		0.187441,0.205481,0.223778,0.234739,0.255657,0.264747,0.275425,0.290036,0.311169,0.338878,
		0.365686,0.384544,0.0516011,0.0841388,0.121743,0.140507,0.159281,0.169571,0.176513,0.190916,
		0.208674,0.228787,0.240078,0.25593,0.271894,0.294544,0.306018,0.340094,0.371927,0.0573174,
		0.0893611,0.121898,0.134304,0.160515,0.181114,0.19993,0.213414,0.242917,0.259057,0.269807,
		0.286408,0.314549,0.329938,0.363807,0.382041,0.402921,0.060685,0.099745,0.128034,0.138859,
		0.151128,0.172454,0.181816,0.204179,0.222811,0.253464,0.265714,0.271894,0.290439,0.310983,
		0.324915,0.337323,0.363962,0.0494473,0.0840409,0.118664,0.141908,0.167875,0.213152,0.22771,
		0.245556,0.263719,0.279533,0.299955,0.313696,0.324772,0.339371,0.357185,0.37155,0.403099,
		0.0187305,0.0896137,0.12895,0.140149,0.160515,0.183142,0.198756,0.216591,0.229377,0.252317,
		0.260753,0.269787,0.285293,0.302927,0.327982,0.35332,0.372877,0.060685,0.0765356,0.0929612,
		0.106286,0.128017,0.142014,0.152775,0.195445,0.215,0.239536,0.253903,0.258472,0.281393,
		0.290256,0.311662,0.329016,0.368441,0.0519857,0.0830867,0.103959,0.13018,0.153039,0.162921,
		0.20159,0.214424,0.229596,0.239761,0.255993,0.268445,0.28385,0.298435,0.313749,0.335703,
		0.355586,0.0360689,0.0880325,0.104567,0.12538,0.136162,0.181269,0.193558,0.222767,0.238223,
		0.246762,0.260265,0.278542,0.288163,0.300571,0.313063,0.334684,0.378436,0.0257914,0.0665855,
		0.091391,0.1399,0.17037,0.191809,0.201439,0.227611,0.236288,0.253083,0.26634,0.278751,
		0.289561,0.300151,0.314464,0.346928,0.380076,0.0267962,0.0734752,0.112609,0.139658,0.157978,
		0.172724,0.206153,0.216329,0.226829,0.244244,0.256716,0.269416,0.283669,0.304598,0.319914,
		0.337525,0.357094,0.0601116,0.0794563,0.0953404,0.110513,0.129084,0.145106,0.161302,0.185252,
		0.215238,0.239516,0.248008,0.27365,0.290439,0.310947,0.327486,0.355117,0.37326,0.0357931,
		0.0687856,0.119622,0.140164,0.154985,0.167103,0.178206,0.198603,0.204774,0.223863,0.229377,
		0.246436,0.256843,0.270462,0.283669,0.323557,0.371663,0.0763898,0.0872851,0.115966,0.135279,
		0.158895,0.183784,0.208755,0.219194,0.233886,0.247292,0.253638,0.277719,0.303529,0.325219,
		0.350728,0.385261,0.409837,0.0291784,0.0623985,0.0902539,0.125302,0.143008,0.164465,0.180279,
		0.201227,0.223693,0.249637,0.265334,0.29746,0.31511,0.325035,0.342783,0.389107,0.428717,
		0.034479,0.0689162,0.0808138,0.10312,0.12895,0.143426,0.17849,0.187099,0.196998,0.212967,
		0.232237,0.260862,0.285187,0.317045,0.341451,0.366044,0.384483,0.0294083,0.0689162,0.0745001,
		0.100233,0.129658,0.162418,0.183342,0.208173,0.233903,0.261836,0.278241,0.295018,0.308237,
		0.32612,0.355423,0.377994,0.394346,0.0667511,0.0945766,0.131892,0.151258,0.160854,0.180758,
		0.198861,0.206438,0.228344,0.239398,0.255364,0.270462,0.290256,0.316288,0.324077,0.344457,
		0.381615,0.0238723,0.0467784,0.0762426,0.102572,0.144796,0.161075,0.17516,0.201822,0.21962,
		0.240588,0.277114,0.298435,0.31475,0.330296,0.348834,0.374607,0.397804,0.0490392,0.0745001,
		0.0965316,0.117427,0.13655,0.163576,0.206881,0.225486,0.254142,0.268643,0.286144,0.307506,
		0.323557,0.338041,0.355117,0.390582,0.427046,0.0739599,0.0857164,0.118591,0.130489,0.155453,
		0.171156,0.200052,0.21066,0.216782,0.221262,0.231992,0.26336,0.268483,0.286689,0.308169,
		0.326993,0.375398,0.0509804,0.0724875,0.106456,0.124534,0.13805,0.17524,0.198384,0.224798,
		0.243483,0.254252,0.275238,0.285407,0.30601,0.329016,0.34971,0.411144,0.432352,0.0649271,
		0.0870339,0.105582,0.126686,0.161622,0.180779,0.199843,0.232237,0.24525,0.250544,0.273206,
		0.288595,0.304864,0.332626,0.353163,0.386162,0.401151,0.0527461,0.0680186,0.0849551,0.111683,
		0.123071,0.142014,0.158762,0.172756,0.190301,0.210425,0.221826,0.2464,0.26541,0.29173,
		0.323812,0.346992,0.370356,0.0371139,0.0745062,0.0836585,0.108602,0.138046,0.152325,0.172471,
		0.207484,0.220601,0.235251,0.260352,0.275622,0.289259,0.312893,0.343983,0.377589,0.397997,
		0.0191161,0.0809807,0.113651,0.146795,0.160939,0.182366,0.199615,0.220172,0.236618,0.24817,
		0.264153,0.285001,0.297566,0.310979,0.322435,0.337145,0.370027,0.0605397,0.0867582,0.108266,
		0.120886,0.147652,0.170453,0.186127,0.197677,0.221581,0.237506,0.253555,0.27881,0.295604,
		0.301346,0.323195,0.353004,0.380305,0.0468986,0.0953353,0.118519,0.138525,0.155171,0.1859,
		0.20991,0.231504,0.238402,0.245657,0.25764,0.27706,0.29994,0.322358,0.332654,0.367859,
		0.389632,0.0294843,0.0701895,0.0770087,0.0944053,0.124263,0.158446,0.168854,0.184525,0.216294,
		0.239285,0.25657,0.275273,0.282131,0.31137,0.331242,0.345663,0.375289,0.0432181,0.096997,
		0.119332,0.129252,0.150753,0.160801,0.169078,0.182372,0.200918,0.208399,0.235204,0.260613,
		0.298814,0.314868,0.327395,0.371565,0.424312,0.0573174,0.0892526,0.0959908,0.124263,0.158719,
		0.172869,0.206022,0.234749,0.257262,0.26992,0.289116,0.29931,0.31323,0.334241,0.342783,
		0.374832,0.441573,0.0287714,0.0624478,0.0723276,0.0980249,0.110737,0.138872,0.183152,0.194421,
		0.207161,0.219645,0.244244,0.255833,0.274424,0.302115,0.325354,0.342624,0.357891,0.0405382,
		0.0697778,0.0887744,0.108371,0.12224,0.148456,0.168489,0.180871,0.199615,0.22727,0.236971,
		0.25218,0.264436,0.297474,0.308389,0.348082,0.357414,0.0452574,0.0796688,0.103979,0.109754,
		0.127034,0.145477,0.162435,0.189112,0.210954,0.223863,0.236978,0.26634,0.279053,0.29881,
		0.344697,0.365205,0.395187,0.0468884,0.0734818,0.102219,0.131344,0.154176,0.161882,0.180871,
		0.191633,0.210361,0.214434,0.244764,0.250523,0.265392,0.282835,0.296982,0.327395,0.362187,
		0.0337862,0.0569165,0.0719394,0.0902539,0.148102,0.178634,0.19414,0.213983,0.223502,0.236974,
		0.248561,0.267171,0.283108,0.291275,0.310351,0.349677,0.374153,0.0429749,0.0770087,0.0979643,
		0.113906,0.134233,0.162418,0.178094,0.18502,0.212445,0.231308,0.253089,0.293103,0.304335,
		0.321223,0.343683,0.369963,0.390484,0.0294843,0.0536117,0.0768307,0.14098,0.155171,0.17379,
		0.183778,0.207812,0.218438,0.227737,0.242969,0.258054,0.276302,0.304549,0.328843,0.348086,
		0.364699,0.0544003,0.0748759,0.0882345,0.108995,0.143438,0.170193,0.203452,0.229856,0.24058,
		0.248297,0.26506,0.276557,0.312623,0.326061,0.341535,0.368443,0.390087,0.0426603,0.0810284,
		0.0948114,0.129712,0.147516,0.163175,0.186248,0.203733,0.228016,0.239398,0.255347,0.27937,
		0.292298,0.303611,0.315859,0.328676,0.340651,0.0249337,0.0600723,0.1067,0.121743,0.166112,
		0.184196,0.201486,0.220325,0.234601,0.244334,0.256964,0.268258,0.286408,0.301215,0.325355,
		0.341025,0.382311,0.0193913,0.0784439,0.09719,0.12322,0.158875,0.168713,0.183924,0.204179,
		0.218973,0.227241,0.237815,0.249637,0.263512,0.279654,0.289616,0.322278,0.347916,0.0533644,
		0.0695697,0.0973604,0.133046,0.15518,0.165633,0.18449,0.200784,0.225534,0.23791,0.248275,
		0.2612,0.282597,0.309733,0.323349,0.340651,0.37155,0.0531367,0.0802278,0.107643,0.120238,
		0.153548,0.183137,0.207241,0.216725,0.231308,0.239665,0.265829,0.279995,0.297266,0.308604,
		0.334854,0.351158,0.375898,0.042439,0.0678254,0.105002,0.130842,0.141349,0.173906,0.193594,
		0.211665,0.224356,0.242569,0.263141,0.278241,0.291406,0.306712,0.326149,0.371733,0.411532,
		0.0502504,0.0827703,0.114938,0.13805,0.172112,0.204713,0.22892,0.238885,0.256904,0.278276,
		0.282465,0.299995,0.315687,0.339021,0.354953,0.37644,0.402809,0.0366766,0.0461328,0.080747,
		0.0875586,0.115282,0.147233,0.17547,0.206022,0.228739,0.250989,0.268561,0.288408,0.309832,
		0.340794,0.367288,0.387134,0.400426,0.0729372,0.0921201,0.113503,0.128606,0.140485,0.153661,
		0.181114,0.198603,0.216262,0.228639,0.243483,0.252638,0.262598,0.284863,0.314868,0.331216,
		0.342597,0.0390168,0.091058,0.110513,0.157178,0.165912,0.184545,0.197216,0.215874,0.225254,
		0.233302,0.252093,0.263674,0.280868,0.290036,0.321223,0.365311,0.380076,0.0509804,0.0998073,
		0.121113,0.138096,0.158895,0.168489,0.197739,0.216553,0.232878,0.235647,0.24645,0.265334,
		0.302501,0.318863,0.332871,0.353004,0.365311,0.0257914,0.0970678,0.119571,0.14743,0.157547,
		0.175262,0.189854,0.19993,0.214442,0.226679,0.242173,0.263468,0.279749,0.292298,0.324121,
		0.345613,0.364242,0.0160165,0.0582143,0.0760776,0.102967,0.131344,0.147701,0.168407,0.199063,
		0.218294,0.22892,0.253083,0.265747,0.288595,0.307789,0.332739,0.361959,0.370454,0.0575998,
		0.0791604,0.104538,0.127093,0.139824,0.161894,0.175744,0.227832,0.240826,0.263662,0.271755,
		0.28523,0.291728,0.32449,0.341339,0.359511,0.390479,0.0289331,0.0558541,0.0803647,0.111065,
		0.12456,0.150629,0.174511,0.183142,0.212572,0.231147,0.246917,0.270959,0.289687,0.304892,
		0.323296,0.366675,0.400116,0.0333046,0.0601116,0.0887744,0.113103,0.140164,0.166392,0.186773,
		0.20422,0.224274,0.232744,0.245657,0.2642,0.292881,0.309293,0.329331,0.348656,0.364618
	};
	
	
	//for every test point:
	unsigned counter = 0;
	for ( size_t pts = 0 ; pts != test_points.size() ; ++pts )
	{
		//for every distance class:
		for ( size_t dist = 0 ; dist != classes_to_compute_knn.size() ; ++dist )
		{
			//std::cout << (*classes_to_compute_knn[dist])( test_points[pts]) << ",";
			//if ( counter % 10 == 9 )  std::cout << std::endl;
			BOOST_CHECK( fabs( (*classes_to_compute_knn[dist])( test_points[pts]) - result[counter] ) <= 5e-07 );
			++counter;
		}
	}		
}	




BOOST_AUTO_TEST_CASE(Distance_to_k_th_nearest_neighbor_periodic_domain_3d__kd_trees)
{
	//and  now the same with kd trees
	
	std::vector< std::pair< double , double > > coordinates_of_grid(3);
	coordinates_of_grid[0] = std::pair< double,double >(0,1);
	coordinates_of_grid[1] = std::pair< double,double >(0,1);
	coordinates_of_grid[2] = std::pair< double,double >(0,1);
	
	
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
	
	std::vector< Gudhi::Topological_inference_with_cubical_complexes::Distance_to_k_th_closest_point_periodic_k_d_tree* >
	classes_to_compute_knn;
	classes_to_compute_knn.reserve(91);
	
		
	for ( size_t i = 5 ; i != 90 ; i=i+5 )
	{
		classes_to_compute_knn.push_back(
		new Gudhi::Topological_inference_with_cubical_complexes::Distance_to_k_th_closest_point_periodic_k_d_tree( point_cloud, coordinates_of_grid , i ) );
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
		0.0183759,0.050156,0.0854197,0.0952759,0.116742,0.153669,0.166477,0.21269,0.221567,0.230418,
		0.239377,0.263733,0.283786,0.293594,0.321605,0.349976,0.375912,0.0433336,0.0849141,0.116313,
		0.135558,0.164437,0.177582,0.18574,0.206727,0.220485,0.247747,0.254062,0.268168,0.284448,
		0.313071,0.342483,0.379331,0.405558,0.0701411,0.0918121,0.109618,0.134716,0.160365,0.174209,
		0.203939,0.22239,0.24407,0.272176,0.281825,0.295286,0.317253,0.345377,0.36605,0.411096,
		0.434179,0.0731621,0.0948533,0.114211,0.125408,0.135625,0.155947,0.173453,0.191111,0.20949,
		0.223264,0.241813,0.253938,0.266445,0.299424,0.313299,0.335877,0.374819,0.0781912,0.0906376,
		0.110597,0.135696,0.160852,0.173948,0.198892,0.212169,0.22744,0.256435,0.275039,0.296375,
		0.320722,0.333974,0.343934,0.375764,0.38845,0.026951,0.0641567,0.0901749,0.114452,0.156911,
		0.171941,0.193983,0.205757,0.221333,0.234941,0.249287,0.265062,0.280321,0.294173,0.313646,
		0.337134,0.368897,0.0633677,0.0898456,0.118775,0.146482,0.171232,0.203563,0.212629,0.23453,
		0.244706,0.254888,0.278418,0.289729,0.303192,0.314339,0.335239,0.361461,0.379228,0.0354172,
		0.0734525,0.115046,0.157256,0.180173,0.196595,0.207793,0.221068,0.24053,0.25401,0.264277,
		0.298375,0.309899,0.357612,0.390137,0.401981,0.443137,0.0524137,0.0924799,0.115075,0.151285,
		0.179952,0.191154,0.205104,0.208745,0.227165,0.243142,0.247776,0.267231,0.278602,0.300889,
		0.342265,0.368338,0.386844,0.082854,0.113557,0.142944,0.170812,0.183613,0.197555,0.210331,
		0.225252,0.235496,0.248688,0.259262,0.279481,0.290262,0.308258,0.329787,0.353001,0.384047,
		0.036045,0.0633767,0.0785472,0.124046,0.131458,0.15699,0.17795,0.197117,0.215986,0.226483,
		0.243657,0.257206,0.279858,0.301454,0.327278,0.37453,0.393458,0.032214,0.0450644,0.0824513,
		0.110081,0.134153,0.149288,0.17288,0.195518,0.21965,0.23522,0.263319,0.284309,0.304114,
		0.32376,0.353991,0.38109,0.396045,0.0404918,0.0645427,0.102641,0.125534,0.141795,0.161711,
		0.183496,0.213371,0.228427,0.250066,0.271525,0.293449,0.314178,0.323826,0.340935,0.358331,
		0.379473,0.0548322,0.0760463,0.0872719,0.0951263,0.128908,0.169091,0.172552,0.191023,0.226863,
		0.234642,0.254567,0.264763,0.287955,0.3114,0.336675,0.373286,0.398503,0.0264551,0.0838497,
		0.12244,0.156027,0.183321,0.196407,0.223131,0.239276,0.260064,0.275484,0.286703,0.291929,
		0.332374,0.350049,0.355976,0.373356,0.406378,0.0518879,0.0731502,0.0957513,0.101714,0.120616,
		0.137544,0.161658,0.190597,0.222431,0.245803,0.260857,0.272118,0.285008,0.315308,0.335909,
		0.380657,0.409299,0.0517743,0.0647756,0.122633,0.147034,0.160533,0.175568,0.213859,0.226917,
		0.243124,0.259374,0.274943,0.309568,0.318433,0.343344,0.376622,0.401656,0.418974,0.0584093,
		0.0873295,0.115267,0.14649,0.16062,0.174997,0.205384,0.218438,0.242034,0.258583,0.279537,
		0.288089,0.295363,0.303829,0.331468,0.371085,0.435779,0.036598,0.0930345,0.114859,0.135668,
		0.150891,0.167559,0.179415,0.189934,0.210879,0.222211,0.260687,0.267593,0.287277,0.298102,
		0.316472,0.351848,0.367152,0.0444286,0.0711179,0.104641,0.124858,0.139523,0.162718,0.201008,
		0.220857,0.246115,0.251656,0.259285,0.276124,0.292382,0.30012,0.330874,0.342834,0.398748,
		0.0474221,0.0698834,0.107054,0.126297,0.156784,0.166107,0.177236,0.189338,0.225009,0.23637,
		0.263587,0.274936,0.285105,0.297695,0.318414,0.340936,0.368789,0.0442658,0.0803653,0.102251,
		0.128093,0.153926,0.182844,0.195337,0.214772,0.220044,0.227997,0.233725,0.255517,0.26965,
		0.28309,0.320436,0.346266,0.360576,0.0525553,0.0951872,0.12616,0.14015,0.152358,0.175526,
		0.185971,0.225609,0.237772,0.246985,0.272088,0.285796,0.298134,0.311944,0.329083,0.345535,
		0.376814,0.0313082,0.0696264,0.0957463,0.115328,0.148675,0.16238,0.18206,0.207722,0.222909,
		0.23237,0.242507,0.263473,0.276149,0.287288,0.30686,0.314378,0.365969,0.0480641,0.0754885,
		0.10108,0.146422,0.17138,0.181572,0.204912,0.224544,0.233708,0.245456,0.282586,0.307301,
		0.313004,0.331311,0.35122,0.392331,0.404508,0.0671232,0.0987809,0.117474,0.135886,0.168395,
		0.192357,0.211643,0.223155,0.241482,0.260738,0.268793,0.290637,0.313142,0.318527,0.338364,
		0.375414,0.410864,0.0371323,0.049734,0.0801157,0.11676,0.152013,0.168053,0.189486,0.22202,
		0.244324,0.264825,0.280685,0.297135,0.313464,0.337977,0.35854,0.375598,0.404634,0.0631982,
		0.0889066,0.100867,0.135869,0.147595,0.154582,0.176302,0.199493,0.214219,0.231054,0.241806,
		0.26401,0.290318,0.303318,0.316758,0.35353,0.376627,0.0672194,0.0801206,0.115851,0.130639,
		0.144166,0.155844,0.191913,0.215407,0.224065,0.238423,0.254188,0.265357,0.278857,0.304519,
		0.358729,0.365763,0.387898,0.06908,0.0777779,0.112258,0.131518,0.151193,0.164883,0.188706,
		0.20113,0.213458,0.230321,0.251221,0.27805,0.304053,0.307862,0.326791,0.369559,0.385289,
		0.0282755,0.0709141,0.111284,0.127006,0.140949,0.166251,0.184391,0.204219,0.216793,0.221529,
		0.234272,0.243547,0.260749,0.277108,0.301137,0.323849,0.340155,0.0430699,0.064071,0.0812654,
		0.0949799,0.122868,0.141125,0.15302,0.16967,0.189385,0.224295,0.256332,0.269294,0.28074,
		0.290134,0.317131,0.338201,0.360822,0.0618022,0.0861668,0.116056,0.125858,0.13695,0.149606,
		0.159664,0.173504,0.183437,0.189956,0.219869,0.238984,0.276521,0.321957,0.351506,0.373734,
		0.427229,0.0501169,0.0758345,0.096246,0.112073,0.136846,0.173078,0.181447,0.192743,0.209459,
		0.22707,0.239461,0.251808,0.278287,0.290159,0.315187,0.356398,0.414609,0.045956,0.0627871,
		0.0938921,0.105928,0.113787,0.131833,0.162519,0.179895,0.207667,0.225152,0.240637,0.252442,
		0.260173,0.283449,0.319892,0.337874,0.363972,0.0230056,0.0643037,0.101121,0.132733,0.153865,
		0.176125,0.185305,0.203859,0.226634,0.249492,0.280769,0.288487,0.297031,0.311581,0.327139,
		0.355252,0.378235,0.0430912,0.0936794,0.118107,0.154825,0.171654,0.192803,0.221002,0.23059,
		0.245322,0.257573,0.269542,0.285429,0.298099,0.310291,0.343752,0.372948,0.403576,0.0540567,
		0.0951568,0.122018,0.133266,0.152758,0.175133,0.19687,0.212454,0.220872,0.235354,0.249563,
		0.254713,0.262467,0.287776,0.303007,0.343136,0.389277,0.0587287,0.0846069,0.121343,0.137246,
		0.149819,0.166565,0.185433,0.200631,0.218634,0.225265,0.244027,0.256073,0.265865,0.300756,
		0.339348,0.389875,0.412834,0.0436849,0.0602848,0.0971883,0.121213,0.136553,0.159103,0.181875,
		0.203753,0.215924,0.22888,0.255348,0.270442,0.278257,0.308987,0.333222,0.349902,0.374102,
		0.0561145,0.0793565,0.103113,0.124414,0.152144,0.156687,0.17095,0.199606,0.219787,0.234407,
		0.249603,0.261792,0.284823,0.312216,0.325437,0.335484,0.374705,0.058914,0.101451,0.128956,
		0.141738,0.168254,0.201891,0.214053,0.224111,0.238482,0.248541,0.288767,0.300093,0.309768,
		0.324322,0.340811,0.351974,0.41013,0.069684,0.117674,0.127092,0.133463,0.154886,0.165166,
		0.191253,0.212156,0.223401,0.241783,0.269149,0.284528,0.293538,0.321312,0.343414,0.371537,
		0.386436,0.0410111,0.0770108,0.10837,0.141674,0.153368,0.172721,0.182769,0.192733,0.210766,
		0.242893,0.248462,0.259458,0.277776,0.300959,0.327088,0.343045,0.397832,0.075746,0.112508,
		0.152994,0.169419,0.183271,0.192182,0.200012,0.229489,0.239542,0.26234,0.278087,0.285547,
		0.299532,0.316097,0.325052,0.342937,0.3798,0.0568505,0.0908555,0.106362,0.138351,0.153333,
		0.169967,0.209412,0.231119,0.246289,0.261582,0.278511,0.290448,0.305389,0.315762,0.370019,
		0.385581,0.416217,0.0469136,0.0746141,0.0952416,0.135664,0.140665,0.178766,0.204778,0.223165,
		0.240898,0.262567,0.274251,0.291377,0.296859,0.308792,0.323376,0.348942,0.387432,0.0576877,
		0.0957264,0.127734,0.132302,0.153498,0.174052,0.18605,0.2021,0.219622,0.254438,0.28521,
		0.293709,0.313651,0.329152,0.355085,0.377538,0.393388,0.0290085,0.0686477,0.0904024,0.117256,
		0.137311,0.167978,0.184821,0.199111,0.231361,0.242537,0.255463,0.272599,0.28407,0.297104,
		0.325695,0.343046,0.370224,0.0397635,0.077511,0.120186,0.132964,0.156778,0.170948,0.19184,
		0.208923,0.225255,0.230447,0.243789,0.257981,0.27377,0.288198,0.317917,0.325527,0.340442,
		0.0825434,0.107836,0.125872,0.157662,0.161542,0.18898,0.199796,0.212102,0.224792,0.238337,
		0.250659,0.264685,0.286394,0.299562,0.314323,0.349633,0.36823,0.048231,0.077672,0.115867,
		0.13166,0.14546,0.174574,0.190212,0.21876,0.231986,0.24791,0.262851,0.283987,0.303624,
		0.340752,0.354301,0.374029,0.424759,0.0733051,0.0980453,0.113196,0.137149,0.143607,0.153135,
		0.189788,0.203928,0.242172,0.257076,0.276031,0.284298,0.301236,0.321448,0.330381,0.383087,
		0.430633,0.0432959,0.0926857,0.115719,0.133813,0.162835,0.175908,0.188539,0.195521,0.223389,
		0.243781,0.2507,0.276083,0.287669,0.312577,0.33307,0.3538,0.394453,0.056986,0.105666,
		0.134688,0.155021,0.163682,0.186907,0.193853,0.206537,0.226872,0.234933,0.255932,0.267824,
		0.290616,0.318909,0.330731,0.362296,0.37305,0.0434465,0.0666833,0.0811692,0.0998344,0.129987,
		0.150075,0.178549,0.202285,0.228999,0.247532,0.271072,0.279786,0.291075,0.327971,0.340329,
		0.365974,0.421268,0.0552682,0.067883,0.0927511,0.108818,0.157368,0.173847,0.192659,0.217133,
		0.228796,0.236,0.259053,0.274255,0.291243,0.309449,0.334397,0.347851,0.385243,0.0463272,
		0.0969957,0.142843,0.167481,0.195948,0.204154,0.214556,0.226546,0.240672,0.257632,0.275752,
		0.294313,0.307846,0.327886,0.345863,0.376465,0.39074,0.0516033,0.0946871,0.13139,0.146485,
		0.161933,0.183074,0.198377,0.218142,0.23701,0.248668,0.258273,0.278289,0.303549,0.328798,
		0.340406,0.376907,0.402387,0.0632555,0.0802736,0.100208,0.122046,0.134417,0.155444,0.172287,
		0.183847,0.210055,0.237385,0.245941,0.268401,0.281511,0.295032,0.306887,0.325572,0.346826,
		0.0457244,0.105376,0.130448,0.158144,0.185897,0.190812,0.199931,0.209953,0.22638,0.235898,
		0.251294,0.256923,0.277295,0.284508,0.296012,0.334649,0.368428,0.058958,0.0707226,0.0862562,
		0.111497,0.129428,0.146319,0.156663,0.170095,0.211333,0.243429,0.250409,0.261702,0.27298,
		0.304688,0.322848,0.365077,0.403764,0.0549889,0.0756352,0.0895365,0.107783,0.136396,0.157928,
		0.176623,0.18712,0.195465,0.218391,0.231926,0.242053,0.253855,0.303037,0.328417,0.337456,
		0.356178,0.0631029,0.0957855,0.114255,0.138061,0.14363,0.160491,0.180536,0.200702,0.216205,
		0.236662,0.259601,0.276716,0.298423,0.320526,0.344438,0.363771,0.378743,0.0319558,0.0591898,
		0.0961042,0.124744,0.149179,0.172466,0.182931,0.207354,0.228077,0.245264,0.268113,0.280612,
		0.299027,0.316322,0.342717,0.374282,0.388452,0.0456138,0.0764367,0.114592,0.131464,0.147096,
		0.163055,0.194683,0.206928,0.229464,0.238752,0.25625,0.276758,0.295696,0.301097,0.32061,
		0.35498,0.405465,0.0728413,0.0943603,0.112523,0.130388,0.154447,0.170106,0.197776,0.209261,
		0.229249,0.247555,0.272033,0.289901,0.303806,0.317909,0.329264,0.349262,0.39963,0.0335719,
		0.065431,0.0932403,0.11392,0.134673,0.145522,0.150511,0.194003,0.215533,0.225058,0.248293,
		0.265224,0.272617,0.297945,0.319411,0.334534,0.379182,0.0347817,0.0550451,0.0936592,0.124887,
		0.155997,0.163167,0.187212,0.207174,0.23939,0.259271,0.277748,0.290344,0.302113,0.329596,
		0.363212,0.392761,0.418568,0.0404042,0.0798069,0.0920877,0.130258,0.143753,0.159461,0.184556,
		0.197297,0.208491,0.226119,0.249759,0.273765,0.295871,0.318998,0.338154,0.369596,0.407124,
		0.0566648,0.0976747,0.108814,0.129526,0.149037,0.175067,0.183224,0.201681,0.241032,0.258876,
		0.283003,0.298701,0.324869,0.34192,0.365862,0.374398,0.403066,0.0425753,0.0614734,0.11446,
		0.128529,0.155394,0.171076,0.194635,0.203035,0.208417,0.236862,0.267307,0.27649,0.29215,
		0.304869,0.322211,0.350851,0.380552,0.0284301,0.0667816,0.105257,0.118228,0.142394,0.16424,
		0.195122,0.212529,0.226995,0.24037,0.264302,0.284574,0.29785,0.317419,0.337048,0.351453,
		0.39342,0.0531317,0.0656493,0.0979867,0.107426,0.152119,0.16687,0.194771,0.214932,0.221505,
		0.254021,0.262993,0.279104,0.298404,0.329369,0.36744,0.397281,0.420932,0.0952979,0.112101,
		0.132949,0.162404,0.175205,0.193996,0.211454,0.21751,0.233648,0.242395,0.269821,0.280955,
		0.291114,0.3015,0.338511,0.370296,0.384848,0.0655936,0.0830657,0.102511,0.137811,0.153873,
		0.167359,0.190081,0.218844,0.231135,0.247114,0.277936,0.302853,0.316698,0.337757,0.360282,
		0.387897,0.408031,0.0696748,0.092167,0.125314,0.146934,0.164461,0.179485,0.200123,0.228424,
		0.250489,0.264964,0.279722,0.296576,0.309304,0.33332,0.337939,0.359677,0.374982,0.0575944,
		0.0892923,0.116085,0.150898,0.164038,0.171576,0.1948,0.220909,0.230884,0.26506,0.279729,
		0.289337,0.299789,0.317562,0.33621,0.370949,0.424944,0.0374426,0.0814392,0.102323,0.113739,
		0.126475,0.140023,0.153616,0.16788,0.191583,0.198805,0.227061,0.246482,0.275369,0.286495,
		0.351043,0.389172,0.434464,0.0474874,0.0888205,0.118908,0.157975,0.188603,0.196854,0.203699,
		0.223392,0.230675,0.243109,0.268619,0.28404,0.300557,0.326597,0.361312,0.367757,0.394013,
		0.0672828,0.0745661,0.092963,0.111597,0.125533,0.144653,0.155148,0.182625,0.201131,0.233315,
		0.241404,0.263513,0.30615,0.328256,0.336288,0.367466,0.373472,0.0713053,0.10961,0.126978,
		0.144005,0.157493,0.164447,0.17508,0.19235,0.210703,0.247828,0.256489,0.28283,0.300482,
		0.327151,0.340344,0.352466,0.39147,0.0879174,0.112189,0.135176,0.144561,0.156346,0.183596,
		0.193246,0.209317,0.225081,0.24034,0.260008,0.273831,0.299575,0.327122,0.361353,0.38525,
		0.412472,0.0446764,0.0705274,0.0967181,0.108961,0.136944,0.146145,0.154098,0.1756,0.191598,
		0.203848,0.223812,0.235063,0.289393,0.296176,0.317556,0.355631,0.38292,0.0686833,0.102484,
		0.123527,0.14176,0.153317,0.165279,0.176563,0.202243,0.225123,0.249096,0.261093,0.280871,
		0.306437,0.325101,0.334561,0.351153,0.382699,0.0386355,0.0615877,0.105781,0.127815,0.173211,
		0.198744,0.211987,0.221692,0.235616,0.251412,0.271586,0.279534,0.288084,0.313007,0.321035,
		0.348882,0.364058,0.0639119,0.0980269,0.115029,0.135102,0.149924,0.180839,0.190455,0.207934,
		0.223521,0.239124,0.259606,0.269464,0.288757,0.306717,0.320385,0.346998,0.373416,0.0367776,
		0.0759555,0.10581,0.131481,0.141885,0.169345,0.19093,0.216963,0.23513,0.247151,0.270766,
		0.282122,0.2962,0.313961,0.330749,0.373204,0.414288,0.0361996,0.0790776,0.0967399,0.115986,
		0.135969,0.173123,0.184563,0.202301,0.239789,0.247554,0.265935,0.283096,0.289404,0.310422,
		0.339714,0.362147,0.399068,0.0645223,0.104762,0.136933,0.15527,0.169919,0.194018,0.217203,
		0.228558,0.244311,0.253957,0.262353,0.277148,0.288172,0.306851,0.327866,0.361736,0.377056,
		0.019434,0.0611104,0.0769851,0.112383,0.133212,0.156786,0.168695,0.188895,0.211367,0.229676,
		0.242844,0.268063,0.290517,0.31524,0.343051,0.350508,0.390929,0.0375679,0.0568989,0.112008,
		0.126088,0.142399,0.16574,0.195719,0.209055,0.223734,0.247832,0.260022,0.286349,0.298334,
		0.315013,0.331174,0.360408,0.402616,0.0438591,0.058596,0.085943,0.102944,0.135411,0.167502,
		0.194606,0.213601,0.243318,0.268107,0.278028,0.298248,0.313023,0.338607,0.36906,0.393389,
		0.410219,0.0279155,0.0543581,0.101675,0.125117,0.161619,0.176894,0.197527,0.207601,0.218565,
		0.234052,0.247047,0.266703,0.279215,0.309739,0.330682,0.346121,0.372843,0.0311365,0.0665587,
		0.0913395,0.109447,0.12486,0.153329,0.172252,0.188736,0.209259,0.245096,0.254932,0.269417,
		0.279777,0.291212,0.330218,0.35204,0.371515,0.0541091,0.0781329,0.11208,0.128634,0.147338,
		0.174048,0.195388,0.211804,0.222972,0.237531,0.250581,0.267409,0.279349,0.294047,0.325958,
		0.355269,0.381029,0.0568813,0.0811195,0.10948,0.123813,0.136097,0.161443,0.184892,0.20637,
		0.220622,0.235694,0.245465,0.265285,0.272836,0.30044,0.310084,0.323531,0.338605,0.0635852,
		0.128975,0.138236,0.152637,0.163128,0.169662,0.181001,0.204523,0.219178,0.232979,0.248811,
		0.283738,0.314527,0.31815,0.332819,0.395547,0.436076,0.0606917,0.0872101,0.118522,0.1437,
		0.16821,0.189876,0.206734,0.226129,0.242668,0.247195,0.261269,0.269441,0.287357,0.298631,
		0.311916,0.344182,0.382126,0.0516573,0.0962368,0.107738,0.125946,0.137595,0.178185,0.197577,
		0.207417,0.21841,0.229219,0.251894,0.287382,0.311368,0.323527,0.338282,0.350829,0.382475,
		0.0429823,0.116589,0.13537,0.161218,0.181873,0.19913,0.215,0.231286,0.23894,0.247342,
		0.255817,0.276557,0.295004,0.336326,0.353093,0.377909,0.40763,0.0437774,0.0873471,0.112953,
		0.129996,0.149958,0.170496,0.181269,0.205413,0.212885,0.222347,0.245556,0.259638,0.289259,
		0.300502,0.345816,0.354053,0.375028,0.0182351,0.0810483,0.0946979,0.105476,0.130595,0.16602,
		0.1859,0.199801,0.214786,0.232254,0.24555,0.262998,0.277686,0.307597,0.320662,0.341339,
		0.358361,0.0238025,0.083126,0.0908245,0.112428,0.125263,0.154985,0.166171,0.185487,0.201225,
		0.233377,0.245572,0.267736,0.290685,0.306254,0.324263,0.351315,0.399721,0.0254217,0.0606709,
		0.100385,0.139113,0.156797,0.176027,0.190301,0.201001,0.224683,0.236598,0.252712,0.260757,
		0.282465,0.305395,0.316905,0.343241,0.370027,0.0656468,0.0758587,0.12022,0.144788,0.14787,
		0.163892,0.186803,0.195486,0.218744,0.237597,0.245772,0.261066,0.283836,0.297696,0.318725,
		0.329331,0.364043,0.0697679,0.100397,0.115345,0.124396,0.160097,0.173061,0.195154,0.204569,
		0.22617,0.249143,0.25764,0.288327,0.303659,0.344634,0.358289,0.374295,0.400307,0.0413221,
		0.0671446,0.0926959,0.119182,0.1411,0.167246,0.184673,0.20159,0.210919,0.2335,0.254968,
		0.271204,0.287514,0.314985,0.328676,0.341006,0.359842,0.0405382,0.0730318,0.111152,0.143112,
		0.150991,0.169206,0.182656,0.19843,0.215006,0.222811,0.233877,0.240588,0.252617,0.274436,
		0.300571,0.320439,0.356682,0.034198,0.092681,0.122671,0.138406,0.153367,0.174449,0.194269,
		0.221541,0.234801,0.259002,0.27679,0.28525,0.299955,0.321325,0.347045,0.36452,0.385617,
		0.032949,0.0762426,0.120847,0.136494,0.159318,0.178064,0.219645,0.235987,0.254968,0.277571,
		0.289337,0.305278,0.324466,0.337103,0.352489,0.38003,0.421366,0.024951,0.0729372,0.103741,
		0.13661,0.160528,0.168567,0.178127,0.212368,0.238717,0.253089,0.262281,0.282942,0.300474,
		0.308604,0.342449,0.356115,0.372193,0.0438464,0.0719709,0.109215,0.13145,0.14956,0.164209,
		0.185244,0.204919,0.223854,0.250523,0.265019,0.298928,0.305712,0.325412,0.348072,0.370968,
		0.391405,0.0825069,0.109226,0.126773,0.136214,0.149318,0.17037,0.189122,0.202935,0.219778,
		0.23138,0.247958,0.283153,0.31312,0.328724,0.359866,0.381204,0.422482,0.0540477,0.0799437,
		0.0995883,0.115345,0.144732,0.172892,0.183303,0.213113,0.236045,0.252728,0.261907,0.27848,
		0.305278,0.341025,0.35196,0.368432,0.413581,0.0429749,0.0649845,0.0923199,0.134886,0.149507,
		0.162055,0.181766,0.228054,0.246436,0.259016,0.292552,0.298096,0.311936,0.338366,0.353093,
		0.368549,0.412401,0.0719394,0.0800354,0.101689,0.127091,0.141908,0.182823,0.205142,0.211237,
		0.223211,0.245174,0.249111,0.263141,0.275812,0.292441,0.341732,0.364242,0.387486,0.0588266,
		0.0896343,0.12081,0.151845,0.179223,0.199778,0.212554,0.22617,0.249877,0.257634,0.275123,
		0.288804,0.303532,0.326695,0.33826,0.371608,0.399948,0.0408899,0.0782837,0.0942494,0.115202,
		0.134635,0.164333,0.178285,0.186248,0.214121,0.231893,0.23893,0.251481,0.282123,0.317101,
		0.339093,0.353034,0.419429,0.0247515,0.0510598,0.0905373,0.10723,0.122324,0.134969,0.178449,
		0.184837,0.203502,0.229405,0.245213,0.268701,0.296626,0.321211,0.353104,0.36781,0.395715,
		0.0406598,0.0788507,0.115703,0.126865,0.148123,0.164925,0.184304,0.214826,0.22361,0.237934,
		0.250989,0.256843,0.264579,0.277264,0.298339,0.328892,0.34901,0.0429823,0.086183,0.143102,
		0.155269,0.179977,0.20849,0.233577,0.241447,0.254625,0.266579,0.271759,0.294998,0.311178,
		0.328886,0.338601,0.361918,0.385212,0.0297326,0.0540998,0.083523,0.109027,0.129129,0.153852,
		0.1869,0.196998,0.21586,0.228344,0.247754,0.263283,0.277571,0.286302,0.291752,0.322435,
		0.365943,0.0539442,0.0897239,0.104567,0.125302,0.146155,0.161041,0.173448,0.222517,0.239487,
		0.252943,0.263264,0.281516,0.304864,0.309293,0.327143,0.349615,0.3876,0.0576386,0.0896137,
		0.11196,0.119349,0.143124,0.156404,0.176243,0.201001,0.223559,0.255407,0.268326,0.276349,
		0.300582,0.309023,0.344806,0.364032,0.397675,0.0417095,0.0634765,0.08016,0.108358,0.123127,
		0.140387,0.147878,0.164312,0.178082,0.228119,0.238212,0.264436,0.278355,0.302501,0.324545,
		0.355423,0.378782,0.0254217,0.0739904,0.102383,0.121792,0.128919,0.140626,0.159228,0.184304,
		0.220421,0.233783,0.252764,0.265752,0.275974,0.296616,0.316556,0.342613,0.389735,0.0342126,
		0.0626945,0.0760776,0.128034,0.142805,0.163876,0.180103,0.18854,0.210951,0.229253,0.23628,
		0.242857,0.254157,0.2936,0.313397,0.340697,0.362404,0.0490687,0.0791758,0.0922137,0.108199,
		0.129416,0.143112,0.152048,0.182894,0.19913,0.216031,0.225534,0.275446,0.299995,0.313188,
		0.347549,0.359866,0.375305,0.0723149,0.106286,0.120741,0.149613,0.158578,0.190934,0.214815,
		0.231095,0.234739,0.249811,0.260344,0.293112,0.307037,0.320104,0.328273,0.347247,0.39001,
		0.0392409,0.0777691,0.124351,0.134304,0.152325,0.178082,0.190839,0.213003,0.228026,0.251366,
		0.255753,0.260753,0.27665,0.295235,0.307923,0.321158,0.338041,0.0211083,0.0472769,0.099423,
		0.118591,0.145187,0.171471,0.210919,0.220592,0.244606,0.258733,0.274242,0.289017,0.301421,
		0.319402,0.325585,0.354953,0.388692,0.0529063,0.0886211,0.107023,0.121299,0.150169,0.163345,
		0.188686,0.208569,0.222327,0.227948,0.243561,0.259059,0.281091,0.303144,0.324196,0.359329,
		0.3711,0.0325257,0.0673999,0.100321,0.128919,0.132222,0.158205,0.173906,0.219239,0.235987,
		0.244658,0.271109,0.291777,0.320241,0.347686,0.3575,0.369995,0.398236,0.0481098,0.0832122,
		0.0966171,0.102572,0.118092,0.149805,0.16948,0.197187,0.206881,0.218021,0.237999,0.26992,
		0.28863,0.306098,0.32708,0.337103,0.40836,0.0349807,0.0519815,0.0788119,0.103741,0.117926,
		0.154594,0.177361,0.201774,0.219277,0.242651,0.25398,0.268596,0.283823,0.310947,0.339371,
		0.351005,0.414781,0.0605683,0.0946979,0.120194,0.140154,0.164058,0.199,0.217363,0.23727,
		0.247754,0.258088,0.273642,0.285289,0.297307,0.321723,0.345586,0.370307,0.38574,0.0290423,
		0.0695697,0.0975931,0.119177,0.153677,0.167246,0.195129,0.214121,0.240195,0.249569,0.260409,
		0.288765,0.303927,0.318506,0.335048,0.352125,0.368061,0.0416555,0.0723032,0.117411,0.139442,
		0.149805,0.16255,0.174698,0.190837,0.215152,0.232367,0.243511,0.252533,0.269066,0.290167,
		0.309034,0.331857,0.351417,0.0210896,0.0671188,0.0867045,0.132566,0.155297,0.171647,0.191729,
		0.222284,0.230467,0.245539,0.275812,0.291409,0.302123,0.310163,0.339788,0.368678,0.383497,
		0.0460161,0.108245,0.125355,0.145685,0.158205,0.180017,0.196589,0.217611,0.229142,0.242125,
		0.250036,0.279375,0.285545,0.302115,0.319679,0.344715,0.391573,0.0522273,0.0644326,0.0878344,
		0.114636,0.137458,0.163777,0.18116,0.205413,0.221581,0.240666,0.26264,0.273508,0.288273,
		0.300582,0.32137,0.330416,0.374469,0.0688875,0.091026,0.138104,0.178345,0.195129,0.206663,
		0.228208,0.234449,0.249421,0.261903,0.274724,0.284741,0.296236,0.308743,0.323814,0.34168,
		0.359488,0.0432181,0.072469,0.0942494,0.128388,0.134969,0.161952,0.193042,0.206527,0.226769,
		0.241584,0.248871,0.265148,0.28629,0.317601,0.331733,0.365457,0.390148,0.0468986,0.0739788,
		0.0975462,0.119229,0.144035,0.160365,0.194193,0.213419,0.223368,0.247286,0.256271,0.276607,
		0.292839,0.304239,0.308167,0.323349,0.343683,0.0400464,0.0872918,0.105679,0.131562,0.176364,
		0.187441,0.205481,0.223778,0.234739,0.255657,0.264747,0.275425,0.290036,0.311169,0.338878,
		0.365686,0.384544,0.0516011,0.0841388,0.121743,0.140507,0.159281,0.169571,0.176513,0.190916,
		0.208674,0.228787,0.240078,0.25593,0.271894,0.294544,0.306018,0.340094,0.371927,0.0573174,
		0.0893611,0.121898,0.134304,0.160515,0.181114,0.19993,0.213414,0.242917,0.259057,0.269807,
		0.286408,0.314549,0.329938,0.363807,0.382041,0.402921,0.060685,0.099745,0.128034,0.138859,
		0.151128,0.172454,0.181816,0.204179,0.222811,0.253464,0.265714,0.271894,0.290439,0.310983,
		0.324915,0.337323,0.363962,0.0494473,0.0840409,0.118664,0.141908,0.167875,0.213152,0.22771,
		0.245556,0.263719,0.279533,0.299955,0.313696,0.324772,0.339371,0.357185,0.37155,0.403099,
		0.0187305,0.0896137,0.12895,0.140149,0.160515,0.183142,0.198756,0.216591,0.229377,0.252317,
		0.260753,0.269787,0.285293,0.302927,0.327982,0.35332,0.372877,0.060685,0.0765356,0.0929612,
		0.106286,0.128017,0.142014,0.152775,0.195445,0.215,0.239536,0.253903,0.258472,0.281393,
		0.290256,0.311662,0.329016,0.368441,0.0519857,0.0830867,0.103959,0.13018,0.153039,0.162921,
		0.20159,0.214424,0.229596,0.239761,0.255993,0.268445,0.28385,0.298435,0.313749,0.335703,
		0.355586,0.0360689,0.0880325,0.104567,0.12538,0.136162,0.181269,0.193558,0.222767,0.238223,
		0.246762,0.260265,0.278542,0.288163,0.300571,0.313063,0.334684,0.378436,0.0257914,0.0665855,
		0.091391,0.1399,0.17037,0.191809,0.201439,0.227611,0.236288,0.253083,0.26634,0.278751,
		0.289561,0.300151,0.314464,0.346928,0.380076,0.0267962,0.0734752,0.112609,0.139658,0.157978,
		0.172724,0.206153,0.216329,0.226829,0.244244,0.256716,0.269416,0.283669,0.304598,0.319914,
		0.337525,0.357094,0.0601116,0.0794563,0.0953404,0.110513,0.129084,0.145106,0.161302,0.185252,
		0.215238,0.239516,0.248008,0.27365,0.290439,0.310947,0.327486,0.355117,0.37326,0.0357931,
		0.0687856,0.119622,0.140164,0.154985,0.167103,0.178206,0.198603,0.204774,0.223863,0.229377,
		0.246436,0.256843,0.270462,0.283669,0.323557,0.371663,0.0763898,0.0872851,0.115966,0.135279,
		0.158895,0.183784,0.208755,0.219194,0.233886,0.247292,0.253638,0.277719,0.303529,0.325219,
		0.350728,0.385261,0.409837,0.0291784,0.0623985,0.0902539,0.125302,0.143008,0.164465,0.180279,
		0.201227,0.223693,0.249637,0.265334,0.29746,0.31511,0.325035,0.342783,0.389107,0.428717,
		0.034479,0.0689162,0.0808138,0.10312,0.12895,0.143426,0.17849,0.187099,0.196998,0.212967,
		0.232237,0.260862,0.285187,0.317045,0.341451,0.366044,0.384483,0.0294083,0.0689162,0.0745001,
		0.100233,0.129658,0.162418,0.183342,0.208173,0.233903,0.261836,0.278241,0.295018,0.308237,
		0.32612,0.355423,0.377994,0.394346,0.0667511,0.0945766,0.131892,0.151258,0.160854,0.180758,
		0.198861,0.206438,0.228344,0.239398,0.255364,0.270462,0.290256,0.316288,0.324077,0.344457,
		0.381615,0.0238723,0.0467784,0.0762426,0.102572,0.144796,0.161075,0.17516,0.201822,0.21962,
		0.240588,0.277114,0.298435,0.31475,0.330296,0.348834,0.374607,0.397804,0.0490392,0.0745001,
		0.0965316,0.117427,0.13655,0.163576,0.206881,0.225486,0.254142,0.268643,0.286144,0.307506,
		0.323557,0.338041,0.355117,0.390582,0.427046,0.0739599,0.0857164,0.118591,0.130489,0.155453,
		0.171156,0.200052,0.21066,0.216782,0.221262,0.231992,0.26336,0.268483,0.286689,0.308169,
		0.326993,0.375398,0.0509804,0.0724875,0.106456,0.124534,0.13805,0.17524,0.198384,0.224798,
		0.243483,0.254252,0.275238,0.285407,0.30601,0.329016,0.34971,0.411144,0.432352,0.0649271,
		0.0870339,0.105582,0.126686,0.161622,0.180779,0.199843,0.232237,0.24525,0.250544,0.273206,
		0.288595,0.304864,0.332626,0.353163,0.386162,0.401151,0.0527461,0.0680186,0.0849551,0.111683,
		0.123071,0.142014,0.158762,0.172756,0.190301,0.210425,0.221826,0.2464,0.26541,0.29173,
		0.323812,0.346992,0.370356,0.0371139,0.0745062,0.0836585,0.108602,0.138046,0.152325,0.172471,
		0.207484,0.220601,0.235251,0.260352,0.275622,0.289259,0.312893,0.343983,0.377589,0.397997,
		0.0191161,0.0809807,0.113651,0.146795,0.160939,0.182366,0.199615,0.220172,0.236618,0.24817,
		0.264153,0.285001,0.297566,0.310979,0.322435,0.337145,0.370027,0.0605397,0.0867582,0.108266,
		0.120886,0.147652,0.170453,0.186127,0.197677,0.221581,0.237506,0.253555,0.27881,0.295604,
		0.301346,0.323195,0.353004,0.380305,0.0468986,0.0953353,0.118519,0.138525,0.155171,0.1859,
		0.20991,0.231504,0.238402,0.245657,0.25764,0.27706,0.29994,0.322358,0.332654,0.367859,
		0.389632,0.0294843,0.0701895,0.0770087,0.0944053,0.124263,0.158446,0.168854,0.184525,0.216294,
		0.239285,0.25657,0.275273,0.282131,0.31137,0.331242,0.345663,0.375289,0.0432181,0.096997,
		0.119332,0.129252,0.150753,0.160801,0.169078,0.182372,0.200918,0.208399,0.235204,0.260613,
		0.298814,0.314868,0.327395,0.371565,0.424312,0.0573174,0.0892526,0.0959908,0.124263,0.158719,
		0.172869,0.206022,0.234749,0.257262,0.26992,0.289116,0.29931,0.31323,0.334241,0.342783,
		0.374832,0.441573,0.0287714,0.0624478,0.0723276,0.0980249,0.110737,0.138872,0.183152,0.194421,
		0.207161,0.219645,0.244244,0.255833,0.274424,0.302115,0.325354,0.342624,0.357891,0.0405382,
		0.0697778,0.0887744,0.108371,0.12224,0.148456,0.168489,0.180871,0.199615,0.22727,0.236971,
		0.25218,0.264436,0.297474,0.308389,0.348082,0.357414,0.0452574,0.0796688,0.103979,0.109754,
		0.127034,0.145477,0.162435,0.189112,0.210954,0.223863,0.236978,0.26634,0.279053,0.29881,
		0.344697,0.365205,0.395187,0.0468884,0.0734818,0.102219,0.131344,0.154176,0.161882,0.180871,
		0.191633,0.210361,0.214434,0.244764,0.250523,0.265392,0.282835,0.296982,0.327395,0.362187,
		0.0337862,0.0569165,0.0719394,0.0902539,0.148102,0.178634,0.19414,0.213983,0.223502,0.236974,
		0.248561,0.267171,0.283108,0.291275,0.310351,0.349677,0.374153,0.0429749,0.0770087,0.0979643,
		0.113906,0.134233,0.162418,0.178094,0.18502,0.212445,0.231308,0.253089,0.293103,0.304335,
		0.321223,0.343683,0.369963,0.390484,0.0294843,0.0536117,0.0768307,0.14098,0.155171,0.17379,
		0.183778,0.207812,0.218438,0.227737,0.242969,0.258054,0.276302,0.304549,0.328843,0.348086,
		0.364699,0.0544003,0.0748759,0.0882345,0.108995,0.143438,0.170193,0.203452,0.229856,0.24058,
		0.248297,0.26506,0.276557,0.312623,0.326061,0.341535,0.368443,0.390087,0.0426603,0.0810284,
		0.0948114,0.129712,0.147516,0.163175,0.186248,0.203733,0.228016,0.239398,0.255347,0.27937,
		0.292298,0.303611,0.315859,0.328676,0.340651,0.0249337,0.0600723,0.1067,0.121743,0.166112,
		0.184196,0.201486,0.220325,0.234601,0.244334,0.256964,0.268258,0.286408,0.301215,0.325355,
		0.341025,0.382311,0.0193913,0.0784439,0.09719,0.12322,0.158875,0.168713,0.183924,0.204179,
		0.218973,0.227241,0.237815,0.249637,0.263512,0.279654,0.289616,0.322278,0.347916,0.0533644,
		0.0695697,0.0973604,0.133046,0.15518,0.165633,0.18449,0.200784,0.225534,0.23791,0.248275,
		0.2612,0.282597,0.309733,0.323349,0.340651,0.37155,0.0531367,0.0802278,0.107643,0.120238,
		0.153548,0.183137,0.207241,0.216725,0.231308,0.239665,0.265829,0.279995,0.297266,0.308604,
		0.334854,0.351158,0.375898,0.042439,0.0678254,0.105002,0.130842,0.141349,0.173906,0.193594,
		0.211665,0.224356,0.242569,0.263141,0.278241,0.291406,0.306712,0.326149,0.371733,0.411532,
		0.0502504,0.0827703,0.114938,0.13805,0.172112,0.204713,0.22892,0.238885,0.256904,0.278276,
		0.282465,0.299995,0.315687,0.339021,0.354953,0.37644,0.402809,0.0366766,0.0461328,0.080747,
		0.0875586,0.115282,0.147233,0.17547,0.206022,0.228739,0.250989,0.268561,0.288408,0.309832,
		0.340794,0.367288,0.387134,0.400426,0.0729372,0.0921201,0.113503,0.128606,0.140485,0.153661,
		0.181114,0.198603,0.216262,0.228639,0.243483,0.252638,0.262598,0.284863,0.314868,0.331216,
		0.342597,0.0390168,0.091058,0.110513,0.157178,0.165912,0.184545,0.197216,0.215874,0.225254,
		0.233302,0.252093,0.263674,0.280868,0.290036,0.321223,0.365311,0.380076,0.0509804,0.0998073,
		0.121113,0.138096,0.158895,0.168489,0.197739,0.216553,0.232878,0.235647,0.24645,0.265334,
		0.302501,0.318863,0.332871,0.353004,0.365311,0.0257914,0.0970678,0.119571,0.14743,0.157547,
		0.175262,0.189854,0.19993,0.214442,0.226679,0.242173,0.263468,0.279749,0.292298,0.324121,
		0.345613,0.364242,0.0160165,0.0582143,0.0760776,0.102967,0.131344,0.147701,0.168407,0.199063,
		0.218294,0.22892,0.253083,0.265747,0.288595,0.307789,0.332739,0.361959,0.370454,0.0575998,
		0.0791604,0.104538,0.127093,0.139824,0.161894,0.175744,0.227832,0.240826,0.263662,0.271755,
		0.28523,0.291728,0.32449,0.341339,0.359511,0.390479,0.0289331,0.0558541,0.0803647,0.111065,
		0.12456,0.150629,0.174511,0.183142,0.212572,0.231147,0.246917,0.270959,0.289687,0.304892,
		0.323296,0.366675,0.400116,0.0333046,0.0601116,0.0887744,0.113103,0.140164,0.166392,0.186773,
		0.20422,0.224274,0.232744,0.245657,0.2642,0.292881,0.309293,0.329331,0.348656,0.364618
	};
	
	
	//for every test point:
	unsigned counter = 0;
	for ( size_t pts = 0 ; pts != test_points.size() ; ++pts )
	{
		//for every distance class:
		for ( size_t dist = 0 ; dist != classes_to_compute_knn.size() ; ++dist )
		{
			BOOST_CHECK( fabs( (*classes_to_compute_knn[dist])( test_points[pts]) - result[counter] ) <= 5e-05 );			
			++counter;
		}
	}		
}	




























BOOST_AUTO_TEST_CASE(Distance_to_k_th_nearest_neighbor_periodic_domain_4d_brute_force)
{
	//first we test the brute force algorithm:

	typedef Gudhi::Topological_inference_with_cubical_complexes::Euclidan_distance_squared Euclidean_distance_squared;		
	typedef Gudhi::Topological_inference_with_cubical_complexes::periodic_domain_distance<Euclidean_distance_squared>
	periodic_Euclidean_distance_squared;
	
	
	std::vector< std::pair< double , double > > coordinates_of_grid(4);
	coordinates_of_grid[0] = std::pair< double,double >(0,1);
	coordinates_of_grid[1] = std::pair< double,double >(0,1);
	coordinates_of_grid[2] = std::pair< double,double >(0,1);
	coordinates_of_grid[3] = std::pair< double,double >(0,1);
	
	Euclidean_distance_squared eu;
	periodic_Euclidean_distance_squared period_eu( coordinates_of_grid , eu );
	
	
	//150 random points from R^4.
	std::vector< std::vector<double> > point_cloud = 
	{
	{0.2161353165,0.4067554909,0.0665765158,0.870283464},{0.7782972094,0.6987265884,0.7606270006,0.7958696354},{0.9662634397,0.1367719397,0.8092096795,0.716037544},
	{0.275288115,0.1794606529,0.8734926581,0.1608733959},{0.3935741177,0.5092781719,0.5718593814,0.355841625},{0.9207663401,0.5319217448,0.258124053,0.1236327693},
	{0.8471451914,0.4675919802,0.4061658739,0.9562745956},{0.9294778437,0.1459124812,0.5308060655,0.4604051244},{0.3011312045,0.4122534681,0.337005893,0.9334948768},
	{0.0638816946,0.1186342924,0.7877057337,0.6814819723},{0.3170833248,0.5444929474,0.15133534,0.7217780657},{0.4718955676,0.0654609054,0.126377013,0.9821652356},
	{0.7545759601,0.0546513223,0.8329449217,0.6777689687},{0.5130777196,0.500694545,0.6574573508,0.7392099628},{0.5893094379,0.0494349443,0.5050142871,0.0808209586},
	{0.4940861093,0.2490698399,0.3612209754,0.7160524682},{0.7133854455,0.0851763973,0.4031818854,0.2452370371},{0.6502694073,0.6256448289,0.0068036765,0.9420920806},
	{0.2771677715,0.5365420913,0.2049486258,0.6807608032},{0.6708589897,0.2624673771,0.5366855075,0.036167918},{0.6042193428,0.3100096218,0.0774625444,0.9480968153},
	{0.8630631317,0.4616721999,0.1888999115,0.4816821855},{0.274362145,0.6687946867,0.9600374054,0.5621335753},{0.9371821445,0.8886038563,0.4903153491,0.5409146673},
	{0.1577256215,0.9017276918,0.6047877711,0.2715029821},{0.4778309672,0.3667936909,0.2353694309,0.2162692775},{0.9561227295,0.2203476694,0.8832569388,0.8253242387},
	{0.3162962156,0.2985790542,0.3671411616,0.1277638588},{0.9136061701,0.1122253346,0.5190847388,0.3615447618},{0.327857591,0.6936850683,0.2818536318,0.5899669132},
	{0.6639132625,0.6278195342,0.8459756977,0.9733827612},{0.8529426146,0.8311685929,0.5352412218,0.5575795863},{0.8291846917,0.7364060413,0.5959779969,0.8778938674},
	{0.7160940336,0.7085216588,0.2321347578,0.2757493216},{0.5377012696,0.9137963476,0.3599800142,0.7714768103},{0.5407209711,0.9124356993,0.0977215923,0.6965975054},
	{0.183090288,0.4487572513,0.3782904518,0.515378186},{0.2529464646,0.2026785347,0.7370024058,0.8433211406},{0.7773119675,0.5135940188,0.0091130023,0.5177475526},
	{0.2755686354,0.5956264965,0.5519859341,0.6927454374},{0.5180953105,0.7024236836,0.873931882,0.3592949654},{0.9604220481,0.627821804,0.0671470871,0.8844678551},
	{0.2330840409,0.819831975,0.1885772266,0.0399145884},{0.6747471029,0.4076069591,0.3782148864,0.2600041842},{0.938375168,0.8666493378,0.9016399393,0.4959568656},
	{0.4837143579,0.5254478094,0.5813259515,0.5528810667},{0.4965164606,0.2691201437,0.9339421515,0.607314809},{0.3963730407,0.0023954718,0.5331149858,0.4707511207},
	{0.0545727671,0.3826355562,0.3722823337,0.7533125849},{0.7794329945,0.6273420472,0.475966346,0.7291254057},{0.302274507,0.0816056114,0.415606502,0.727463753},
	{0.8516243717,0.5476219286,0.1300290304,0.8781453874},{0.7446079317,0.4378378205,0.212476521,0.6074005007},{0.0543138136,0.0176691897,0.7037446832,0.3270085847},
	{0.9654455474,0.6369031293,0.4479166996,0.1196174838},{0.1637513947,0.8642580179,0.2845635938,0.7220461073},{0.2128266117,0.9449356911,0.4176520605,0.591265711},
	{0.0476004935,0.286355197,0.2003003603,0.3544409689},{0.4617424291,0.34739498,0.5184426282,0.7242208878},{0.9389569708,0.1110231471,0.1960562568,0.8751947973},
	{0.1297257203,0.1082481672,0.9613296099,0.6461059342},{0.6314383652,0.863803653,0.4156785782,0.1647425564},{0.8742997588,0.7209185674,0.1300667289,0.9857225746},
	{0.4410702195,0.5473736799,0.3925960511,0.073340355},{0.2864084891,0.0529552582,0.4401347351,0.1997130518},{0.2518973188,0.1512412846,0.0166082208,0.244082951},
	{0.2608718651,0.3542067721,0.1776841788,0.5686277715},{0.0090282573,0.6212397611,0.6485787611,0.9783152305},{0.3855588885,0.0695148581,0.036164256,0.549218531},
	{0.1092119429,0.1176601963,0.9561683217,0.0636541741},{0.0621516255,0.4523660035,0.7733616829,0.5127586231},{0.6562764596,0.2179020974,0.5065915762,0.3336082227},
	{0.9608324207,0.2805596716,0.9524878454,0.4502769466},{0.4204976272,0.9747959129,0.6945173966,0.0523253432},{0.9322453286,0.9045769402,0.4068475547,0.7923901044},
	{0.0912914628,0.4110324632,0.1179300488,0.4090027225},{0.9372338173,0.0228120044,0.1099257688,0.7709252352},{0.2473172548,0.9430129114,0.9248211253,0.2755537038},
	{0.3151711877,0.3909966208,0.5856341401,0.8289346434},{0.6222296427,0.7623051493,0.7437775929,0.5864307864},{0.2865205696,0.4823481457,0.2886831271,0.16002433},
	{0.2867317577,0.3187225992,0.6152163229,0.9980220243},{0.814303647,0.6177073771,0.0959672085,0.5412626567},{0.5198940996,0.8805264027,0.5142570818,0.9733662861},
	{0.0274844745,0.193287927,0.969342147,0.9436305743},{0.9917424994,0.652558957,0.6949542321,0.9143778963},{0.158117207,0.4643659701,0.8407568894,0.79356269},
	{0.6943920886,0.5211990806,0.2490926124,0.58702025},{0.6541512765,0.5538034292,0.2002501043,0.499721779},{0.2718869599,0.667435063,0.1980533334,0.411654298},
	{0.2582200989,0.838053823,0.5876435062,0.2448927639},{0.452425757,0.8820036445,0.9339749734,0.4761730228},{0.6869763259,0.8356626607,0.4592519393,0.1507639568},
	{0.1855826536,0.1115194834,0.9367056177,0.6076889979},{0.2906259801,0.5802919739,0.2736881424,0.3915352314},{0.12288481,0.4239620366,0.0535635338,0.5429183329},
	{0.6517687093,0.0795862924,0.1951899431,0.8655872962},{0.4832760945,0.1621366784,0.7051915845,0.2026124394},{0.9915547213,0.7069917477,0.5824110755,0.843934183},
	{0.0752148305,0.6797068669,0.1407429664,0.2335030579},{0.602499583,0.2250958043,0.7167737996,0.5247004263},{0.0792115591,0.0035277302,0.9021162488,0.9836529621},
	{0.0630885197,0.6648809183,0.2074113868,0.986199873},{0.7841980574,0.054547969,0.3682132296,0.3902004855},{0.5905348728,0.179866985,0.3925852075,0.8659482158},
	{0.7789685093,0.7006078116,0.2386029409,0.5079281065},{0.8778327638,0.1721387503,0.4325646544,0.8477176272},{0.9655417092,0.8003545383,0.2018559808,0.447240192},
	{0.3959916683,0.716806635,0.3991919297,0.4022070405},{0.3342622917,0.6906662339,0.6621937822,0.1680860706},{0.6058353547,0.8673584729,0.7881982641,0.4259585107},
	{0.5841680258,0.7195332507,0.0491419374,0.1035793731},{0.4424925682,0.1699164014,0.0806449894,0.8657210511},{0.4598111401,0.7129467162,0.1730200783,0.6224096753},
	{0.9347412391,0.8530473537,0.3221795461,0.3831341839},{0.515093598,0.9728921819,0.3625707706,0.6163540559},{0.4335091764,0.9474623057,0.3618741892,0.6030456028},
	{0.6046422494,0.2836575734,0.6610661431,0.2691357303},{0.2177470808,0.5447009606,0.8826168361,0.1196069277},{0.3651480607,0.7555859217,0.0177540325,0.5050244427},
	{0.7203258723,0.3371161658,0.4509870093,0.9302381405},{0.0989933989,0.2220814053,0.050527609,0.5632122078},{0.369650627,0.5319386604,0.9319623525,0.9647439241},
	{0.8442949792,0.6522289219,0.7757588418,0.3020221314},{0.9597743074,0.558198594,0.9786667374,0.5753302046},{0.2922401265,0.3846208933,0.748712206,0.0844571677},
	{0.3025243522,0.1353433062,0.3826868525,0.1765479192},{0.7594928641,0.8685317831,0.2470071204,0.111033147},{0.9367945509,0.9874413582,0.6886330319,0.0756251651},
	{0.7691251929,0.0561589014,0.5874789746,0.006362404},{0.7453925973,0.2688428422,0.930551802,0.7137127921},{0.4403857084,0.0046569624,0.3807356679,0.8867194413},
	{0.5359529725,0.5732941926,0.6598083961,0.8664365266},{0.4278839084,0.705274327,0.7008550165,0.1929006591},{0.0876745624,0.7284803414,0.1906603929,0.3610600757},
	{0.9920528498,0.2349508412,0.8696436284,0.8481443948},{0.0071836323,0.377837914,0.4632825549,0.4406018984},{0.2740049621,0.4837526211,0.6996950118,0.2481896421},
	{0.5023594771,0.6058867699,0.5802190339,0.0271982334},{0.6278990563,0.8451684192,0.6360527,0.8441247668},{0.9253878356,0.4954766021,0.4400199461,0.7489505075},
	{0.5729722457,0.3896644767,0.3697294721,0.0447440601},{0.9223849133,0.5377308403,0.3215464076,0.7362179228},{0.9915811059,0.863216745,0.5659185979,0.913221871},
	{0.1882577084,0.9678615758,0.3086162433,0.2597987708},{0.8596413189,0.2405901959,0.8497658458,0.6602620669},{0.259099799,0.8367184943,0.1531046529,0.1532220987},
	{0.1685487258,0.7992893311,0.2069046353,0.5789962241},{0.8992731075,0.4039457724,0.5378789764,0.682111514},{0.1674966314,0.6672669824,0.7791987706,0.5535525433}
	};	
	
	std::vector< Gudhi::Topological_inference_with_cubical_complexes::Distance_to_k_th_closest_point<periodic_Euclidean_distance_squared>* >
	classes_to_compute_knn;	
	for ( size_t i = 5 ; i != 90 ; i=i+5 )
	{
		classes_to_compute_knn.push_back(
		new Gudhi::Topological_inference_with_cubical_complexes::Distance_to_k_th_closest_point<periodic_Euclidean_distance_squared>
		( point_cloud, period_eu , i ) );
	}  
	
	
	
	//150 random points to test the distance:
	std::vector< std::vector<double> > test_points =
	{
    {0.8212883368,0.9323835692,0.1940174669,0.1572616117},{0.8278029691,0.1360385891,0.319820496,0.5313325708},{0.3791655428,0.6726754943,0.4926699875,0.4564590848},
	{0.4690400108,0.2233169088,0.8851099976,0.962797866},{0.9169740682,0.7425048782,0.273605769,0.1857918468},{0.0603475252,0.5917361956,0.4607867487,0.6087088764},
	{0.5774184729,0.5252969777,0.3485118898,0.7196434641},{0.2960761792,0.1822913771,0.095581162,0.5431221263},{0.8644298147,0.2528241375,0.0786190066,0.7023359034},
	{0.7118016153,0.7744334091,0.4906539726,0.5824714976},{0.3051843431,0.6273342867,0.6362808293,0.6466759269},{0.0039319792,0.2000689865,0.9937741337,0.3070294219},
	{0.0044878707,0.5990754166,0.8037203881,0.5964336854},{0.4894565104,0.2266998396,0.4111016593,0.2347550471},{0.1824641675,0.1050906335,0.9326545962,0.9113309076},
	{0.2736623602,0.3678764997,0.1222360081,0.6234545144},{0.4156279133,0.567760671,0.0257940979,0.5717319997},{0.7782721445,0.2349191937,0.509351237,0.3980594773},
	{0.2643903233,0.4370372372,0.652751063,0.3665730103},{0.8815587896,0.0884902957,0.0368820874,0.3579230048},{0.7030115363,0.8036939281,0.2909499747,0.0831794378},
	{0.7496797275,0.0440072753,0.6118866967,0.0836046506},{0.547318938,0.7143913424,0.9479668851,0.8093799036},{0.3125984911,0.4593533373,0.539115974,0.8845099781},
	{0.6280563318,0.7782411044,0.9531585327,0.9004038912},{0.1372506365,0.5994290735,0.6813214177,0.2109666273},{0.2932690966,0.356553782,0.5277665989,0.409365017},
	{0.2584309243,0.0374846419,0.8003505901,0.5980306529},{0.4126640663,0.842117093,0.0805273314,0.4085733381},{0.6963340475,0.0788598373,0.5735952742,0.0623137136},
	{0.4921945077,0.0705234474,0.3346862411,0.7895697765},{0.4298108651,0.4701639407,0.9527220742,0.2533458215},{0.7194097794,0.0143105998,0.1430793703,0.6942974012},
	{0.0974518694,0.1685294074,0.995312826,0.8557458415},{0.1206104693,0.3122800179,0.9435293309,0.8277134791},{0.4957749804,0.556660573,0.0333856768,0.9662811328},
	{0.5047073457,0.8086313256,0.1043228635,0.0894169754},{0.1364972428,0.8070002738,0.6322037117,0.378678191},{0.4142133044,0.9022612546,0.6024251701,0.2312628429},
	{0.568076449,0.4356079088,0.3142526348,0.168163758},{0.5078567564,0.5690632688,0.3501578458,0.9163810769},{0.844418278,0.6727674475,0.0808001885,0.2138062643},
	{0.9504417959,0.2526346936,0.5921653241,0.4655109756},{0.5316570103,0.0480489123,0.9929036549,0.5633992346},{0.4578632398,0.9795822352,0.2667927761,0.7587279216},
	{0.4388266846,0.7511075051,0.8202521054,0.7180580611},{0.7534749946,0.952236762,0.5843379819,0.7115072433},{0.082605399,0.4273217756,0.3213982,0.2979257607},
	{0.2389182923,0.9083379556,0.3061810983,0.1688461353},{0.7554514646,0.2075481357,0.7696146015,0.0587229319},{0.1392995564,0.1688323335,0.4324479382,0.7296311036},
	{0.4620787606,0.6012671189,0.2642026399,0.289015061},{0.6163531339,0.305217291,0.0370885001,0.8070358008},{0.4411920793,0.6116947453,0.7699594698,0.7750890346},
	{0.6547697294,0.2096259503,0.7220228808,0.5579349836},{0.5650564132,0.9872080924,0.3335652947,0.4187995577},{0.3477798183,0.5968151875,0.9324476945,0.5770690262},
	{0.6453652238,0.2499848837,0.3044766523,0.0086027875},{0.9733923627,0.7567299451,0.930025806,0.8176017015},{0.7080920769,0.5054853631,0.9846395506,0.0078840693},
	{0.3508421381,0.3803836836,0.7224787688,0.6663569745},{0.6061761128,0.2027737375,0.5866615586,0.1747405585},{0.14487503,0.5968572146,0.4438704676,0.8298940249},
	{0.7310764417,0.1569997829,0.0215511122,0.6790861671},{0.972724688,0.8266677347,0.6363947629,0.2776682025},{0.6350256647,0.5276081464,0.3509130494,0.5236328584},
	{0.7530832237,0.2824680347,0.297566067,0.4665729476},{0.5298101068,0.3094308714,0.9016442844,0.2776958565},{0.4552057874,0.9791901188,0.7813639792,0.8012753038},
	{0.2698913116,0.2485755188,0.026619998,0.4234410296},{0.9437902363,0.1250614747,0.3625405419,0.5881153196},{0.0210727728,0.1865349938,0.5801751288,0.7010335349},
	{0.379988872,0.2493978483,0.0100528188,0.9897437978},{0.3330588217,0.0681916452,0.5993071343,0.310666007},{0.0678583481,0.7531240345,0.6921427962,0.9757962844},
	{0.3633407755,0.7788240199,0.0492652224,0.1798770491},{0.3644540526,0.666081354,0.4692665741,0.8864048931},{0.3809871299,0.963774354,0.6546022166,0.0812799898},
	{0.7600035591,0.9536345575,0.9952779922,0.9471432264},{0.7484435302,0.7401210254,0.6538912386,0.5343515719},{0.5182223762,0.2501847574,0.5248613374,0.2511267415},
	{0.870385376,0.9125930276,0.7259652191,0.298885145},{0.3968986368,0.2084553535,0.252257138,0.8872261166},{0.6226723371,0.1969264294,0.3852901196,0.4376316178},
	{0.597791444,0.9029960178,0.711057466,0.0935030601},{0.6070077359,0.83265277,0.1798961195,0.0366028517},{0.4236347917,0.7429085399,0.6370544124,0.9928254755},
	{0.7458736077,0.8505542614,0.0402770936,0.0593888904},{0.7101504609,0.155220282,0.7446765681,0.136858959},{0.7322202548,0.2679031657,0.0821328175,0.509475555},
	{0.3194076915,0.5857978296,0.142673234,0.4252319071},{0.8345975818,0.1324339285,0.1419702419,0.9519022231},{0.4677787514,0.8552821742,0.0701958707,0.6066120334},
	{0.0781411156,0.4336681184,0.7544679798,0.5528892693},{0.4690446341,0.7194810561,0.8873615207,0.4407964002},{0.6952171784,0.2912570324,0.4385715786,0.7828665404},
	{0.6395390311,0.6555431942,0.2789093792,0.0409351103},{0.8470924867,0.0225327257,0.8786042873,0.0315075989},{0.1136108062,0.1157168602,0.1151479227,0.7284565477},
	{0.2400109852,0.3520848716,0.4474164532,0.5024195253},{0.9868204636,0.2810682349,0.7242138409,0.1066441021},{0.5960941922,0.8775974365,0.3824219184,0.3542780208},
	{0.8545002637,0.2566439696,0.8116236029,0.8005692856},{0.0528010298,0.3028914153,0.1698182218,0.9578424878},{0.30687282,0.0163802253,0.5300991598,0.3148234084},
	{0.844853862,0.9723400001,0.9709698244,0.3836798114},{0.2996011388,0.3647791226,0.3967912565,0.804769835},{0.4388404957,0.0744674248,0.2168146712,0.5686755246},
	{0.3596344753,0.5750163288,0.827843966,0.7845171792},{0.2247795686,0.1692882418,0.8996852383,0.7610411681},{0.0615380472,0.3320678554,0.361293945,0.5495825282},
	{0.3598536041,0.3687895206,0.1466682267,0.4858825665},{0.1583273448,0.5435633028,0.842119056,0.8990433763},{0.2906893739,0.9992338961,0.0365902244,0.2515490062},
	{0.708700492,0.6068974198,0.2160222987,0.9092870934},{0.039755008,0.4774281885,0.8254781992,0.0544538819},{0.1710351154,0.9192502983,0.6021166332,0.2370720017},
	{0.5406602423,0.265420089,0.2110569589,0.8268453497},{0.8766309114,0.8712419572,0.4400364414,0.707734108},{0.4772796489,0.9419202821,0.0377061213,0.3069938442},
	{0.6366287735,0.5056157932,0.6075186774,0.0095591163},{0.9636435816,0.2276716309,0.7912478447,0.4505429803},{0.44003758,0.0816964463,0.7389633986,0.4780616283},
	{0.1236324313,0.3426056569,0.6533704132,0.716917925},{0.0089956555,0.4026452634,0.51142625,0.9524126311},{0.4314487611,0.8592355724,0.1587632687,0.3773702965},
	{0.1030082908,0.757335756,0.7177376123,0.5277655818},{0.3001537488,0.9666405551,0.6284318171,0.9788859088},{0.7184872362,0.411033236,0.3672082569,0.5313134938},
	{0.1407736824,0.686925529,0.0366793531,0.6887775373},{0.2582462602,0.8335112447,0.5239669241,0.7044642316},{0.0834274916,0.4183128239,0.9547965054,0.1445864537},
	{0.8125402953,0.194330411,0.0731127004,0.2454144256},{0.8005590437,0.3136733051,0.9915199603,0.4998389932},{0.757094119,0.2729283979,0.3822575053,0.7096403379},
	{0.229019637,0.0695607793,0.2090823546,0.2863386087},{0.9698746274,0.1490741572,0.5285111498,0.9226275445},{0.7017279847,0.8202711521,0.5225879175,0.9342761824},
	{0.504911415,0.0338309484,0.6523185109,0.2559272265},{0.0569188015,0.8451312648,0.7274005236,0.1772637055},{0.7712833968,0.252321901,0.0951467424,0.5321299012},
	{0.4056106803,0.5530326331,0.5664660218,0.3138524885},{0.6589405159,0.2750555868,0.5445726956,0.9707992866},{0.3487046889,0.2726835969,0.5993727229,0.2305923773},
	{0.7604714078,0.7029063362,0.3828110378,0.5409116901},{0.3817545557,0.1703660444,0.8216343466,0.7968471635},{0.281624191,0.506147047,0.5236769335,0.5702272665},
	{0.1489132931,0.6962481788,0.6821822962,0.5757386542},{0.2034474257,0.807975766,0.1412100238,0.6455826147},{0.7724329121,0.5498912982,0.2408004289,0.0549012059}
	};	
	
	
	std::vector<double> result = 
	{
	0.0866097,0.125324,0.159467,0.175506,0.188266,0.202629,0.218229,0.238788,0.258367,0.272335,
	0.287535,0.291796,0.302103,0.314302,0.328816,0.338456,0.359513,0.104478,0.116617,0.140207,
	0.154483,0.174251,0.19794,0.209221,0.227695,0.241821,0.255328,0.260923,0.279312,0.282396,
	0.300593,0.312468,0.326268,0.339791,0.0685421,0.110851,0.132986,0.151099,0.174988,0.188392,
	0.200116,0.215879,0.227822,0.236895,0.252518,0.265997,0.287909,0.296851,0.305085,0.310382,
	0.327059,0.0835115,0.116422,0.155861,0.180695,0.200998,0.207592,0.221693,0.241037,0.256161,
	0.268959,0.301436,0.316509,0.321373,0.328739,0.355433,0.362623,0.374216,0.0513198,0.0792087,
	0.134415,0.165734,0.184031,0.200341,0.208617,0.225573,0.244771,0.258952,0.271968,0.277888,
	0.286844,0.293381,0.305178,0.324316,0.342208,0.0724997,0.108502,0.121614,0.145896,0.154156,
	0.174737,0.184184,0.198854,0.212056,0.240021,0.244794,0.257261,0.283923,0.299432,0.303161,
	0.315952,0.326779,0.0770476,0.107026,0.122477,0.135401,0.154982,0.167743,0.194236,0.212707,
	0.221344,0.237592,0.248993,0.263536,0.274802,0.278451,0.289309,0.306458,0.329846,0.061784,
	0.119136,0.148535,0.164464,0.171897,0.194357,0.219205,0.238269,0.248606,0.253468,0.263734,
	0.268021,0.276119,0.284281,0.303345,0.331981,0.34947,0.0693335,0.0966073,0.114494,0.131379,
	0.152204,0.170513,0.19538,0.217221,0.223733,0.241719,0.269773,0.288419,0.302165,0.308418,
	0.3335,0.35408,0.365201,0.0790481,0.116622,0.128577,0.152609,0.166748,0.191723,0.197587,
	0.206596,0.231318,0.252562,0.265415,0.274281,0.28668,0.292747,0.302854,0.313834,0.32478,
	0.0917386,0.122776,0.13392,0.185631,0.196152,0.20375,0.218107,0.228042,0.234576,0.246406,
	0.25741,0.263371,0.279396,0.292965,0.322036,0.340515,0.35218,0.0783718,0.131059,0.157269,
	0.192033,0.206777,0.224203,0.240961,0.246639,0.264559,0.271092,0.278815,0.306389,0.318983,
	0.332846,0.342972,0.359236,0.368122,0.0956505,0.115942,0.146143,0.170414,0.176753,0.204491,
	0.218034,0.229521,0.242055,0.260202,0.262535,0.282024,0.288318,0.291844,0.303969,0.313169,
	0.322649,0.068779,0.0801762,0.129558,0.156775,0.173762,0.199728,0.213599,0.220743,0.235786,
	0.254517,0.267364,0.283058,0.305141,0.317227,0.319803,0.367818,0.382489,0.0610831,0.0922659,
	0.124723,0.153995,0.165965,0.218019,0.238821,0.250126,0.279189,0.288199,0.292874,0.307344,
	0.319509,0.333863,0.341756,0.345347,0.360549,0.0605367,0.108147,0.127604,0.136212,0.164292,
	0.200737,0.206311,0.220285,0.224466,0.246824,0.259217,0.267648,0.277226,0.290026,0.303104,
	0.312308,0.325163,0.0641294,0.0968253,0.129973,0.156516,0.16666,0.195934,0.205461,0.235305,
	0.246166,0.256025,0.262578,0.269232,0.289435,0.303173,0.316643,0.337324,0.343012,0.0612598,
	0.124711,0.168873,0.191569,0.205565,0.215769,0.22248,0.234508,0.243369,0.264701,0.279279,
	0.298818,0.310758,0.322746,0.340821,0.357411,0.372742,0.0957317,0.127598,0.153114,0.177751,
	0.192934,0.211943,0.222981,0.231499,0.247166,0.255986,0.261922,0.275213,0.295218,0.314247,
	0.326658,0.338533,0.355751,0.121453,0.150036,0.165379,0.18332,0.195748,0.214231,0.228823,
	0.251365,0.262165,0.265568,0.282378,0.29385,0.302009,0.317558,0.323095,0.352809,0.363861,
	0.0715727,0.122657,0.144742,0.159225,0.166199,0.196014,0.202384,0.224946,0.240898,0.26239,
	0.269074,0.282638,0.293598,0.297921,0.315352,0.327804,0.342365,0.0726949,0.112302,0.120957,
	0.15468,0.191265,0.214065,0.224953,0.238154,0.247458,0.25768,0.283219,0.296907,0.31284,
	0.318608,0.324269,0.343991,0.361648,0.0892491,0.122093,0.148321,0.187194,0.197991,0.208214,
	0.220901,0.227742,0.240416,0.250237,0.274866,0.288959,0.298797,0.326771,0.334757,0.340789,
	0.360899,0.0608986,0.0899089,0.127018,0.149539,0.169289,0.184124,0.208084,0.227737,0.238195,
	0.253864,0.267916,0.278967,0.288704,0.302056,0.313724,0.325431,0.340671,0.0880712,0.137694,
	0.15614,0.194605,0.208665,0.222524,0.237535,0.250299,0.263831,0.269494,0.275921,0.289822,
	0.301681,0.313393,0.321955,0.34341,0.350783,0.0815019,0.10582,0.150856,0.192749,0.197319,
	0.221869,0.231888,0.24647,0.263957,0.275979,0.282388,0.297067,0.311492,0.321032,0.330129,
	0.336379,0.348431,0.0882601,0.133602,0.143852,0.15695,0.178324,0.189598,0.20582,0.211729,
	0.232506,0.252056,0.266006,0.276394,0.290782,0.298351,0.317209,0.332221,0.348685,0.0914995,
	0.124842,0.152299,0.173437,0.192333,0.209723,0.22475,0.231043,0.2421,0.259526,0.274164,
	0.288814,0.307066,0.317301,0.325682,0.331918,0.346532,0.0741936,0.10266,0.123711,0.150315,
	0.18033,0.200827,0.217878,0.233333,0.242617,0.257915,0.277601,0.290513,0.295838,0.307851,
	0.323902,0.335572,0.349632,0.0795895,0.09272,0.107584,0.15596,0.183878,0.213134,0.231215,
	0.244493,0.257193,0.268749,0.286588,0.287938,0.293456,0.305612,0.321109,0.332637,0.351873,
	0.0408375,0.0826852,0.140027,0.161198,0.174718,0.198565,0.214005,0.227787,0.246141,0.267658,
	0.286096,0.292201,0.296739,0.311791,0.317159,0.329718,0.351254,0.0942555,0.137532,0.153218,
	0.164805,0.179295,0.209377,0.220127,0.233618,0.244996,0.256981,0.269251,0.276568,0.292969,
	0.308634,0.320397,0.334288,0.354947,0.0961236,0.135735,0.147026,0.164276,0.171454,0.192377,
	0.212539,0.225045,0.244028,0.257454,0.266532,0.273357,0.286257,0.301753,0.306998,0.323628,
	0.333987,0.0497792,0.0759842,0.118954,0.167147,0.188571,0.201003,0.223686,0.239499,0.249619,
	0.265449,0.278572,0.296129,0.309439,0.317208,0.33409,0.348384,0.365242,0.0391462,0.0900126,
	0.123074,0.143729,0.151276,0.164739,0.198676,0.217583,0.235066,0.265741,0.278013,0.292401,
	0.309013,0.322635,0.336504,0.347827,0.365956,0.0748703,0.143573,0.153131,0.167772,0.196315,
	0.211537,0.229611,0.236906,0.256497,0.270613,0.275117,0.28644,0.290044,0.316725,0.329796,
	0.338824,0.35595,0.085887,0.138247,0.155707,0.184998,0.19304,0.205036,0.226871,0.233774,
	0.248536,0.271459,0.279219,0.292937,0.29832,0.31757,0.33432,0.346091,0.360739,0.0928386,
	0.12704,0.138961,0.15889,0.174679,0.206007,0.21826,0.227383,0.240309,0.251842,0.270802,
	0.282871,0.297702,0.306995,0.321059,0.330588,0.342693,0.0663763,0.0859251,0.11808,0.16228,
	0.177882,0.199928,0.220991,0.240196,0.252693,0.263487,0.273606,0.284867,0.296141,0.317397,
	0.32462,0.348309,0.354792,0.0821784,0.119541,0.137463,0.155327,0.188503,0.202968,0.218645,
	0.234524,0.240787,0.244127,0.273046,0.280023,0.285655,0.300656,0.309052,0.323299,0.338604,
	0.0991852,0.12733,0.131931,0.144924,0.155779,0.169712,0.187696,0.196109,0.206093,0.224401,
	0.24703,0.264852,0.272773,0.286893,0.299019,0.309315,0.319852,0.0830694,0.111396,0.127359,
	0.147548,0.171328,0.207739,0.225412,0.249429,0.263104,0.276446,0.287832,0.310363,0.320783,
	0.328947,0.344325,0.34986,0.383047,0.0874361,0.122702,0.148753,0.176472,0.203901,0.2169,
	0.231378,0.245832,0.257375,0.264989,0.274648,0.282829,0.300013,0.308324,0.329985,0.333259,
	0.343251,0.0884037,0.121899,0.153379,0.179372,0.201228,0.228302,0.245323,0.253796,0.265163,
	0.279425,0.290852,0.308346,0.31918,0.333683,0.345473,0.35975,0.369886,0.0438193,0.0846735,
	0.112045,0.161229,0.200156,0.223314,0.230548,0.248852,0.258213,0.263938,0.270567,0.281429,
	0.295858,0.311941,0.322963,0.331735,0.33723,0.0898389,0.123435,0.13549,0.167141,0.193621,
	0.223138,0.245095,0.252018,0.261089,0.266054,0.274821,0.280248,0.291952,0.312171,0.319247,
	0.337419,0.340997,0.0757395,0.103078,0.12992,0.149682,0.167497,0.188656,0.210744,0.220612,
	0.236916,0.25018,0.260414,0.270016,0.283781,0.312304,0.322618,0.335296,0.352874,0.0646964,
	0.102234,0.13339,0.171403,0.1854,0.198939,0.212749,0.221991,0.229961,0.248018,0.256368,
	0.265681,0.277281,0.295467,0.311712,0.321506,0.339707,0.0614896,0.124465,0.146145,0.163679,
	0.183217,0.194214,0.225156,0.229859,0.233289,0.247264,0.262113,0.279266,0.284534,0.294258,
	0.306311,0.337329,0.342652,0.100996,0.136054,0.168452,0.177143,0.207277,0.233349,0.240245,
	0.247784,0.256557,0.264866,0.280801,0.295236,0.316373,0.335641,0.345588,0.356674,0.364501,
	0.113615,0.126265,0.137569,0.149641,0.167878,0.197135,0.216156,0.229975,0.236185,0.250339,
	0.259254,0.266057,0.278863,0.294072,0.303504,0.325873,0.339713,0.0622003,0.10953,0.131944,
	0.158866,0.170058,0.191859,0.206908,0.212745,0.22427,0.240502,0.256132,0.273809,0.294221,
	0.299023,0.312715,0.318466,0.328725,0.080588,0.127807,0.146645,0.171701,0.176168,0.198319,
	0.21268,0.229767,0.243313,0.260777,0.272898,0.293153,0.302359,0.312851,0.32506,0.329043,
	0.343662,0.0917355,0.107191,0.127893,0.175457,0.188276,0.203196,0.214822,0.235494,0.242891,
	0.271935,0.297286,0.304847,0.312831,0.32011,0.332003,0.3389,0.350738,0.0794702,0.134934,
	0.153001,0.188422,0.199888,0.218363,0.229559,0.240265,0.250144,0.269981,0.280699,0.304677,
	0.316568,0.330702,0.349111,0.355207,0.362931,0.0712038,0.131217,0.14486,0.156068,0.174391,
	0.198625,0.213459,0.220099,0.237021,0.248379,0.268483,0.274627,0.285641,0.298111,0.312112,
	0.331654,0.35263,0.0859684,0.108655,0.132032,0.154179,0.163348,0.188725,0.219401,0.231379,
	0.249295,0.262776,0.268252,0.291761,0.296496,0.308749,0.326946,0.332517,0.34272,0.0604923,
	0.0965592,0.118084,0.142133,0.170715,0.197492,0.21756,0.231135,0.237452,0.25676,0.265766,
	0.285056,0.300225,0.305928,0.329366,0.346901,0.356555,0.100459,0.12185,0.138151,0.16448,
	0.178396,0.19812,0.210938,0.217487,0.232068,0.244832,0.254612,0.278257,0.288287,0.302087,
	0.317163,0.335666,0.346261,0.0744918,0.137403,0.178545,0.20112,0.213921,0.229158,0.235943,
	0.249066,0.261907,0.270161,0.288761,0.303677,0.310612,0.327844,0.338527,0.350215,0.357326,
	0.0726899,0.114705,0.148934,0.174076,0.189532,0.213648,0.220725,0.232535,0.248909,0.265304,
	0.27706,0.281436,0.290919,0.307221,0.314537,0.327124,0.336118,0.0392844,0.109267,0.142175,
	0.147805,0.172161,0.192495,0.208501,0.249398,0.262134,0.277359,0.301428,0.314384,0.324151,
	0.336354,0.343642,0.355436,0.37063,0.0767362,0.0967328,0.122746,0.14268,0.160642,0.179364,
	0.189245,0.211358,0.218854,0.242638,0.256464,0.262695,0.266503,0.283954,0.304517,0.314057,
	0.327732,0.0772157,0.114239,0.125842,0.161884,0.178399,0.202112,0.219612,0.244349,0.265477,
	0.28005,0.293964,0.301829,0.32143,0.330037,0.338636,0.343926,0.360834,0.0850882,0.111993,
	0.148329,0.172057,0.189524,0.213236,0.22661,0.251363,0.263,0.277996,0.285864,0.301574,
	0.307207,0.324508,0.331168,0.345006,0.349815,0.0768453,0.106453,0.131103,0.1448,0.174137,
	0.202092,0.216306,0.221784,0.238828,0.247637,0.252568,0.268753,0.288929,0.293761,0.303363,
	0.306932,0.317041,0.0749116,0.104201,0.140335,0.167106,0.191231,0.205081,0.220393,0.22497,
	0.237117,0.255584,0.274871,0.293,0.308767,0.313705,0.325572,0.33344,0.345504,0.112425,
	0.145085,0.174019,0.21364,0.218983,0.229121,0.23907,0.258361,0.267786,0.281671,0.288335,
	0.300495,0.310846,0.333876,0.347286,0.353972,0.360197,0.114881,0.141465,0.168672,0.179058,
	0.189807,0.201754,0.213917,0.233079,0.239043,0.246625,0.256755,0.28447,0.297504,0.310206,
	0.322354,0.330466,0.346788,0.0668362,0.0941798,0.156804,0.180992,0.194353,0.205086,0.217438,
	0.240646,0.259589,0.268828,0.276649,0.282394,0.294358,0.30961,0.322069,0.333912,0.345657,
	0.0788634,0.107873,0.125384,0.160394,0.176808,0.195176,0.207734,0.218926,0.240365,0.260119,
	0.264112,0.289572,0.298589,0.307487,0.320737,0.330111,0.339929,0.0703792,0.112669,0.129513,
	0.150632,0.172436,0.184785,0.205699,0.220392,0.242551,0.254562,0.270966,0.283158,0.286289,
	0.315834,0.329296,0.344201,0.351778,0.0690751,0.114319,0.152391,0.166012,0.18931,0.20815,
	0.216895,0.230242,0.245848,0.263272,0.277203,0.299106,0.311008,0.32847,0.346993,0.355618,
	0.367335,0.0630264,0.118122,0.140714,0.164881,0.185818,0.19113,0.208594,0.237461,0.250509,
	0.263443,0.279506,0.284324,0.292117,0.319068,0.327881,0.336775,0.347627,0.0760772,0.122869,
	0.15636,0.179355,0.201149,0.212025,0.225639,0.242882,0.252217,0.265287,0.281016,0.289248,
	0.292147,0.298521,0.311098,0.324454,0.33762,0.0927185,0.119665,0.139867,0.164137,0.183787,
	0.19952,0.22158,0.23188,0.256252,0.263369,0.280217,0.293577,0.310009,0.316349,0.332449,
	0.339481,0.349364,0.0797347,0.118081,0.139748,0.153132,0.168745,0.183387,0.19849,0.204668,
	0.22541,0.237048,0.241634,0.256276,0.265419,0.281252,0.293998,0.300715,0.315721,0.0731127,
	0.111941,0.13677,0.163443,0.188737,0.197857,0.209361,0.225942,0.251,0.266321,0.27914,
	0.299717,0.317422,0.320806,0.330218,0.348269,0.359229,0.102284,0.119777,0.142937,0.15857,
	0.180401,0.217232,0.235617,0.25172,0.272078,0.287679,0.29556,0.305355,0.319068,0.330144,
	0.345311,0.350405,0.365992,0.0832736,0.121844,0.161151,0.182954,0.195955,0.2071,0.216066,
	0.222344,0.243316,0.261911,0.272174,0.28418,0.299429,0.316116,0.324338,0.345983,0.350908,
	0.0708657,0.0891118,0.12455,0.152897,0.171665,0.193708,0.223837,0.242997,0.249632,0.267399,
	0.293357,0.301859,0.310979,0.322583,0.335521,0.354081,0.370786,0.0884469,0.122106,0.167167,
	0.183091,0.194101,0.216149,0.230086,0.250652,0.270133,0.278237,0.284311,0.289981,0.299481,
	0.323254,0.337247,0.351064,0.369168,0.0599319,0.0875535,0.131651,0.15685,0.166842,0.179865,
	0.201158,0.218273,0.233975,0.257614,0.27169,0.279539,0.294794,0.301501,0.325893,0.33156,
	0.336832,0.0942224,0.118426,0.148402,0.174459,0.18667,0.199476,0.215677,0.229703,0.253635,
	0.27137,0.27869,0.295702,0.30324,0.309264,0.32085,0.333331,0.33955,0.0756793,0.112726,
	0.141026,0.165351,0.186983,0.210219,0.229916,0.250261,0.263287,0.272302,0.287889,0.296794,
	0.302923,0.31654,0.330686,0.345761,0.356344,0.0836171,0.11411,0.14007,0.168675,0.2076,
	0.222609,0.239705,0.245646,0.252493,0.264371,0.28003,0.29083,0.298695,0.307116,0.316116,
	0.329705,0.341579,0.0578763,0.115005,0.1409,0.165111,0.181799,0.200973,0.213138,0.225617,
	0.23601,0.256892,0.273024,0.291231,0.296728,0.306326,0.319511,0.33514,0.346151,0.101478,
	0.143856,0.16839,0.180834,0.215203,0.22731,0.240905,0.250313,0.260393,0.2641,0.284022,
	0.301775,0.309147,0.319524,0.337031,0.351068,0.367746,0.0863714,0.144618,0.181002,0.191477,
	0.207476,0.219968,0.230643,0.258247,0.278749,0.287187,0.295456,0.303397,0.315763,0.326518,
	0.335684,0.345821,0.354184,0.072736,0.130302,0.151094,0.172721,0.189127,0.213904,0.220457,
	0.237709,0.251279,0.267584,0.278921,0.293285,0.309073,0.322461,0.334776,0.346181,0.366824,
	0.0610155,0.0834556,0.101004,0.129648,0.175243,0.200372,0.208401,0.223487,0.240048,0.241287,
	0.253385,0.266108,0.275223,0.297898,0.3094,0.318562,0.347239,0.088783,0.122642,0.150914,
	0.178254,0.195109,0.211548,0.232988,0.237384,0.252011,0.26643,0.276273,0.2956,0.308303,
	0.322157,0.341267,0.353028,0.360762,0.0571079,0.112126,0.139146,0.161584,0.181456,0.191683,
	0.202602,0.235365,0.259493,0.270633,0.277032,0.293519,0.319303,0.331497,0.336299,0.343372,
	0.354646,0.0869447,0.11709,0.136713,0.154159,0.164471,0.189594,0.219064,0.231908,0.237789,
	0.256803,0.267706,0.274035,0.296756,0.31176,0.319487,0.330123,0.338544,0.0604749,0.127838,
	0.153142,0.165924,0.196884,0.21114,0.228728,0.243872,0.255087,0.266564,0.286064,0.289592,
	0.298367,0.314313,0.329439,0.351998,0.37327,0.0674806,0.0979469,0.127132,0.159232,0.167303,
	0.194463,0.210262,0.216675,0.238042,0.24526,0.257838,0.279103,0.289669,0.294223,0.298902,
	0.306579,0.316529,0.0774721,0.101642,0.124886,0.156422,0.187589,0.195694,0.210388,0.221973,
	0.227985,0.247221,0.257485,0.273696,0.282306,0.301981,0.314224,0.322761,0.331448,0.0925958,
	0.136807,0.169163,0.185664,0.195215,0.2112,0.225809,0.24335,0.250636,0.260341,0.286231,
	0.303753,0.31758,0.334553,0.35033,0.35641,0.373486,0.0586043,0.108019,0.124277,0.137655,
	0.16665,0.192384,0.207246,0.214843,0.23039,0.257324,0.264713,0.277807,0.280874,0.288152,
	0.297467,0.327447,0.333905,0.0971172,0.124306,0.14288,0.155493,0.16972,0.180381,0.198003,
	0.220574,0.233381,0.251591,0.268809,0.281429,0.287686,0.309842,0.31516,0.322481,0.33073,
	0.105097,0.126728,0.144192,0.175846,0.181042,0.202362,0.217117,0.237092,0.255407,0.263605,
	0.284819,0.296108,0.30284,0.315827,0.330997,0.344201,0.357134,0.0691672,0.104287,0.135275,
	0.168289,0.179279,0.188852,0.208246,0.215645,0.224918,0.236623,0.253709,0.269839,0.279728,
	0.299555,0.307055,0.33142,0.341521,0.0340118,0.13623,0.153594,0.165725,0.174385,0.190972,
	0.203867,0.223174,0.235996,0.253889,0.268241,0.285812,0.302672,0.316463,0.327054,0.339961,
	0.351844,0.0943352,0.115839,0.13336,0.162647,0.182941,0.191529,0.208967,0.219414,0.231203,
	0.245279,0.256459,0.272936,0.286895,0.292727,0.304123,0.328313,0.340691,0.0550216,0.110578,
	0.145821,0.166149,0.180678,0.2028,0.220144,0.240504,0.248595,0.261185,0.269743,0.283433,
	0.299313,0.302738,0.308811,0.325226,0.337577,0.120465,0.165244,0.168282,0.178406,0.189637,
	0.207124,0.220658,0.23599,0.249548,0.270185,0.276632,0.296876,0.313442,0.32621,0.335114,
	0.349813,0.366306,0.0636063,0.0997389,0.122595,0.151887,0.168329,0.196454,0.218305,0.227159,
	0.233997,0.24372,0.261896,0.276295,0.287373,0.292355,0.301144,0.310144,0.318296,0.0761112,
	0.115902,0.133989,0.150793,0.174443,0.19073,0.214835,0.234738,0.257601,0.266254,0.281151,
	0.296429,0.303394,0.326666,0.338239,0.344329,0.353129,0.0829953,0.113253,0.141536,0.161254,
	0.175881,0.194755,0.219118,0.227507,0.253025,0.261269,0.280762,0.293031,0.306985,0.309005,
	0.334526,0.342887,0.349116,0.0669628,0.0911037,0.108613,0.151058,0.17831,0.203073,0.234575,
	0.245622,0.252097,0.270338,0.285851,0.297895,0.309903,0.328699,0.346095,0.361071,0.375646,
	0.0742994,0.0911847,0.130471,0.172131,0.185628,0.208557,0.215593,0.223525,0.232455,0.243552,
	0.251835,0.260496,0.269395,0.276546,0.283039,0.299985,0.319843,0.0806465,0.104794,0.131362,
	0.14777,0.176275,0.2053,0.215809,0.226753,0.235692,0.244471,0.265634,0.273161,0.285807,
	0.293786,0.31317,0.329131,0.340454,0.0720657,0.125007,0.143192,0.157986,0.177412,0.193465,
	0.207258,0.226482,0.240785,0.265872,0.277719,0.281239,0.294046,0.300999,0.312408,0.319344,
	0.33827,0.0855427,0.117851,0.1597,0.174151,0.186415,0.208825,0.21485,0.240745,0.251324,
	0.263184,0.273975,0.290397,0.304306,0.31719,0.325133,0.348546,0.366005,0.0865814,0.107587,
	0.128564,0.144936,0.170277,0.180803,0.203571,0.219222,0.229753,0.23537,0.254037,0.268164,
	0.282044,0.306533,0.317725,0.329872,0.344548,0.0824779,0.128152,0.151422,0.171779,0.186215,
	0.213889,0.223892,0.24451,0.266028,0.27521,0.289762,0.294296,0.303598,0.3182,0.323909,
	0.343504,0.356392,0.0872621,0.115798,0.158029,0.161887,0.195026,0.213633,0.225745,0.232213,
	0.247867,0.260275,0.267328,0.279441,0.293108,0.308397,0.322491,0.333769,0.343857,0.048632,
	0.110428,0.142333,0.153489,0.166466,0.185041,0.197186,0.228679,0.253917,0.261973,0.273064,
	0.283013,0.29222,0.301292,0.308364,0.319135,0.324149,0.0713495,0.119152,0.132524,0.148131,
	0.161317,0.195961,0.206772,0.214898,0.229108,0.236556,0.257722,0.269625,0.27832,0.295376,
	0.312609,0.317899,0.332201,0.0885819,0.128667,0.160287,0.184927,0.197021,0.21447,0.223629,
	0.232221,0.249758,0.269039,0.281548,0.287262,0.295114,0.304859,0.324047,0.337538,0.366614,
	0.0738487,0.107811,0.131159,0.153789,0.164193,0.179242,0.19139,0.216516,0.226188,0.246075,
	0.263701,0.279958,0.29435,0.310588,0.331371,0.33784,0.362286,0.0752822,0.108601,0.14123,
	0.158362,0.182345,0.209964,0.214395,0.245556,0.267633,0.27568,0.294184,0.305045,0.322739,
	0.348439,0.358415,0.371849,0.386037,0.085353,0.131927,0.169685,0.17694,0.18574,0.195036,
	0.22927,0.240342,0.253878,0.266071,0.282581,0.301264,0.315549,0.325435,0.342666,0.354399,
	0.378471,0.0719061,0.0979673,0.127286,0.154186,0.159123,0.183349,0.19705,0.231169,0.237668,
	0.249969,0.264387,0.270985,0.294968,0.314376,0.317838,0.334414,0.345607,0.0858032,0.0978837,
	0.117975,0.15683,0.17609,0.196763,0.224105,0.233408,0.24756,0.257568,0.265986,0.282148,
	0.295271,0.303874,0.320034,0.325939,0.339522,0.0804856,0.106019,0.118565,0.136561,0.167915,
	0.180156,0.206613,0.226454,0.23552,0.261712,0.275658,0.287299,0.302089,0.312008,0.323509,
	0.331444,0.338995,0.0979951,0.127548,0.159247,0.182576,0.190378,0.206484,0.217948,0.224542,
	0.235132,0.247965,0.261603,0.277685,0.282364,0.294322,0.312957,0.331067,0.345518,0.0909583,
	0.116089,0.127927,0.173863,0.189903,0.200383,0.204621,0.224951,0.242193,0.257417,0.269236,
	0.276154,0.286126,0.310862,0.32457,0.334272,0.34765,0.0756548,0.101906,0.125566,0.163069,
	0.177698,0.186845,0.212928,0.232779,0.24488,0.25757,0.262906,0.285118,0.291149,0.299354,
	0.308125,0.328257,0.337387,0.0695973,0.0945282,0.12413,0.13665,0.163441,0.178353,0.193922,
	0.2146,0.22641,0.249348,0.256519,0.266675,0.282196,0.309411,0.318305,0.326,0.333543,
	0.0802617,0.116198,0.132775,0.172006,0.189522,0.199803,0.218355,0.226147,0.246105,0.25493,
	0.262287,0.268627,0.27468,0.288259,0.301248,0.314731,0.333959,0.0976088,0.116413,0.1371,
	0.162443,0.172151,0.200739,0.223263,0.242068,0.257668,0.270821,0.287911,0.295735,0.312491,
	0.322849,0.329495,0.341309,0.349157,0.140605,0.157256,0.181432,0.205278,0.221546,0.233098,
	0.237935,0.255805,0.265905,0.271016,0.27509,0.290212,0.298278,0.308765,0.335595,0.341964,
	0.357522,0.0650989,0.109293,0.126012,0.150571,0.201708,0.221332,0.242948,0.256907,0.264695,
	0.2785,0.289707,0.295681,0.311086,0.315735,0.341543,0.356904,0.381846,0.0666296,0.102536,
	0.130034,0.152498,0.172291,0.177294,0.203633,0.227108,0.240743,0.250142,0.257903,0.280665,
	0.29245,0.297218,0.304761,0.316983,0.330537,0.0759739,0.123526,0.152275,0.162357,0.180373,
	0.196809,0.203063,0.22699,0.238048,0.255064,0.263501,0.290106,0.30256,0.311425,0.331981,
	0.342248,0.355956,0.0929587,0.121635,0.135097,0.166535,0.192567,0.200047,0.215297,0.220377,
	0.230892,0.242567,0.265364,0.274777,0.287731,0.293379,0.313749,0.328445,0.343168,0.0695927,
	0.0885967,0.111911,0.143098,0.163242,0.181265,0.196372,0.215883,0.237569,0.25676,0.276172,
	0.284884,0.297006,0.303887,0.322557,0.336902,0.350339,0.0731276,0.108281,0.122628,0.141478,
	0.177175,0.205087,0.217407,0.240658,0.254903,0.266899,0.284993,0.301768,0.321723,0.337408,
	0.350338,0.364817,0.384992,0.0935971,0.105115,0.14644,0.170578,0.18814,0.219826,0.2352,
	0.262944,0.269516,0.280621,0.285375,0.294756,0.302603,0.315299,0.324174,0.34691,0.359126,
	0.0759081,0.116708,0.137456,0.166048,0.181472,0.199758,0.215416,0.222324,0.231229,0.254284,
	0.267913,0.293953,0.311728,0.324537,0.330346,0.347818,0.366289,0.0627016,0.108165,0.131921,
	0.165824,0.179121,0.196457,0.213357,0.228059,0.238969,0.259755,0.274208,0.281844,0.293119,
	0.307027,0.321349,0.337843,0.355811,0.0631621,0.105604,0.134326,0.153431,0.174078,0.191667,
	0.215092,0.223349,0.233288,0.240478,0.251417,0.262748,0.274282,0.284032,0.299785,0.323678,
	0.331263,0.0662261,0.094083,0.145618,0.166959,0.17664,0.186313,0.204052,0.229969,0.248943,
	0.265448,0.283015,0.304228,0.314081,0.322447,0.328711,0.343626,0.375055,0.0685608,0.0950153,
	0.133996,0.147783,0.171011,0.183074,0.199839,0.221071,0.231131,0.242096,0.251254,0.266953,
	0.290816,0.293721,0.307981,0.319196,0.324223,0.109615,0.121177,0.146341,0.160705,0.169396,
	0.188269,0.217159,0.230534,0.249595,0.261825,0.280702,0.301303,0.311353,0.326552,0.334763,
	0.350327,0.36518,0.0813862,0.101185,0.123248,0.164334,0.176006,0.208952,0.220644,0.230336,
	0.246282,0.261614,0.266886,0.272053,0.290027,0.301495,0.313338,0.325855,0.332749,0.10676,
	0.136148,0.151959,0.171667,0.180044,0.192493,0.202568,0.213818,0.2323,0.239598,0.266103,
	0.280737,0.290539,0.300109,0.314691,0.326531,0.340611,0.0641875,0.0982165,0.127956,0.141531,
	0.160579,0.174176,0.197402,0.203905,0.217523,0.231737,0.246484,0.264925,0.287368,0.299377,
	0.330067,0.338877,0.342968,0.0771866,0.100616,0.123925,0.150439,0.180829,0.214389,0.228747,
	0.236533,0.255469,0.266065,0.273667,0.282974,0.29297,0.315278,0.323968,0.333982,0.347466
	};
	
	
	//for every test point:
	unsigned counter = 0;
	for ( size_t pts = 0 ; pts != test_points.size() ; ++pts )
	{
		//for every distance class:
		for ( size_t dist = 0 ; dist != classes_to_compute_knn.size() ; ++dist )
		{
			//std::cout << (*classes_to_compute_knn[dist])( test_points[pts]) << ",";
			//if ( counter % 10 == 9 )  std::cout << std::endl;
			BOOST_CHECK( fabs( (*classes_to_compute_knn[dist])( test_points[pts]) - result[counter] ) <= 5e-07 );
			++counter;
		}
	}		
}	





BOOST_AUTO_TEST_CASE(Distance_to_k_th_nearest_neighbor_periodic_domain_4d_kd_trees)
{
	//now the k_d_trees:
	
	std::vector< std::pair< double , double > > coordinates_of_grid(4);
	coordinates_of_grid[0] = std::pair< double,double >(0,1);
	coordinates_of_grid[1] = std::pair< double,double >(0,1);
	coordinates_of_grid[2] = std::pair< double,double >(0,1);
	coordinates_of_grid[3] = std::pair< double,double >(0,1);
	
	
	
	
	//150 random points from R^4.
	std::vector< std::vector<double> > point_cloud = 
	{
	{0.2161353165,0.4067554909,0.0665765158,0.870283464},{0.7782972094,0.6987265884,0.7606270006,0.7958696354},{0.9662634397,0.1367719397,0.8092096795,0.716037544},
	{0.275288115,0.1794606529,0.8734926581,0.1608733959},{0.3935741177,0.5092781719,0.5718593814,0.355841625},{0.9207663401,0.5319217448,0.258124053,0.1236327693},
	{0.8471451914,0.4675919802,0.4061658739,0.9562745956},{0.9294778437,0.1459124812,0.5308060655,0.4604051244},{0.3011312045,0.4122534681,0.337005893,0.9334948768},
	{0.0638816946,0.1186342924,0.7877057337,0.6814819723},{0.3170833248,0.5444929474,0.15133534,0.7217780657},{0.4718955676,0.0654609054,0.126377013,0.9821652356},
	{0.7545759601,0.0546513223,0.8329449217,0.6777689687},{0.5130777196,0.500694545,0.6574573508,0.7392099628},{0.5893094379,0.0494349443,0.5050142871,0.0808209586},
	{0.4940861093,0.2490698399,0.3612209754,0.7160524682},{0.7133854455,0.0851763973,0.4031818854,0.2452370371},{0.6502694073,0.6256448289,0.0068036765,0.9420920806},
	{0.2771677715,0.5365420913,0.2049486258,0.6807608032},{0.6708589897,0.2624673771,0.5366855075,0.036167918},{0.6042193428,0.3100096218,0.0774625444,0.9480968153},
	{0.8630631317,0.4616721999,0.1888999115,0.4816821855},{0.274362145,0.6687946867,0.9600374054,0.5621335753},{0.9371821445,0.8886038563,0.4903153491,0.5409146673},
	{0.1577256215,0.9017276918,0.6047877711,0.2715029821},{0.4778309672,0.3667936909,0.2353694309,0.2162692775},{0.9561227295,0.2203476694,0.8832569388,0.8253242387},
	{0.3162962156,0.2985790542,0.3671411616,0.1277638588},{0.9136061701,0.1122253346,0.5190847388,0.3615447618},{0.327857591,0.6936850683,0.2818536318,0.5899669132},
	{0.6639132625,0.6278195342,0.8459756977,0.9733827612},{0.8529426146,0.8311685929,0.5352412218,0.5575795863},{0.8291846917,0.7364060413,0.5959779969,0.8778938674},
	{0.7160940336,0.7085216588,0.2321347578,0.2757493216},{0.5377012696,0.9137963476,0.3599800142,0.7714768103},{0.5407209711,0.9124356993,0.0977215923,0.6965975054},
	{0.183090288,0.4487572513,0.3782904518,0.515378186},{0.2529464646,0.2026785347,0.7370024058,0.8433211406},{0.7773119675,0.5135940188,0.0091130023,0.5177475526},
	{0.2755686354,0.5956264965,0.5519859341,0.6927454374},{0.5180953105,0.7024236836,0.873931882,0.3592949654},{0.9604220481,0.627821804,0.0671470871,0.8844678551},
	{0.2330840409,0.819831975,0.1885772266,0.0399145884},{0.6747471029,0.4076069591,0.3782148864,0.2600041842},{0.938375168,0.8666493378,0.9016399393,0.4959568656},
	{0.4837143579,0.5254478094,0.5813259515,0.5528810667},{0.4965164606,0.2691201437,0.9339421515,0.607314809},{0.3963730407,0.0023954718,0.5331149858,0.4707511207},
	{0.0545727671,0.3826355562,0.3722823337,0.7533125849},{0.7794329945,0.6273420472,0.475966346,0.7291254057},{0.302274507,0.0816056114,0.415606502,0.727463753},
	{0.8516243717,0.5476219286,0.1300290304,0.8781453874},{0.7446079317,0.4378378205,0.212476521,0.6074005007},{0.0543138136,0.0176691897,0.7037446832,0.3270085847},
	{0.9654455474,0.6369031293,0.4479166996,0.1196174838},{0.1637513947,0.8642580179,0.2845635938,0.7220461073},{0.2128266117,0.9449356911,0.4176520605,0.591265711},
	{0.0476004935,0.286355197,0.2003003603,0.3544409689},{0.4617424291,0.34739498,0.5184426282,0.7242208878},{0.9389569708,0.1110231471,0.1960562568,0.8751947973},
	{0.1297257203,0.1082481672,0.9613296099,0.6461059342},{0.6314383652,0.863803653,0.4156785782,0.1647425564},{0.8742997588,0.7209185674,0.1300667289,0.9857225746},
	{0.4410702195,0.5473736799,0.3925960511,0.073340355},{0.2864084891,0.0529552582,0.4401347351,0.1997130518},{0.2518973188,0.1512412846,0.0166082208,0.244082951},
	{0.2608718651,0.3542067721,0.1776841788,0.5686277715},{0.0090282573,0.6212397611,0.6485787611,0.9783152305},{0.3855588885,0.0695148581,0.036164256,0.549218531},
	{0.1092119429,0.1176601963,0.9561683217,0.0636541741},{0.0621516255,0.4523660035,0.7733616829,0.5127586231},{0.6562764596,0.2179020974,0.5065915762,0.3336082227},
	{0.9608324207,0.2805596716,0.9524878454,0.4502769466},{0.4204976272,0.9747959129,0.6945173966,0.0523253432},{0.9322453286,0.9045769402,0.4068475547,0.7923901044},
	{0.0912914628,0.4110324632,0.1179300488,0.4090027225},{0.9372338173,0.0228120044,0.1099257688,0.7709252352},{0.2473172548,0.9430129114,0.9248211253,0.2755537038},
	{0.3151711877,0.3909966208,0.5856341401,0.8289346434},{0.6222296427,0.7623051493,0.7437775929,0.5864307864},{0.2865205696,0.4823481457,0.2886831271,0.16002433},
	{0.2867317577,0.3187225992,0.6152163229,0.9980220243},{0.814303647,0.6177073771,0.0959672085,0.5412626567},{0.5198940996,0.8805264027,0.5142570818,0.9733662861},
	{0.0274844745,0.193287927,0.969342147,0.9436305743},{0.9917424994,0.652558957,0.6949542321,0.9143778963},{0.158117207,0.4643659701,0.8407568894,0.79356269},
	{0.6943920886,0.5211990806,0.2490926124,0.58702025},{0.6541512765,0.5538034292,0.2002501043,0.499721779},{0.2718869599,0.667435063,0.1980533334,0.411654298},
	{0.2582200989,0.838053823,0.5876435062,0.2448927639},{0.452425757,0.8820036445,0.9339749734,0.4761730228},{0.6869763259,0.8356626607,0.4592519393,0.1507639568},
	{0.1855826536,0.1115194834,0.9367056177,0.6076889979},{0.2906259801,0.5802919739,0.2736881424,0.3915352314},{0.12288481,0.4239620366,0.0535635338,0.5429183329},
	{0.6517687093,0.0795862924,0.1951899431,0.8655872962},{0.4832760945,0.1621366784,0.7051915845,0.2026124394},{0.9915547213,0.7069917477,0.5824110755,0.843934183},
	{0.0752148305,0.6797068669,0.1407429664,0.2335030579},{0.602499583,0.2250958043,0.7167737996,0.5247004263},{0.0792115591,0.0035277302,0.9021162488,0.9836529621},
	{0.0630885197,0.6648809183,0.2074113868,0.986199873},{0.7841980574,0.054547969,0.3682132296,0.3902004855},{0.5905348728,0.179866985,0.3925852075,0.8659482158},
	{0.7789685093,0.7006078116,0.2386029409,0.5079281065},{0.8778327638,0.1721387503,0.4325646544,0.8477176272},{0.9655417092,0.8003545383,0.2018559808,0.447240192},
	{0.3959916683,0.716806635,0.3991919297,0.4022070405},{0.3342622917,0.6906662339,0.6621937822,0.1680860706},{0.6058353547,0.8673584729,0.7881982641,0.4259585107},
	{0.5841680258,0.7195332507,0.0491419374,0.1035793731},{0.4424925682,0.1699164014,0.0806449894,0.8657210511},{0.4598111401,0.7129467162,0.1730200783,0.6224096753},
	{0.9347412391,0.8530473537,0.3221795461,0.3831341839},{0.515093598,0.9728921819,0.3625707706,0.6163540559},{0.4335091764,0.9474623057,0.3618741892,0.6030456028},
	{0.6046422494,0.2836575734,0.6610661431,0.2691357303},{0.2177470808,0.5447009606,0.8826168361,0.1196069277},{0.3651480607,0.7555859217,0.0177540325,0.5050244427},
	{0.7203258723,0.3371161658,0.4509870093,0.9302381405},{0.0989933989,0.2220814053,0.050527609,0.5632122078},{0.369650627,0.5319386604,0.9319623525,0.9647439241},
	{0.8442949792,0.6522289219,0.7757588418,0.3020221314},{0.9597743074,0.558198594,0.9786667374,0.5753302046},{0.2922401265,0.3846208933,0.748712206,0.0844571677},
	{0.3025243522,0.1353433062,0.3826868525,0.1765479192},{0.7594928641,0.8685317831,0.2470071204,0.111033147},{0.9367945509,0.9874413582,0.6886330319,0.0756251651},
	{0.7691251929,0.0561589014,0.5874789746,0.006362404},{0.7453925973,0.2688428422,0.930551802,0.7137127921},{0.4403857084,0.0046569624,0.3807356679,0.8867194413},
	{0.5359529725,0.5732941926,0.6598083961,0.8664365266},{0.4278839084,0.705274327,0.7008550165,0.1929006591},{0.0876745624,0.7284803414,0.1906603929,0.3610600757},
	{0.9920528498,0.2349508412,0.8696436284,0.8481443948},{0.0071836323,0.377837914,0.4632825549,0.4406018984},{0.2740049621,0.4837526211,0.6996950118,0.2481896421},
	{0.5023594771,0.6058867699,0.5802190339,0.0271982334},{0.6278990563,0.8451684192,0.6360527,0.8441247668},{0.9253878356,0.4954766021,0.4400199461,0.7489505075},
	{0.5729722457,0.3896644767,0.3697294721,0.0447440601},{0.9223849133,0.5377308403,0.3215464076,0.7362179228},{0.9915811059,0.863216745,0.5659185979,0.913221871},
	{0.1882577084,0.9678615758,0.3086162433,0.2597987708},{0.8596413189,0.2405901959,0.8497658458,0.6602620669},{0.259099799,0.8367184943,0.1531046529,0.1532220987},
	{0.1685487258,0.7992893311,0.2069046353,0.5789962241},{0.8992731075,0.4039457724,0.5378789764,0.682111514},{0.1674966314,0.6672669824,0.7791987706,0.5535525433}
	};	
	
	std::vector< Gudhi::Topological_inference_with_cubical_complexes::Distance_to_k_th_closest_point_periodic_k_d_tree* >
	classes_to_compute_knn;	
	for ( size_t i = 5 ; i != 90 ; i=i+5 )
	{
		classes_to_compute_knn.push_back(
		new Gudhi::Topological_inference_with_cubical_complexes::Distance_to_k_th_closest_point_periodic_k_d_tree
		( point_cloud, coordinates_of_grid , i ) );
	}  
	
	
	
	//150 random points to test the distance:
	std::vector< std::vector<double> > test_points =
	{
    {0.8212883368,0.9323835692,0.1940174669,0.1572616117},{0.8278029691,0.1360385891,0.319820496,0.5313325708},{0.3791655428,0.6726754943,0.4926699875,0.4564590848},
	{0.4690400108,0.2233169088,0.8851099976,0.962797866},{0.9169740682,0.7425048782,0.273605769,0.1857918468},{0.0603475252,0.5917361956,0.4607867487,0.6087088764},
	{0.5774184729,0.5252969777,0.3485118898,0.7196434641},{0.2960761792,0.1822913771,0.095581162,0.5431221263},{0.8644298147,0.2528241375,0.0786190066,0.7023359034},
	{0.7118016153,0.7744334091,0.4906539726,0.5824714976},{0.3051843431,0.6273342867,0.6362808293,0.6466759269},{0.0039319792,0.2000689865,0.9937741337,0.3070294219},
	{0.0044878707,0.5990754166,0.8037203881,0.5964336854},{0.4894565104,0.2266998396,0.4111016593,0.2347550471},{0.1824641675,0.1050906335,0.9326545962,0.9113309076},
	{0.2736623602,0.3678764997,0.1222360081,0.6234545144},{0.4156279133,0.567760671,0.0257940979,0.5717319997},{0.7782721445,0.2349191937,0.509351237,0.3980594773},
	{0.2643903233,0.4370372372,0.652751063,0.3665730103},{0.8815587896,0.0884902957,0.0368820874,0.3579230048},{0.7030115363,0.8036939281,0.2909499747,0.0831794378},
	{0.7496797275,0.0440072753,0.6118866967,0.0836046506},{0.547318938,0.7143913424,0.9479668851,0.8093799036},{0.3125984911,0.4593533373,0.539115974,0.8845099781},
	{0.6280563318,0.7782411044,0.9531585327,0.9004038912},{0.1372506365,0.5994290735,0.6813214177,0.2109666273},{0.2932690966,0.356553782,0.5277665989,0.409365017},
	{0.2584309243,0.0374846419,0.8003505901,0.5980306529},{0.4126640663,0.842117093,0.0805273314,0.4085733381},{0.6963340475,0.0788598373,0.5735952742,0.0623137136},
	{0.4921945077,0.0705234474,0.3346862411,0.7895697765},{0.4298108651,0.4701639407,0.9527220742,0.2533458215},{0.7194097794,0.0143105998,0.1430793703,0.6942974012},
	{0.0974518694,0.1685294074,0.995312826,0.8557458415},{0.1206104693,0.3122800179,0.9435293309,0.8277134791},{0.4957749804,0.556660573,0.0333856768,0.9662811328},
	{0.5047073457,0.8086313256,0.1043228635,0.0894169754},{0.1364972428,0.8070002738,0.6322037117,0.378678191},{0.4142133044,0.9022612546,0.6024251701,0.2312628429},
	{0.568076449,0.4356079088,0.3142526348,0.168163758},{0.5078567564,0.5690632688,0.3501578458,0.9163810769},{0.844418278,0.6727674475,0.0808001885,0.2138062643},
	{0.9504417959,0.2526346936,0.5921653241,0.4655109756},{0.5316570103,0.0480489123,0.9929036549,0.5633992346},{0.4578632398,0.9795822352,0.2667927761,0.7587279216},
	{0.4388266846,0.7511075051,0.8202521054,0.7180580611},{0.7534749946,0.952236762,0.5843379819,0.7115072433},{0.082605399,0.4273217756,0.3213982,0.2979257607},
	{0.2389182923,0.9083379556,0.3061810983,0.1688461353},{0.7554514646,0.2075481357,0.7696146015,0.0587229319},{0.1392995564,0.1688323335,0.4324479382,0.7296311036},
	{0.4620787606,0.6012671189,0.2642026399,0.289015061},{0.6163531339,0.305217291,0.0370885001,0.8070358008},{0.4411920793,0.6116947453,0.7699594698,0.7750890346},
	{0.6547697294,0.2096259503,0.7220228808,0.5579349836},{0.5650564132,0.9872080924,0.3335652947,0.4187995577},{0.3477798183,0.5968151875,0.9324476945,0.5770690262},
	{0.6453652238,0.2499848837,0.3044766523,0.0086027875},{0.9733923627,0.7567299451,0.930025806,0.8176017015},{0.7080920769,0.5054853631,0.9846395506,0.0078840693},
	{0.3508421381,0.3803836836,0.7224787688,0.6663569745},{0.6061761128,0.2027737375,0.5866615586,0.1747405585},{0.14487503,0.5968572146,0.4438704676,0.8298940249},
	{0.7310764417,0.1569997829,0.0215511122,0.6790861671},{0.972724688,0.8266677347,0.6363947629,0.2776682025},{0.6350256647,0.5276081464,0.3509130494,0.5236328584},
	{0.7530832237,0.2824680347,0.297566067,0.4665729476},{0.5298101068,0.3094308714,0.9016442844,0.2776958565},{0.4552057874,0.9791901188,0.7813639792,0.8012753038},
	{0.2698913116,0.2485755188,0.026619998,0.4234410296},{0.9437902363,0.1250614747,0.3625405419,0.5881153196},{0.0210727728,0.1865349938,0.5801751288,0.7010335349},
	{0.379988872,0.2493978483,0.0100528188,0.9897437978},{0.3330588217,0.0681916452,0.5993071343,0.310666007},{0.0678583481,0.7531240345,0.6921427962,0.9757962844},
	{0.3633407755,0.7788240199,0.0492652224,0.1798770491},{0.3644540526,0.666081354,0.4692665741,0.8864048931},{0.3809871299,0.963774354,0.6546022166,0.0812799898},
	{0.7600035591,0.9536345575,0.9952779922,0.9471432264},{0.7484435302,0.7401210254,0.6538912386,0.5343515719},{0.5182223762,0.2501847574,0.5248613374,0.2511267415},
	{0.870385376,0.9125930276,0.7259652191,0.298885145},{0.3968986368,0.2084553535,0.252257138,0.8872261166},{0.6226723371,0.1969264294,0.3852901196,0.4376316178},
	{0.597791444,0.9029960178,0.711057466,0.0935030601},{0.6070077359,0.83265277,0.1798961195,0.0366028517},{0.4236347917,0.7429085399,0.6370544124,0.9928254755},
	{0.7458736077,0.8505542614,0.0402770936,0.0593888904},{0.7101504609,0.155220282,0.7446765681,0.136858959},{0.7322202548,0.2679031657,0.0821328175,0.509475555},
	{0.3194076915,0.5857978296,0.142673234,0.4252319071},{0.8345975818,0.1324339285,0.1419702419,0.9519022231},{0.4677787514,0.8552821742,0.0701958707,0.6066120334},
	{0.0781411156,0.4336681184,0.7544679798,0.5528892693},{0.4690446341,0.7194810561,0.8873615207,0.4407964002},{0.6952171784,0.2912570324,0.4385715786,0.7828665404},
	{0.6395390311,0.6555431942,0.2789093792,0.0409351103},{0.8470924867,0.0225327257,0.8786042873,0.0315075989},{0.1136108062,0.1157168602,0.1151479227,0.7284565477},
	{0.2400109852,0.3520848716,0.4474164532,0.5024195253},{0.9868204636,0.2810682349,0.7242138409,0.1066441021},{0.5960941922,0.8775974365,0.3824219184,0.3542780208},
	{0.8545002637,0.2566439696,0.8116236029,0.8005692856},{0.0528010298,0.3028914153,0.1698182218,0.9578424878},{0.30687282,0.0163802253,0.5300991598,0.3148234084},
	{0.844853862,0.9723400001,0.9709698244,0.3836798114},{0.2996011388,0.3647791226,0.3967912565,0.804769835},{0.4388404957,0.0744674248,0.2168146712,0.5686755246},
	{0.3596344753,0.5750163288,0.827843966,0.7845171792},{0.2247795686,0.1692882418,0.8996852383,0.7610411681},{0.0615380472,0.3320678554,0.361293945,0.5495825282},
	{0.3598536041,0.3687895206,0.1466682267,0.4858825665},{0.1583273448,0.5435633028,0.842119056,0.8990433763},{0.2906893739,0.9992338961,0.0365902244,0.2515490062},
	{0.708700492,0.6068974198,0.2160222987,0.9092870934},{0.039755008,0.4774281885,0.8254781992,0.0544538819},{0.1710351154,0.9192502983,0.6021166332,0.2370720017},
	{0.5406602423,0.265420089,0.2110569589,0.8268453497},{0.8766309114,0.8712419572,0.4400364414,0.707734108},{0.4772796489,0.9419202821,0.0377061213,0.3069938442},
	{0.6366287735,0.5056157932,0.6075186774,0.0095591163},{0.9636435816,0.2276716309,0.7912478447,0.4505429803},{0.44003758,0.0816964463,0.7389633986,0.4780616283},
	{0.1236324313,0.3426056569,0.6533704132,0.716917925},{0.0089956555,0.4026452634,0.51142625,0.9524126311},{0.4314487611,0.8592355724,0.1587632687,0.3773702965},
	{0.1030082908,0.757335756,0.7177376123,0.5277655818},{0.3001537488,0.9666405551,0.6284318171,0.9788859088},{0.7184872362,0.411033236,0.3672082569,0.5313134938},
	{0.1407736824,0.686925529,0.0366793531,0.6887775373},{0.2582462602,0.8335112447,0.5239669241,0.7044642316},{0.0834274916,0.4183128239,0.9547965054,0.1445864537},
	{0.8125402953,0.194330411,0.0731127004,0.2454144256},{0.8005590437,0.3136733051,0.9915199603,0.4998389932},{0.757094119,0.2729283979,0.3822575053,0.7096403379},
	{0.229019637,0.0695607793,0.2090823546,0.2863386087},{0.9698746274,0.1490741572,0.5285111498,0.9226275445},{0.7017279847,0.8202711521,0.5225879175,0.9342761824},
	{0.504911415,0.0338309484,0.6523185109,0.2559272265},{0.0569188015,0.8451312648,0.7274005236,0.1772637055},{0.7712833968,0.252321901,0.0951467424,0.5321299012},
	{0.4056106803,0.5530326331,0.5664660218,0.3138524885},{0.6589405159,0.2750555868,0.5445726956,0.9707992866},{0.3487046889,0.2726835969,0.5993727229,0.2305923773},
	{0.7604714078,0.7029063362,0.3828110378,0.5409116901},{0.3817545557,0.1703660444,0.8216343466,0.7968471635},{0.281624191,0.506147047,0.5236769335,0.5702272665},
	{0.1489132931,0.6962481788,0.6821822962,0.5757386542},{0.2034474257,0.807975766,0.1412100238,0.6455826147},{0.7724329121,0.5498912982,0.2408004289,0.0549012059}
	};	
	
	
	std::vector<double> result = 
	{
	0.0866097,0.125324,0.159467,0.175506,0.188266,0.202629,0.218229,0.238788,0.258367,0.272335,
	0.287535,0.291796,0.302103,0.314302,0.328816,0.338456,0.359513,0.104478,0.116617,0.140207,
	0.154483,0.174251,0.19794,0.209221,0.227695,0.241821,0.255328,0.260923,0.279312,0.282396,
	0.300593,0.312468,0.326268,0.339791,0.0685421,0.110851,0.132986,0.151099,0.174988,0.188392,
	0.200116,0.215879,0.227822,0.236895,0.252518,0.265997,0.287909,0.296851,0.305085,0.310382,
	0.327059,0.0835115,0.116422,0.155861,0.180695,0.200998,0.207592,0.221693,0.241037,0.256161,
	0.268959,0.301436,0.316509,0.321373,0.328739,0.355433,0.362623,0.374216,0.0513198,0.0792087,
	0.134415,0.165734,0.184031,0.200341,0.208617,0.225573,0.244771,0.258952,0.271968,0.277888,
	0.286844,0.293381,0.305178,0.324316,0.342208,0.0724997,0.108502,0.121614,0.145896,0.154156,
	0.174737,0.184184,0.198854,0.212056,0.240021,0.244794,0.257261,0.283923,0.299432,0.303161,
	0.315952,0.326779,0.0770476,0.107026,0.122477,0.135401,0.154982,0.167743,0.194236,0.212707,
	0.221344,0.237592,0.248993,0.263536,0.274802,0.278451,0.289309,0.306458,0.329846,0.061784,
	0.119136,0.148535,0.164464,0.171897,0.194357,0.219205,0.238269,0.248606,0.253468,0.263734,
	0.268021,0.276119,0.284281,0.303345,0.331981,0.34947,0.0693335,0.0966073,0.114494,0.131379,
	0.152204,0.170513,0.19538,0.217221,0.223733,0.241719,0.269773,0.288419,0.302165,0.308418,
	0.3335,0.35408,0.365201,0.0790481,0.116622,0.128577,0.152609,0.166748,0.191723,0.197587,
	0.206596,0.231318,0.252562,0.265415,0.274281,0.28668,0.292747,0.302854,0.313834,0.32478,
	0.0917386,0.122776,0.13392,0.185631,0.196152,0.20375,0.218107,0.228042,0.234576,0.246406,
	0.25741,0.263371,0.279396,0.292965,0.322036,0.340515,0.35218,0.0783718,0.131059,0.157269,
	0.192033,0.206777,0.224203,0.240961,0.246639,0.264559,0.271092,0.278815,0.306389,0.318983,
	0.332846,0.342972,0.359236,0.368122,0.0956505,0.115942,0.146143,0.170414,0.176753,0.204491,
	0.218034,0.229521,0.242055,0.260202,0.262535,0.282024,0.288318,0.291844,0.303969,0.313169,
	0.322649,0.068779,0.0801762,0.129558,0.156775,0.173762,0.199728,0.213599,0.220743,0.235786,
	0.254517,0.267364,0.283058,0.305141,0.317227,0.319803,0.367818,0.382489,0.0610831,0.0922659,
	0.124723,0.153995,0.165965,0.218019,0.238821,0.250126,0.279189,0.288199,0.292874,0.307344,
	0.319509,0.333863,0.341756,0.345347,0.360549,0.0605367,0.108147,0.127604,0.136212,0.164292,
	0.200737,0.206311,0.220285,0.224466,0.246824,0.259217,0.267648,0.277226,0.290026,0.303104,
	0.312308,0.325163,0.0641294,0.0968253,0.129973,0.156516,0.16666,0.195934,0.205461,0.235305,
	0.246166,0.256025,0.262578,0.269232,0.289435,0.303173,0.316643,0.337324,0.343012,0.0612598,
	0.124711,0.168873,0.191569,0.205565,0.215769,0.22248,0.234508,0.243369,0.264701,0.279279,
	0.298818,0.310758,0.322746,0.340821,0.357411,0.372742,0.0957317,0.127598,0.153114,0.177751,
	0.192934,0.211943,0.222981,0.231499,0.247166,0.255986,0.261922,0.275213,0.295218,0.314247,
	0.326658,0.338533,0.355751,0.121453,0.150036,0.165379,0.18332,0.195748,0.214231,0.228823,
	0.251365,0.262165,0.265568,0.282378,0.29385,0.302009,0.317558,0.323095,0.352809,0.363861,
	0.0715727,0.122657,0.144742,0.159225,0.166199,0.196014,0.202384,0.224946,0.240898,0.26239,
	0.269074,0.282638,0.293598,0.297921,0.315352,0.327804,0.342365,0.0726949,0.112302,0.120957,
	0.15468,0.191265,0.214065,0.224953,0.238154,0.247458,0.25768,0.283219,0.296907,0.31284,
	0.318608,0.324269,0.343991,0.361648,0.0892491,0.122093,0.148321,0.187194,0.197991,0.208214,
	0.220901,0.227742,0.240416,0.250237,0.274866,0.288959,0.298797,0.326771,0.334757,0.340789,
	0.360899,0.0608986,0.0899089,0.127018,0.149539,0.169289,0.184124,0.208084,0.227737,0.238195,
	0.253864,0.267916,0.278967,0.288704,0.302056,0.313724,0.325431,0.340671,0.0880712,0.137694,
	0.15614,0.194605,0.208665,0.222524,0.237535,0.250299,0.263831,0.269494,0.275921,0.289822,
	0.301681,0.313393,0.321955,0.34341,0.350783,0.0815019,0.10582,0.150856,0.192749,0.197319,
	0.221869,0.231888,0.24647,0.263957,0.275979,0.282388,0.297067,0.311492,0.321032,0.330129,
	0.336379,0.348431,0.0882601,0.133602,0.143852,0.15695,0.178324,0.189598,0.20582,0.211729,
	0.232506,0.252056,0.266006,0.276394,0.290782,0.298351,0.317209,0.332221,0.348685,0.0914995,
	0.124842,0.152299,0.173437,0.192333,0.209723,0.22475,0.231043,0.2421,0.259526,0.274164,
	0.288814,0.307066,0.317301,0.325682,0.331918,0.346532,0.0741936,0.10266,0.123711,0.150315,
	0.18033,0.200827,0.217878,0.233333,0.242617,0.257915,0.277601,0.290513,0.295838,0.307851,
	0.323902,0.335572,0.349632,0.0795895,0.09272,0.107584,0.15596,0.183878,0.213134,0.231215,
	0.244493,0.257193,0.268749,0.286588,0.287938,0.293456,0.305612,0.321109,0.332637,0.351873,
	0.0408375,0.0826852,0.140027,0.161198,0.174718,0.198565,0.214005,0.227787,0.246141,0.267658,
	0.286096,0.292201,0.296739,0.311791,0.317159,0.329718,0.351254,0.0942555,0.137532,0.153218,
	0.164805,0.179295,0.209377,0.220127,0.233618,0.244996,0.256981,0.269251,0.276568,0.292969,
	0.308634,0.320397,0.334288,0.354947,0.0961236,0.135735,0.147026,0.164276,0.171454,0.192377,
	0.212539,0.225045,0.244028,0.257454,0.266532,0.273357,0.286257,0.301753,0.306998,0.323628,
	0.333987,0.0497792,0.0759842,0.118954,0.167147,0.188571,0.201003,0.223686,0.239499,0.249619,
	0.265449,0.278572,0.296129,0.309439,0.317208,0.33409,0.348384,0.365242,0.0391462,0.0900126,
	0.123074,0.143729,0.151276,0.164739,0.198676,0.217583,0.235066,0.265741,0.278013,0.292401,
	0.309013,0.322635,0.336504,0.347827,0.365956,0.0748703,0.143573,0.153131,0.167772,0.196315,
	0.211537,0.229611,0.236906,0.256497,0.270613,0.275117,0.28644,0.290044,0.316725,0.329796,
	0.338824,0.35595,0.085887,0.138247,0.155707,0.184998,0.19304,0.205036,0.226871,0.233774,
	0.248536,0.271459,0.279219,0.292937,0.29832,0.31757,0.33432,0.346091,0.360739,0.0928386,
	0.12704,0.138961,0.15889,0.174679,0.206007,0.21826,0.227383,0.240309,0.251842,0.270802,
	0.282871,0.297702,0.306995,0.321059,0.330588,0.342693,0.0663763,0.0859251,0.11808,0.16228,
	0.177882,0.199928,0.220991,0.240196,0.252693,0.263487,0.273606,0.284867,0.296141,0.317397,
	0.32462,0.348309,0.354792,0.0821784,0.119541,0.137463,0.155327,0.188503,0.202968,0.218645,
	0.234524,0.240787,0.244127,0.273046,0.280023,0.285655,0.300656,0.309052,0.323299,0.338604,
	0.0991852,0.12733,0.131931,0.144924,0.155779,0.169712,0.187696,0.196109,0.206093,0.224401,
	0.24703,0.264852,0.272773,0.286893,0.299019,0.309315,0.319852,0.0830694,0.111396,0.127359,
	0.147548,0.171328,0.207739,0.225412,0.249429,0.263104,0.276446,0.287832,0.310363,0.320783,
	0.328947,0.344325,0.34986,0.383047,0.0874361,0.122702,0.148753,0.176472,0.203901,0.2169,
	0.231378,0.245832,0.257375,0.264989,0.274648,0.282829,0.300013,0.308324,0.329985,0.333259,
	0.343251,0.0884037,0.121899,0.153379,0.179372,0.201228,0.228302,0.245323,0.253796,0.265163,
	0.279425,0.290852,0.308346,0.31918,0.333683,0.345473,0.35975,0.369886,0.0438193,0.0846735,
	0.112045,0.161229,0.200156,0.223314,0.230548,0.248852,0.258213,0.263938,0.270567,0.281429,
	0.295858,0.311941,0.322963,0.331735,0.33723,0.0898389,0.123435,0.13549,0.167141,0.193621,
	0.223138,0.245095,0.252018,0.261089,0.266054,0.274821,0.280248,0.291952,0.312171,0.319247,
	0.337419,0.340997,0.0757395,0.103078,0.12992,0.149682,0.167497,0.188656,0.210744,0.220612,
	0.236916,0.25018,0.260414,0.270016,0.283781,0.312304,0.322618,0.335296,0.352874,0.0646964,
	0.102234,0.13339,0.171403,0.1854,0.198939,0.212749,0.221991,0.229961,0.248018,0.256368,
	0.265681,0.277281,0.295467,0.311712,0.321506,0.339707,0.0614896,0.124465,0.146145,0.163679,
	0.183217,0.194214,0.225156,0.229859,0.233289,0.247264,0.262113,0.279266,0.284534,0.294258,
	0.306311,0.337329,0.342652,0.100996,0.136054,0.168452,0.177143,0.207277,0.233349,0.240245,
	0.247784,0.256557,0.264866,0.280801,0.295236,0.316373,0.335641,0.345588,0.356674,0.364501,
	0.113615,0.126265,0.137569,0.149641,0.167878,0.197135,0.216156,0.229975,0.236185,0.250339,
	0.259254,0.266057,0.278863,0.294072,0.303504,0.325873,0.339713,0.0622003,0.10953,0.131944,
	0.158866,0.170058,0.191859,0.206908,0.212745,0.22427,0.240502,0.256132,0.273809,0.294221,
	0.299023,0.312715,0.318466,0.328725,0.080588,0.127807,0.146645,0.171701,0.176168,0.198319,
	0.21268,0.229767,0.243313,0.260777,0.272898,0.293153,0.302359,0.312851,0.32506,0.329043,
	0.343662,0.0917355,0.107191,0.127893,0.175457,0.188276,0.203196,0.214822,0.235494,0.242891,
	0.271935,0.297286,0.304847,0.312831,0.32011,0.332003,0.3389,0.350738,0.0794702,0.134934,
	0.153001,0.188422,0.199888,0.218363,0.229559,0.240265,0.250144,0.269981,0.280699,0.304677,
	0.316568,0.330702,0.349111,0.355207,0.362931,0.0712038,0.131217,0.14486,0.156068,0.174391,
	0.198625,0.213459,0.220099,0.237021,0.248379,0.268483,0.274627,0.285641,0.298111,0.312112,
	0.331654,0.35263,0.0859684,0.108655,0.132032,0.154179,0.163348,0.188725,0.219401,0.231379,
	0.249295,0.262776,0.268252,0.291761,0.296496,0.308749,0.326946,0.332517,0.34272,0.0604923,
	0.0965592,0.118084,0.142133,0.170715,0.197492,0.21756,0.231135,0.237452,0.25676,0.265766,
	0.285056,0.300225,0.305928,0.329366,0.346901,0.356555,0.100459,0.12185,0.138151,0.16448,
	0.178396,0.19812,0.210938,0.217487,0.232068,0.244832,0.254612,0.278257,0.288287,0.302087,
	0.317163,0.335666,0.346261,0.0744918,0.137403,0.178545,0.20112,0.213921,0.229158,0.235943,
	0.249066,0.261907,0.270161,0.288761,0.303677,0.310612,0.327844,0.338527,0.350215,0.357326,
	0.0726899,0.114705,0.148934,0.174076,0.189532,0.213648,0.220725,0.232535,0.248909,0.265304,
	0.27706,0.281436,0.290919,0.307221,0.314537,0.327124,0.336118,0.0392844,0.109267,0.142175,
	0.147805,0.172161,0.192495,0.208501,0.249398,0.262134,0.277359,0.301428,0.314384,0.324151,
	0.336354,0.343642,0.355436,0.37063,0.0767362,0.0967328,0.122746,0.14268,0.160642,0.179364,
	0.189245,0.211358,0.218854,0.242638,0.256464,0.262695,0.266503,0.283954,0.304517,0.314057,
	0.327732,0.0772157,0.114239,0.125842,0.161884,0.178399,0.202112,0.219612,0.244349,0.265477,
	0.28005,0.293964,0.301829,0.32143,0.330037,0.338636,0.343926,0.360834,0.0850882,0.111993,
	0.148329,0.172057,0.189524,0.213236,0.22661,0.251363,0.263,0.277996,0.285864,0.301574,
	0.307207,0.324508,0.331168,0.345006,0.349815,0.0768453,0.106453,0.131103,0.1448,0.174137,
	0.202092,0.216306,0.221784,0.238828,0.247637,0.252568,0.268753,0.288929,0.293761,0.303363,
	0.306932,0.317041,0.0749116,0.104201,0.140335,0.167106,0.191231,0.205081,0.220393,0.22497,
	0.237117,0.255584,0.274871,0.293,0.308767,0.313705,0.325572,0.33344,0.345504,0.112425,
	0.145085,0.174019,0.21364,0.218983,0.229121,0.23907,0.258361,0.267786,0.281671,0.288335,
	0.300495,0.310846,0.333876,0.347286,0.353972,0.360197,0.114881,0.141465,0.168672,0.179058,
	0.189807,0.201754,0.213917,0.233079,0.239043,0.246625,0.256755,0.28447,0.297504,0.310206,
	0.322354,0.330466,0.346788,0.0668362,0.0941798,0.156804,0.180992,0.194353,0.205086,0.217438,
	0.240646,0.259589,0.268828,0.276649,0.282394,0.294358,0.30961,0.322069,0.333912,0.345657,
	0.0788634,0.107873,0.125384,0.160394,0.176808,0.195176,0.207734,0.218926,0.240365,0.260119,
	0.264112,0.289572,0.298589,0.307487,0.320737,0.330111,0.339929,0.0703792,0.112669,0.129513,
	0.150632,0.172436,0.184785,0.205699,0.220392,0.242551,0.254562,0.270966,0.283158,0.286289,
	0.315834,0.329296,0.344201,0.351778,0.0690751,0.114319,0.152391,0.166012,0.18931,0.20815,
	0.216895,0.230242,0.245848,0.263272,0.277203,0.299106,0.311008,0.32847,0.346993,0.355618,
	0.367335,0.0630264,0.118122,0.140714,0.164881,0.185818,0.19113,0.208594,0.237461,0.250509,
	0.263443,0.279506,0.284324,0.292117,0.319068,0.327881,0.336775,0.347627,0.0760772,0.122869,
	0.15636,0.179355,0.201149,0.212025,0.225639,0.242882,0.252217,0.265287,0.281016,0.289248,
	0.292147,0.298521,0.311098,0.324454,0.33762,0.0927185,0.119665,0.139867,0.164137,0.183787,
	0.19952,0.22158,0.23188,0.256252,0.263369,0.280217,0.293577,0.310009,0.316349,0.332449,
	0.339481,0.349364,0.0797347,0.118081,0.139748,0.153132,0.168745,0.183387,0.19849,0.204668,
	0.22541,0.237048,0.241634,0.256276,0.265419,0.281252,0.293998,0.300715,0.315721,0.0731127,
	0.111941,0.13677,0.163443,0.188737,0.197857,0.209361,0.225942,0.251,0.266321,0.27914,
	0.299717,0.317422,0.320806,0.330218,0.348269,0.359229,0.102284,0.119777,0.142937,0.15857,
	0.180401,0.217232,0.235617,0.25172,0.272078,0.287679,0.29556,0.305355,0.319068,0.330144,
	0.345311,0.350405,0.365992,0.0832736,0.121844,0.161151,0.182954,0.195955,0.2071,0.216066,
	0.222344,0.243316,0.261911,0.272174,0.28418,0.299429,0.316116,0.324338,0.345983,0.350908,
	0.0708657,0.0891118,0.12455,0.152897,0.171665,0.193708,0.223837,0.242997,0.249632,0.267399,
	0.293357,0.301859,0.310979,0.322583,0.335521,0.354081,0.370786,0.0884469,0.122106,0.167167,
	0.183091,0.194101,0.216149,0.230086,0.250652,0.270133,0.278237,0.284311,0.289981,0.299481,
	0.323254,0.337247,0.351064,0.369168,0.0599319,0.0875535,0.131651,0.15685,0.166842,0.179865,
	0.201158,0.218273,0.233975,0.257614,0.27169,0.279539,0.294794,0.301501,0.325893,0.33156,
	0.336832,0.0942224,0.118426,0.148402,0.174459,0.18667,0.199476,0.215677,0.229703,0.253635,
	0.27137,0.27869,0.295702,0.30324,0.309264,0.32085,0.333331,0.33955,0.0756793,0.112726,
	0.141026,0.165351,0.186983,0.210219,0.229916,0.250261,0.263287,0.272302,0.287889,0.296794,
	0.302923,0.31654,0.330686,0.345761,0.356344,0.0836171,0.11411,0.14007,0.168675,0.2076,
	0.222609,0.239705,0.245646,0.252493,0.264371,0.28003,0.29083,0.298695,0.307116,0.316116,
	0.329705,0.341579,0.0578763,0.115005,0.1409,0.165111,0.181799,0.200973,0.213138,0.225617,
	0.23601,0.256892,0.273024,0.291231,0.296728,0.306326,0.319511,0.33514,0.346151,0.101478,
	0.143856,0.16839,0.180834,0.215203,0.22731,0.240905,0.250313,0.260393,0.2641,0.284022,
	0.301775,0.309147,0.319524,0.337031,0.351068,0.367746,0.0863714,0.144618,0.181002,0.191477,
	0.207476,0.219968,0.230643,0.258247,0.278749,0.287187,0.295456,0.303397,0.315763,0.326518,
	0.335684,0.345821,0.354184,0.072736,0.130302,0.151094,0.172721,0.189127,0.213904,0.220457,
	0.237709,0.251279,0.267584,0.278921,0.293285,0.309073,0.322461,0.334776,0.346181,0.366824,
	0.0610155,0.0834556,0.101004,0.129648,0.175243,0.200372,0.208401,0.223487,0.240048,0.241287,
	0.253385,0.266108,0.275223,0.297898,0.3094,0.318562,0.347239,0.088783,0.122642,0.150914,
	0.178254,0.195109,0.211548,0.232988,0.237384,0.252011,0.26643,0.276273,0.2956,0.308303,
	0.322157,0.341267,0.353028,0.360762,0.0571079,0.112126,0.139146,0.161584,0.181456,0.191683,
	0.202602,0.235365,0.259493,0.270633,0.277032,0.293519,0.319303,0.331497,0.336299,0.343372,
	0.354646,0.0869447,0.11709,0.136713,0.154159,0.164471,0.189594,0.219064,0.231908,0.237789,
	0.256803,0.267706,0.274035,0.296756,0.31176,0.319487,0.330123,0.338544,0.0604749,0.127838,
	0.153142,0.165924,0.196884,0.21114,0.228728,0.243872,0.255087,0.266564,0.286064,0.289592,
	0.298367,0.314313,0.329439,0.351998,0.37327,0.0674806,0.0979469,0.127132,0.159232,0.167303,
	0.194463,0.210262,0.216675,0.238042,0.24526,0.257838,0.279103,0.289669,0.294223,0.298902,
	0.306579,0.316529,0.0774721,0.101642,0.124886,0.156422,0.187589,0.195694,0.210388,0.221973,
	0.227985,0.247221,0.257485,0.273696,0.282306,0.301981,0.314224,0.322761,0.331448,0.0925958,
	0.136807,0.169163,0.185664,0.195215,0.2112,0.225809,0.24335,0.250636,0.260341,0.286231,
	0.303753,0.31758,0.334553,0.35033,0.35641,0.373486,0.0586043,0.108019,0.124277,0.137655,
	0.16665,0.192384,0.207246,0.214843,0.23039,0.257324,0.264713,0.277807,0.280874,0.288152,
	0.297467,0.327447,0.333905,0.0971172,0.124306,0.14288,0.155493,0.16972,0.180381,0.198003,
	0.220574,0.233381,0.251591,0.268809,0.281429,0.287686,0.309842,0.31516,0.322481,0.33073,
	0.105097,0.126728,0.144192,0.175846,0.181042,0.202362,0.217117,0.237092,0.255407,0.263605,
	0.284819,0.296108,0.30284,0.315827,0.330997,0.344201,0.357134,0.0691672,0.104287,0.135275,
	0.168289,0.179279,0.188852,0.208246,0.215645,0.224918,0.236623,0.253709,0.269839,0.279728,
	0.299555,0.307055,0.33142,0.341521,0.0340118,0.13623,0.153594,0.165725,0.174385,0.190972,
	0.203867,0.223174,0.235996,0.253889,0.268241,0.285812,0.302672,0.316463,0.327054,0.339961,
	0.351844,0.0943352,0.115839,0.13336,0.162647,0.182941,0.191529,0.208967,0.219414,0.231203,
	0.245279,0.256459,0.272936,0.286895,0.292727,0.304123,0.328313,0.340691,0.0550216,0.110578,
	0.145821,0.166149,0.180678,0.2028,0.220144,0.240504,0.248595,0.261185,0.269743,0.283433,
	0.299313,0.302738,0.308811,0.325226,0.337577,0.120465,0.165244,0.168282,0.178406,0.189637,
	0.207124,0.220658,0.23599,0.249548,0.270185,0.276632,0.296876,0.313442,0.32621,0.335114,
	0.349813,0.366306,0.0636063,0.0997389,0.122595,0.151887,0.168329,0.196454,0.218305,0.227159,
	0.233997,0.24372,0.261896,0.276295,0.287373,0.292355,0.301144,0.310144,0.318296,0.0761112,
	0.115902,0.133989,0.150793,0.174443,0.19073,0.214835,0.234738,0.257601,0.266254,0.281151,
	0.296429,0.303394,0.326666,0.338239,0.344329,0.353129,0.0829953,0.113253,0.141536,0.161254,
	0.175881,0.194755,0.219118,0.227507,0.253025,0.261269,0.280762,0.293031,0.306985,0.309005,
	0.334526,0.342887,0.349116,0.0669628,0.0911037,0.108613,0.151058,0.17831,0.203073,0.234575,
	0.245622,0.252097,0.270338,0.285851,0.297895,0.309903,0.328699,0.346095,0.361071,0.375646,
	0.0742994,0.0911847,0.130471,0.172131,0.185628,0.208557,0.215593,0.223525,0.232455,0.243552,
	0.251835,0.260496,0.269395,0.276546,0.283039,0.299985,0.319843,0.0806465,0.104794,0.131362,
	0.14777,0.176275,0.2053,0.215809,0.226753,0.235692,0.244471,0.265634,0.273161,0.285807,
	0.293786,0.31317,0.329131,0.340454,0.0720657,0.125007,0.143192,0.157986,0.177412,0.193465,
	0.207258,0.226482,0.240785,0.265872,0.277719,0.281239,0.294046,0.300999,0.312408,0.319344,
	0.33827,0.0855427,0.117851,0.1597,0.174151,0.186415,0.208825,0.21485,0.240745,0.251324,
	0.263184,0.273975,0.290397,0.304306,0.31719,0.325133,0.348546,0.366005,0.0865814,0.107587,
	0.128564,0.144936,0.170277,0.180803,0.203571,0.219222,0.229753,0.23537,0.254037,0.268164,
	0.282044,0.306533,0.317725,0.329872,0.344548,0.0824779,0.128152,0.151422,0.171779,0.186215,
	0.213889,0.223892,0.24451,0.266028,0.27521,0.289762,0.294296,0.303598,0.3182,0.323909,
	0.343504,0.356392,0.0872621,0.115798,0.158029,0.161887,0.195026,0.213633,0.225745,0.232213,
	0.247867,0.260275,0.267328,0.279441,0.293108,0.308397,0.322491,0.333769,0.343857,0.048632,
	0.110428,0.142333,0.153489,0.166466,0.185041,0.197186,0.228679,0.253917,0.261973,0.273064,
	0.283013,0.29222,0.301292,0.308364,0.319135,0.324149,0.0713495,0.119152,0.132524,0.148131,
	0.161317,0.195961,0.206772,0.214898,0.229108,0.236556,0.257722,0.269625,0.27832,0.295376,
	0.312609,0.317899,0.332201,0.0885819,0.128667,0.160287,0.184927,0.197021,0.21447,0.223629,
	0.232221,0.249758,0.269039,0.281548,0.287262,0.295114,0.304859,0.324047,0.337538,0.366614,
	0.0738487,0.107811,0.131159,0.153789,0.164193,0.179242,0.19139,0.216516,0.226188,0.246075,
	0.263701,0.279958,0.29435,0.310588,0.331371,0.33784,0.362286,0.0752822,0.108601,0.14123,
	0.158362,0.182345,0.209964,0.214395,0.245556,0.267633,0.27568,0.294184,0.305045,0.322739,
	0.348439,0.358415,0.371849,0.386037,0.085353,0.131927,0.169685,0.17694,0.18574,0.195036,
	0.22927,0.240342,0.253878,0.266071,0.282581,0.301264,0.315549,0.325435,0.342666,0.354399,
	0.378471,0.0719061,0.0979673,0.127286,0.154186,0.159123,0.183349,0.19705,0.231169,0.237668,
	0.249969,0.264387,0.270985,0.294968,0.314376,0.317838,0.334414,0.345607,0.0858032,0.0978837,
	0.117975,0.15683,0.17609,0.196763,0.224105,0.233408,0.24756,0.257568,0.265986,0.282148,
	0.295271,0.303874,0.320034,0.325939,0.339522,0.0804856,0.106019,0.118565,0.136561,0.167915,
	0.180156,0.206613,0.226454,0.23552,0.261712,0.275658,0.287299,0.302089,0.312008,0.323509,
	0.331444,0.338995,0.0979951,0.127548,0.159247,0.182576,0.190378,0.206484,0.217948,0.224542,
	0.235132,0.247965,0.261603,0.277685,0.282364,0.294322,0.312957,0.331067,0.345518,0.0909583,
	0.116089,0.127927,0.173863,0.189903,0.200383,0.204621,0.224951,0.242193,0.257417,0.269236,
	0.276154,0.286126,0.310862,0.32457,0.334272,0.34765,0.0756548,0.101906,0.125566,0.163069,
	0.177698,0.186845,0.212928,0.232779,0.24488,0.25757,0.262906,0.285118,0.291149,0.299354,
	0.308125,0.328257,0.337387,0.0695973,0.0945282,0.12413,0.13665,0.163441,0.178353,0.193922,
	0.2146,0.22641,0.249348,0.256519,0.266675,0.282196,0.309411,0.318305,0.326,0.333543,
	0.0802617,0.116198,0.132775,0.172006,0.189522,0.199803,0.218355,0.226147,0.246105,0.25493,
	0.262287,0.268627,0.27468,0.288259,0.301248,0.314731,0.333959,0.0976088,0.116413,0.1371,
	0.162443,0.172151,0.200739,0.223263,0.242068,0.257668,0.270821,0.287911,0.295735,0.312491,
	0.322849,0.329495,0.341309,0.349157,0.140605,0.157256,0.181432,0.205278,0.221546,0.233098,
	0.237935,0.255805,0.265905,0.271016,0.27509,0.290212,0.298278,0.308765,0.335595,0.341964,
	0.357522,0.0650989,0.109293,0.126012,0.150571,0.201708,0.221332,0.242948,0.256907,0.264695,
	0.2785,0.289707,0.295681,0.311086,0.315735,0.341543,0.356904,0.381846,0.0666296,0.102536,
	0.130034,0.152498,0.172291,0.177294,0.203633,0.227108,0.240743,0.250142,0.257903,0.280665,
	0.29245,0.297218,0.304761,0.316983,0.330537,0.0759739,0.123526,0.152275,0.162357,0.180373,
	0.196809,0.203063,0.22699,0.238048,0.255064,0.263501,0.290106,0.30256,0.311425,0.331981,
	0.342248,0.355956,0.0929587,0.121635,0.135097,0.166535,0.192567,0.200047,0.215297,0.220377,
	0.230892,0.242567,0.265364,0.274777,0.287731,0.293379,0.313749,0.328445,0.343168,0.0695927,
	0.0885967,0.111911,0.143098,0.163242,0.181265,0.196372,0.215883,0.237569,0.25676,0.276172,
	0.284884,0.297006,0.303887,0.322557,0.336902,0.350339,0.0731276,0.108281,0.122628,0.141478,
	0.177175,0.205087,0.217407,0.240658,0.254903,0.266899,0.284993,0.301768,0.321723,0.337408,
	0.350338,0.364817,0.384992,0.0935971,0.105115,0.14644,0.170578,0.18814,0.219826,0.2352,
	0.262944,0.269516,0.280621,0.285375,0.294756,0.302603,0.315299,0.324174,0.34691,0.359126,
	0.0759081,0.116708,0.137456,0.166048,0.181472,0.199758,0.215416,0.222324,0.231229,0.254284,
	0.267913,0.293953,0.311728,0.324537,0.330346,0.347818,0.366289,0.0627016,0.108165,0.131921,
	0.165824,0.179121,0.196457,0.213357,0.228059,0.238969,0.259755,0.274208,0.281844,0.293119,
	0.307027,0.321349,0.337843,0.355811,0.0631621,0.105604,0.134326,0.153431,0.174078,0.191667,
	0.215092,0.223349,0.233288,0.240478,0.251417,0.262748,0.274282,0.284032,0.299785,0.323678,
	0.331263,0.0662261,0.094083,0.145618,0.166959,0.17664,0.186313,0.204052,0.229969,0.248943,
	0.265448,0.283015,0.304228,0.314081,0.322447,0.328711,0.343626,0.375055,0.0685608,0.0950153,
	0.133996,0.147783,0.171011,0.183074,0.199839,0.221071,0.231131,0.242096,0.251254,0.266953,
	0.290816,0.293721,0.307981,0.319196,0.324223,0.109615,0.121177,0.146341,0.160705,0.169396,
	0.188269,0.217159,0.230534,0.249595,0.261825,0.280702,0.301303,0.311353,0.326552,0.334763,
	0.350327,0.36518,0.0813862,0.101185,0.123248,0.164334,0.176006,0.208952,0.220644,0.230336,
	0.246282,0.261614,0.266886,0.272053,0.290027,0.301495,0.313338,0.325855,0.332749,0.10676,
	0.136148,0.151959,0.171667,0.180044,0.192493,0.202568,0.213818,0.2323,0.239598,0.266103,
	0.280737,0.290539,0.300109,0.314691,0.326531,0.340611,0.0641875,0.0982165,0.127956,0.141531,
	0.160579,0.174176,0.197402,0.203905,0.217523,0.231737,0.246484,0.264925,0.287368,0.299377,
	0.330067,0.338877,0.342968,0.0771866,0.100616,0.123925,0.150439,0.180829,0.214389,0.228747,
	0.236533,0.255469,0.266065,0.273667,0.282974,0.29297,0.315278,0.323968,0.333982,0.347466
	};
	
	
	//for every test point:
	unsigned counter = 0;
	for ( size_t pts = 0 ; pts != test_points.size() ; ++pts )
	{
		//for every distance class:
		for ( size_t dist = 0 ; dist != classes_to_compute_knn.size() ; ++dist )
		{
			//std::cout << (*classes_to_compute_knn[dist])( test_points[pts]) << ",";
			//if ( counter % 10 == 9 )  std::cout << std::endl;
			BOOST_CHECK( fabs( (*classes_to_compute_knn[dist])( test_points[pts]) - result[counter] ) <= 5e-07 );
			++counter;
		}
	}		
}













BOOST_AUTO_TEST_CASE(Distance_to_k_th_nearest_neighbor_periodic_domain_5d_brute_force)
{
	//first we test the brute force algorithm:

	typedef Gudhi::Topological_inference_with_cubical_complexes::Euclidan_distance_squared Euclidean_distance_squared;		
	typedef Gudhi::Topological_inference_with_cubical_complexes::periodic_domain_distance<Euclidean_distance_squared>
	periodic_Euclidean_distance_squared;
	
	
	std::vector< std::pair< double , double > > coordinates_of_grid(5);
	coordinates_of_grid[0] = std::pair< double,double >(0,1);
	coordinates_of_grid[1] = std::pair< double,double >(0,1);
	coordinates_of_grid[2] = std::pair< double,double >(0,1);
	coordinates_of_grid[3] = std::pair< double,double >(0,1);
	coordinates_of_grid[4] = std::pair< double,double >(0,1);
	
	Euclidean_distance_squared eu;
	periodic_Euclidean_distance_squared period_eu( coordinates_of_grid , eu );
	
	
	//150 random points from R^5.
	std::vector< std::vector<double> > point_cloud = 
	{
	{0.0834491136,0.1779654133,0.4660556072,0.9012534379,0.5278197562},{0.4226206406,0.9520258962,0.758530349,0.8480201229,0.9766116729},{0.5977266638,0.1149142182,0.5232524716,0.1020889988,0.5119932909},
	{0.2711605153,0.8722113823,0.1650474358,0.7590726558,0.6219455686},{0.2334552624,0.1426430494,0.3078001516,0.5127559074,0.7537710918},{0.5213451183,0.048260008,0.1171943273,0.8425379661,0.4401132523},
	{0.4934588228,0.9310606949,0.6055611819,0.8583622978,0.1591549711},{0.6904272269,0.3707491688,0.9402776032,0.2538115475,0.4506233654},{0.6859491069,0.0036608037,0.10070653,0.9493759079,0.467630194},
	{0.8316467148,0.5517980717,0.0530057724,0.541103472,0.6380860347},{0.5978446342,0.0052102418,0.1664767258,0.1295384325,0.5667432889},{0.1400902977,0.3020685925,0.9834268689,0.4287117575,0.7183386472},
	{0.7368947901,0.085055118,0.8047293324,0.3464451593,0.5628888928},{0.7584376174,0.6761060369,0.2481096517,0.4506944369,0.2866351088},{0.0174080082,0.5826915498,0.2013254741,0.0690666109,0.6637215966},
	{0.8963396226,0.4346969973,0.5366920386,0.1864755638,0.793938949},{0.4703953473,0.8998906589,0.2146033144,0.7245412155,0.5346764219},{0.8069159831,0.0554156622,0.435896911,0.5484183766,0.5096576582},
	{0.798520593,0.0641915055,0.3239880204,0.9386689118,0.0557170159},{0.530292036,0.2768214454,0.1465934187,0.2411683465,0.6128252845},{0.3530370982,0.1426775174,0.6023191165,0.561994422,0.4964279565},
	{0.324631419,0.7056863343,0.7003368,0.390044505,0.7481189952},{0.6095938848,0.2573646542,0.7254240075,0.8075247274,0.1642136844},{0.393126196,0.5145787923,0.094276082,0.7294868501,0.9212118047},
	{0.5256770824,0.9961912902,0.6040280494,0.7194164298,0.6646895669},{0.5757870409,0.9228846987,0.5187003482,0.5765369246,0.0548731403},{0.01892799,0.7966554496,0.4025369862,0.9533785735,0.657591725},
	{0.6390500749,0.2819597092,0.0462486721,0.3264356456,0.0678714598},{0.9890820112,0.6545676773,0.9895538145,0.6786556507,0.6729139483},{0.3993956982,0.230341801,0.6619966628,0.0914732311,0.8457622584},
	{0.1578961639,0.395656131,0.7078005325,0.31785637,0.799688499},{0.1970842553,0.2578141042,0.6121606366,0.895383196,0.0099467547},{0.4440435269,0.8098267191,0.3113613564,0.9574479449,0.590720393},
	{0.9892045073,0.6384359652,0.113328327,0.0851192193,0.2706557543},{0.4174217242,0.0269506972,0.3091899455,0.5228551389,0.9364830989},{0.3051545946,0.8718159075,0.8864578567,0.8777936108,0.2914440557},
	{0.5721974459,0.255364249,0.2955540344,0.2383654471,0.4111704903},{0.5746684466,0.8488304173,0.3052570869,0.5827505253,0.0440383067},{0.8940656511,0.6254615316,0.7737742225,0.7862985716,0.5748524205},
	{0.1149175356,0.1231047418,0.0073562327,0.1489178708,0.1661646403},{0.8169913939,0.0589694376,0.6719363823,0.706263863,0.6720062813},{0.1250505799,0.6658676548,0.7596249031,0.2439054565,0.8981688344},
	{0.7740507799,0.9550943908,0.5599841925,0.8664209775,0.2164539213},{0.7468787492,0.3491825978,0.3749921625,0.5969773906,0.5457275284},{0.5082658448,0.3032837098,0.6396805949,0.0767936446,0.6775225908},
	{0.3115240003,0.3442898365,0.4645590561,0.6598785846,0.0875235305},{0.5962873036,0.5788370646,0.866793555,0.8045623263,0.8786653674},{0.2912724817,0.1452073315,0.1192471741,0.3550970335,0.6413987856},
	{0.6897958643,0.7324675771,0.2894716042,0.1806752072,0.5615636073},{0.0443081178,0.4395110959,0.3479123823,0.4973605119,0.0052013861},{0.2848336392,0.6381691515,0.9037866145,0.5798630451,0.8606929414},
	{0.9962989811,0.6965007521,0.7684573706,0.6817337216,0.0348774616},{0.4469122239,0.2762809929,0.7385759256,0.0328260458,0.9326410054},{0.5957374696,0.1654470081,0.3903057303,0.081382965,0.8578045396},
	{0.2540016014,0.6546578028,0.8285822051,0.4711142455,0.5774126945},{0.4869810222,0.3995124626,0.6691552743,0.7868390677,0.8419497332},{0.3616856213,0.6670947331,0.0732981151,0.7414509852,0.5219372443},
	{0.2986048725,0.2752836884,0.9118515207,0.3655172742,0.1339146204},{0.9688093858,0.9839632511,0.0486876701,0.6362382285,0.57918212},{0.6837973881,0.9770498336,0.9009042082,0.9714206185,0.6519823587},
	{0.7795881715,0.0317756881,0.8326280683,0.8814528899,0.8136011944},{0.5725898354,0.2234702245,0.5974976968,0.4569877312,0.213307021},{0.7687702996,0.4392909671,0.8704199663,0.921301692,0.8155898936},
	{0.8652025636,0.4180430397,0.9787512729,0.6994813962,0.6923105665},{0.0511887514,0.4375224516,0.9453371717,0.0047599326,0.0097769545},{0.724004366,0.6219634346,0.1438660068,0.7190580415,0.5052032052},
	{0.4562779004,0.3600603256,0.9501040571,0.7044810567,0.8574787069},{0.8937145416,0.164020116,0.058241178,0.106201407,0.7696755296},{0.4888281797,0.766221792,0.936441195,0.062458213,0.0697754715},
	{0.9425974668,0.8562304049,0.7708064395,0.2568905996,0.6197091083},{0.335584769,0.1298263362,0.4283242184,0.8620519785,0.7777455822},{0.8629977736,0.8917749887,0.3128077728,0.0634608201,0.926209972},
	{0.7719056602,0.861920055,0.0497314511,0.6510132572,0.5403517021},{0.3530974325,0.0193289635,0.221168828,0.3852439316,0.0185912373},{0.0121225882,0.0861171798,0.5727441111,0.1638555448,0.6454807117},
	{0.8268414487,0.6949828749,0.8489062793,0.0978258594,0.0601115541},{0.1035729912,0.560648279,0.9121483676,0.2378044515,0.6134148247},{0.7951907262,0.9540303422,0.166627543,0.1205986922,0.240749388},
	{0.2657136505,0.4355481011,0.1904063555,0.574914173,0.8218553839},{0.2939299098,0.2204667355,0.9821966698,0.679141239,0.5713087774},{0.4602067114,0.4914364009,0.4681642714,0.1970379374,0.7063334736},
	{0.0375891852,0.0569641995,0.7455471558,0.2937743538,0.4793669109},{0.3807092302,0.5626776498,0.0227227262,0.3776322922,0.3102924996},{0.1073450749,0.1354784286,0.6111545078,0.0798917336,0.883374803},
	{0.8246012561,0.9048718265,0.0842031094,0.922818,0.1667780017},{0.3447176681,0.8322002466,0.6765862678,0.1630718487,0.4693403409},{0.1308114366,0.1080321325,0.1855266048,0.1161769815,0.547545319},
	{0.7396875334,0.4003160442,0.9041407737,0.4025102544,0.3810156607},{0.6799464941,0.1193180687,0.6461838551,0.203010628,0.9633627844},{0.6320310263,0.8174625379,0.2596935623,0.5561422193,0.5005034981},
	{0.2930128623,0.9211458443,0.8391498474,0.1270901843,0.6022140398},{0.0262923602,0.2269908616,0.423874124,0.6669442728,0.7403559329},{0.6004114933,0.0432240716,0.6893737733,0.7959725389,0.3219275989},
	{0.8859459318,0.6765051291,0.2261610297,0.9232959126,0.3750975686},{0.291880592,0.6774714992,0.4957603388,0.8096024855,0.9660385114},{0.0824616721,0.606248806,0.1800767379,0.0973097328,0.6039455077},
	{0.9669107667,0.8577861548,0.7781087584,0.9469715115,0.6491487199},{0.7758012626,0.4867858763,0.9445895557,0.1278092652,0.2407444329},{0.1551746733,0.9757911062,0.7507881073,0.1180643593,0.1485485181},
	{0.3832001633,0.3472629497,0.5424197274,0.2347819048,0.5460624669},{0.0356127378,0.3928180397,0.3578011249,0.4599245598,0.5875645839},{0.0915482768,0.8792861651,0.5032183493,0.4352760403,0.4537544101},
	{0.0198585391,0.2713452806,0.2556442935,0.1045154808,0.9693496483},{0.2023309676,0.9850838145,0.0787245242,0.7884351299,0.4749226195},{0.903540008,0.8357364505,0.0291661206,0.6786283345,0.6844293238},
	{0.6532926024,0.9651191607,0.0164092269,0.5536394559,0.20937135},{0.4534949544,0.2433945769,0.7650616835,0.8789803525,0.1802423056},{0.0479354102,0.5243833917,0.1338649723,0.3890133824,0.7519316291},
	{0.9194105798,0.1130673101,0.6573031303,0.8761956429,0.8871361001},{0.4819012247,0.0394615314,0.2906131288,0.5298820648,0.139980295},{0.9022430899,0.9868422102,0.7773927888,0.054481644,0.1240634422},
	{0.7944725945,0.8888423275,0.152707455,0.4338914778,0.9359518671},{0.1997224065,0.1918049152,0.1200931163,0.5627134435,0.752902451},{0.3878603228,0.8341636786,0.6820615984,0.2949337328,0.9880458449},
	{0.7945821437,0.5599427966,0.0107656661,0.3007575846,0.0132312805},{0.4316681421,0.2866483754,0.7283170826,0.9506913973,0.781435441},{0.4010647747,0.8797198473,0.1713188484,0.8730948223,0.4066781153},
	{0.7693045647,0.315029734,0.0817150525,0.7694409494,0.0741932406},{0.7115833587,0.2275419647,0.9498028783,0.1129892636,0.8711559433},{0.1673052812,0.8056030611,0.5349432651,0.2845978811,0.3322721664},
	{0.7902489975,0.2670829475,0.3567338747,0.7903393579,0.8519300122},{0.9817763602,0.5290622578,0.8170695412,0.1410530538,0.174742341},{0.5457802389,0.5899946659,0.0248836742,0.2630379079,0.9055795665},
	{0.6633629678,0.3279535186,0.8295126562,0.308469286,0.6487535806},{0.6586975821,0.0229490069,0.0493111487,0.8900570092,0.2659982985},{0.170515669,0.2176624145,0.3272540492,0.2422923984,0.5008245346},
	{0.5240949828,0.7313395606,0.8073427537,0.143435156,0.0210109602},{0.6617908955,0.9260073218,0.9820222741,0.8356835872,0.5093474893},{0.5819300218,0.2114108466,0.2355436478,0.4126429521,0.8040479794},
	{0.8829782095,0.0962476581,0.0428924561,0.0863333354,0.4965009892},{0.2279743734,0.8029537937,0.1028554635,0.6660963602,0.6046901501},{0.431163084,0.645723691,0.9890625442,0.7445918641,0.948092517},
	{0.0918637037,0.7363437202,0.930700952,0.3582612039,0.4598332951},{0.2691765854,0.1872886773,0.9881781302,0.3337876112,0.5859420574},{0.9778241748,0.9903361972,0.2781383186,0.3943724416,0.2876809081},
	{0.3436939148,0.829814577,0.8227764482,0.3955216904,0.6047144763},{0.8864226528,0.8145404733,0.9444013245,0.4600515633,0.7838328138},{0.711984874,0.5752081149,0.2798226338,0.9386363029,0.9097414555},
	{0.177447621,0.5581884892,0.7912682011,0.6638984308,0.5882288776},{0.5586297519,0.3793644228,0.3251419647,0.4594679768,0.1841192399},{0.6058601916,0.8781331608,0.4892053111,0.4373561433,0.3042695038},
	{0.052442651,0.3735250605,0.2064282887,0.8883306035,0.2979985315},{0.4378396841,0.4180250121,0.9897906918,0.2218166715,0.3110161072},{0.5262017339,0.6203828678,0.8560292989,0.8867694151,0.0002818264},
	{0.6696474531,0.3079569174,0.3395057747,0.2152758248,0.7913365713},{0.2165289272,0.3357743523,0.6355468514,0.2182823461,0.9831528417},{0.3455453599,0.4667630245,0.9083769515,0.7256064662,0.292692615},
	{0.2951296815,0.8800065387,0.8613584242,0.5988326801,0.2341764281},{0.0876358571,0.6400671883,0.0284326931,0.5182432309,0.0632727058},{0.1320868398,0.512312592,0.4902214163,0.6152805286,0.3113408913}
    };	
	
	std::vector< Gudhi::Topological_inference_with_cubical_complexes::Distance_to_k_th_closest_point<periodic_Euclidean_distance_squared>* >
	classes_to_compute_knn;	
	for ( size_t i = 5 ; i != 90 ; i=i+5 )
	{
		classes_to_compute_knn.push_back(
		new Gudhi::Topological_inference_with_cubical_complexes::Distance_to_k_th_closest_point<periodic_Euclidean_distance_squared>
		( point_cloud, period_eu , i ) );
	}  
	
	
	
	//150 random points to test the distance:
	std::vector< std::vector<double> > test_points =
	{
        {0.6156253421,0.0866876692,0.6546369994,0.7980181805,0.3036833229},{0.3805250507,0.5532490518,0.7087723003,0.8031594998,0.4108131663},{0.3181996481,0.8681068709,0.5779479195,0.1869301829,0.7758894945},
		{0.7863677139,0.3242641296,0.4091913456,0.5515330643,0.8118765294},{0.9259262064,0.4886787878,0.5708910702,0.64194575,0.8815862499},{0.7040850376,0.7995467752,0.5784777766,0.8222728742,0.62802323},
		{0.3415478489,0.1666721171,0.2690072777,0.1288172284,0.7955359623},{0.5544780709,0.6510210049,0.7063087481,0.5599074461,0.2960705815},{0.5480834248,0.7081139782,0.119110408,0.2898809859,0.5940711163},
		{0.403030508,0.4703414962,0.9957842538,0.7432314595,0.2526669232},{0.3263766619,0.6720029712,0.0591752841,0.3556350819,0.4690573232},{0.2064622687,0.5762362943,0.6930124161,0.8952994619,0.0152675807},
		{0.1605475396,0.0204712148,0.1160777733,0.247559764,0.8003176553},{0.6224387698,0.455551483,0.8776516232,0.3776043199,0.2195397872},{0.9540606579,0.1067592113,0.2777793263,0.1421314587,0.5134446463},
		{0.0406719297,0.6811534585,0.8866318022,0.4287882582,0.5006276881},{0.870287016,0.2856612143,0.1356323606,0.8231363078,0.335178002},{0.7992206207,0.2248386785,0.652452287,0.4542188642,0.1205765211},
		{0.5416745862,0.4066060479,0.8934159388,0.9497496316,0.7745553942},{0.5069284546,0.7684237019,0.5580311327,0.9252927683,0.7101417952},{0.0541560685,0.3102081812,0.5665760315,0.4848779158,0.047492543},
		{0.7745181869,0.5853341455,0.5212711582,0.7912951233,0.8805664999},{0.2830229441,0.032924467,0.8854985246,0.2631558124,0.6652764447},{0.6711386403,0.7294949514,0.078212844,0.8218062653,0.4102030019},
		{0.1846682311,0.2522296603,0.0682438477,0.3138512317,0.5468119136},{0.1068816297,0.8906914974,0.3227836362,0.9410896939,0.7084375294},{0.6661021973,0.2722590868,0.9509407876,0.3421589313,0.0829368159},
		{0.2404863781,0.9042210197,0.214130345,0.7180377392,0.5869330307},{0.4962969481,0.7940989283,0.9015325028,0.8254819594,0.1567669387},{0.7736787216,0.1744941329,0.021256563,0.6872016888,0.2290441715},
		{0.176827118,0.524485471,0.7884605827,0.9387880943,0.8023405797},{0.9890728234,0.6826795805,0.712908695,0.592580826,0.7100058955},{0.8278585384,0.3474175616,0.7225242155,0.4581460671,0.041690419},
		{0.3977338392,0.7302346188,0.6874800553,0.5466448646,0.0409471621},{0.2029349247,0.9303811735,0.3612236425,0.4019850001,0.7212867679},{0.0979896556,0.5700703752,0.8520184332,0.7930803259,0.1817028737},
		{0.4875848652,0.2855914799,0.1294624861,0.6426335031,0.2070947611},{0.9098975097,0.2684065762,0.4937620726,0.2055649131,0.5361320146},{0.4482507212,0.4339038245,0.4895149507,0.9056274116,0.6298295192},
		{0.3551141147,0.7516602229,0.0647157019,0.6934965344,0.9270651848},{0.9786594613,0.4435110099,0.3609527228,0.8899996865,0.3251161473},{0.153831613,0.4260658224,0.334146712,0.7129775598,0.853001711},
		{0.9990624285,0.3813134674,0.3835411465,0.098076452,0.4698095596},{0.9865754729,0.0481281504,0.9224064483,0.7275167136,0.2618650477},{0.0151539787,0.7301042816,0.1246206658,0.0785468051,0.6520250721},
		{0.1044145077,0.2527292348,0.3791852833,0.8313468685,0.6523172415},{0.677102146,0.5618522577,0.0822143261,0.3621671796,0.9020784714},{0.9044152494,0.487927774,0.4113524007,0.5217605885,0.4028691752},
		{0.9884912123,0.1933964493,0.8975249475,0.9582569813,0.0643633602},{0.6265765063,0.6546486341,0.6538947937,0.7916958395,0.2581015776},{0.1925964975,0.8758495869,0.7928991041,0.8404584907,0.0501623263},
		{0.5288344054,0.0940175443,0.3880097063,0.5240216907,0.2688257839},{0.151839297,0.4815489671,0.1200810564,0.5240271261,0.4553921684},{0.6643557919,0.7364383612,0.6186252583,0.9226146308,0.033159988},
		{0.7722630778,0.9259055518,0.8249877302,0.7974182731,0.7244633394},{0.1482382049,0.1489458592,0.760812395,0.4457491287,0.1411920183},{0.8833106682,0.1271617629,0.7506208313,0.9937048177,0.6780842175},
		{0.0076348747,0.578255268,0.6369468854,0.1451479641,0.5575621051},{0.4688744638,0.067973051,0.9584744396,0.4920711045,0.5704484717},{0.7850542171,0.0124715585,0.3579944114,0.5550846129,0.9165706104},
		{0.783577254,0.4116507513,0.7414916439,0.8980796521,0.4050635439},{0.1215615212,0.7875910646,0.141919357,0.5981509653,0.5387970002},{0.7539457558,0.4804520637,0.4587839125,0.2325882488,0.6723778283},
		{0.5750644468,0.3006371737,0.1820672811,0.0022302219,0.9543300846},{0.1849272884,0.3658491941,0.5859630774,0.2220380711,0.1284696721},{0.6724957384,0.4000277666,0.9780439236,0.0553804729,0.3112557554},
		{0.9480736756,0.5819517735,0.4126889601,0.2862954626,0.9022503761},{0.3177717375,0.2309198242,0.1540829227,0.7880887291,0.8507864},{0.7116010217,0.6213103663,0.2930113161,0.3727468587,0.9516895316},
		{0.3819570281,0.8180342214,0.2295541067,0.7412023214,0.8717434809},{0.7886661366,0.2495642588,0.048790921,0.1196943321,0.3640789606},{0.4985892794,0.0430643624,0.1038529542,0.9094649146,0.3767474603},
		{0.4307315028,0.4978282363,0.0282152942,0.2946561221,0.888703462},{0.6723954086,0.784596333,0.0863898948,0.7806103502,0.2525719581},{0.5714831352,0.4371523738,0.4938304564,0.8903684383,0.567739527},
		{0.3035034607,0.9636900444,0.9626759095,0.9843311792,0.2204016687},{0.4772450773,0.2324201055,0.7633715817,0.3358677777,0.9510003324},{0.3538140105,0.9728750708,0.5895036547,0.9765105813,0.6960510097},
		{0.857701164,0.8197515106,0.1978863701,0.0605267917,0.5568545174},{0.1860871471,0.5490134645,0.7524937373,0.3706550032,0.4462244327},{0.4302063505,0.2829570582,0.6203799653,0.6574771008,0.9612034687},
		{0.8402685211,0.9488183272,0.3951937449,0.683283203,0.6873453115},{0.5700790631,0.1601365875,0.7559559469,0.1201780995,0.2192969187},{0.5646118473,0.1211143576,0.7025673755,0.6807703031,0.0820197742},
		{0.7155782501,0.3898127335,0.8637648029,0.2135735601,0.130055618},{0.5476660505,0.0127128572,0.0949698025,0.1445218846,0.094173088},{0.7984791892,0.6017809003,0.4891143411,0.4428788116,0.5716396598},
		{0.2536687476,0.3673323095,0.6265329041,0.6520072613,0.8490074677},{0.7154158063,0.3787585783,0.1292777357,0.3725556275,0.3312357822},{0.4494496752,0.0517784276,0.2138768274,0.3533489655,0.6084807601},
		{0.2920235356,0.1409710939,0.1904358806,0.4535624995,0.9371193126},{0.5526888377,0.5882126573,0.827937352,0.858154837,0.7465898003},{0.5247145414,0.6173376355,0.8801819091,0.2626326918,0.0044825338},
		{0.3793205258,0.2487422219,0.5722060839,0.229396201,0.0556458414},{0.3603424882,0.3157115481,0.4432009482,0.318166235,0.0832488397},{0.1691172069,0.4944495703,0.1648998465,0.7261980609,0.8623710414},
		{0.3742564206,0.7017267803,0.7860315635,0.5872677292,0.7537247543},{0.7428908623,0.0345434777,0.9549517429,0.8074096772,0.849237503},{0.7985906939,0.0946825973,0.3014169701,0.8058131211,0.1560831016},
		{0.2603009364,0.1148938027,0.8615724845,0.4697891294,0.9368353917},{0.5338982816,0.2307557308,0.6977123844,0.2302658497,0.1019402801},{0.9767520463,0.718872142,0.6688464943,0.4160435961,0.9022799341},
		{0.0166666615,0.4389940493,0.7579744349,0.5370429766,0.1495361354},{0.3979732869,0.4755268127,0.5795001904,0.9908350941,0.1669989298},{0.8949154501,0.2682687736,0.6063345135,0.2275951863,0.017111354},
		{0.4043473324,0.3614203415,0.5158000104,0.1085771453,0.8999612131},{0.977441187,0.1501589406,0.8518329735,0.3219707806,0.8540179215},{0.4475535215,0.0862195999,0.8618425971,0.9331469764,0.4519766094},
		{0.4441091421,0.9066383392,0.0553346425,0.5744697268,0.9810071255},{0.9267698682,0.9928793737,0.5727517998,0.1823302133,0.5895151296},{0.7495447183,0.6460728368,0.470270294,0.6661383023,0.0136451058},
		{0.6532815052,0.325947958,0.675558131,0.2983020572,0.5749601708},{0.59615145,0.4472527723,0.581177111,0.5676356531,0.77461539},{0.7588686878,0.4715215371,0.7138359281,0.6454586445,0.6092199853},
		{0.6003141189,0.7735408389,0.7334262682,0.5612569719,0.5469313441},{0.8622437229,0.4840239387,0.3084369798,0.2880281513,0.3305151053},{0.7834043312,0.2515344019,0.8712158261,0.3956733951,0.6839187194},
		{0.7728605138,0.4649363579,0.1062444423,0.2089061183,0.2042179685},{0.2607685262,0.5400420581,0.7663626072,0.5831070796,0.1841536507},{0.4371812299,0.9649007791,0.4027777456,0.0520443097,0.6338951529},
		{0.5583380917,0.7672120824,0.3233393638,0.3555802056,0.5658726953},{0.5354220215,0.7534636182,0.9579864799,0.6235753896,0.9860881409},{0.9109844903,0.9209281902,0.2032730938,0.0754138096,0.4649585146},
		{0.20800173,0.8517249846,0.9164571092,0.2328442468,0.1823428017},{0.2348916477,0.4451530215,0.5246549759,0.7947810211,0.6921255561},{0.6184620375,0.1530659061,0.9847622805,0.7646336413,0.6839535253},
		{0.4769351536,0.1649140997,0.5587477549,0.4871179559,0.1689623857},{0.0984163808,0.6076200083,0.966948214,0.1425744521,0.2599090473},{0.6092678127,0.14414063,0.6357082536,0.4362932767,0.9267203519},
		{0.2324069291,0.3263996779,0.682590036,0.591015572,0.0428221659},{0.2328919342,0.7355384333,0.8440155704,0.7060465526,0.3706862081},{0.8768993784,0.8748575826,0.4019506322,0.9472513744,0.7443623699},
		{0.0835958309,0.2113429978,0.4997533069,0.3181976748,0.5457607873},{0.4077742845,0.1392330623,0.8563316758,0.1835489655,0.016605665},{0.9483211499,0.9293324859,0.939470334,0.891758071,0.5508040495},
		{0.1962504159,0.1287407137,0.7581467074,0.7718736557,0.280741273},{0.4436504522,0.8687367935,0.5447200509,0.848678587,0.4843223449},{0.5696217748,0.9954605456,0.6056541163,0.2167235308,0.562121178},
		{0.9663550276,0.1291832027,0.860948727,0.2084195951,0.1163226201},{0.1134767649,0.3865185466,0.9239129042,0.7176613885,0.3527642083},{0.3190164822,0.7672213989,0.2879216047,0.6133680425,0.4904669598},
		{0.2613077373,0.1525334122,0.048044354,0.678694438,0.5231821293},{0.1799443669,0.7424147304,0.896189518,0.6574717208,0.9911259802},{0.043235773,0.0527958591,0.0389118574,0.0312990851,0.419792705},
		{0.1913606562,0.6254830866,0.5446361115,0.6739698576,0.0345029861},{0.9619858956,0.8867313419,0.8675114349,0.7746892928,0.0028468464},{0.4198780232,0.6695927838,0.4730970562,0.8541963291,0.0656680651},
		{0.4705536903,0.2660579658,0.9417342963,0.9129116719,0.1820742511},{0.2058038623,0.0234297065,0.9110244056,0.4768905714,0.1417649831},{0.2174516551,0.0631806597,0.2748815217,0.5283225405,0.3984943791}
	};	
	
	
	std::vector<double> result = 
	{
	0.0848295,0.167682,0.190143,0.224555,0.260605,0.276135,0.299343,0.308722,0.325974,0.332576,
	0.355922,0.363254,0.382531,0.395273,0.411649,0.428977,0.445709,0.162351,0.213094,0.233111,
	0.273276,0.291066,0.299465,0.307219,0.316869,0.329778,0.349699,0.363402,0.369242,0.379291,
	0.393963,0.399038,0.416935,0.437603,0.129412,0.175697,0.205706,0.222322,0.237218,0.257294,
	0.281272,0.290291,0.30679,0.32004,0.334267,0.356391,0.375863,0.399707,0.420824,0.42818,
	0.437744,0.123879,0.183859,0.223911,0.250866,0.272731,0.284377,0.29067,0.316663,0.330722,
	0.346474,0.360469,0.373744,0.384353,0.390481,0.416428,0.424699,0.435283,0.17581,0.203497,
	0.220094,0.234968,0.271108,0.301817,0.316558,0.327927,0.34733,0.360667,0.371183,0.385407,
	0.396205,0.410337,0.442214,0.456112,0.460307,0.148153,0.185379,0.206594,0.246517,0.264187,
	0.281161,0.293867,0.312138,0.321497,0.331995,0.347012,0.355599,0.36766,0.375588,0.405274,
	0.415481,0.440765,0.116477,0.145421,0.171826,0.210955,0.23101,0.248492,0.271623,0.302032,
	0.308633,0.314454,0.327844,0.339908,0.353454,0.367661,0.385579,0.404055,0.420144,0.168283,
	0.212246,0.237069,0.250381,0.264192,0.266797,0.2842,0.291756,0.305796,0.324195,0.346134,
	0.365502,0.377756,0.394817,0.413001,0.428632,0.442215,0.146687,0.187936,0.211868,0.235742,
	0.247644,0.265255,0.272803,0.282406,0.296366,0.310016,0.326382,0.334602,0.347457,0.370799,
	0.38372,0.394413,0.399228,0.130958,0.180165,0.212301,0.226185,0.251889,0.278232,0.292023,
	0.31947,0.334075,0.342909,0.362938,0.37578,0.377176,0.386597,0.407513,0.414678,0.438542,
	0.118381,0.18156,0.209455,0.233214,0.244773,0.256898,0.269794,0.29571,0.301167,0.306994,
	0.322214,0.330463,0.341341,0.35721,0.385369,0.406063,0.417629,0.131054,0.17019,0.193571,
	0.209064,0.233001,0.26,0.275177,0.305442,0.315445,0.343536,0.362786,0.381563,0.398439,
	0.428502,0.440018,0.444011,0.461674,0.116072,0.15786,0.197647,0.238552,0.25556,0.265478,
	0.28922,0.296169,0.30705,0.316529,0.32188,0.333712,0.347416,0.36699,0.375887,0.383202,
	0.40527,0.0918235,0.157325,0.209278,0.256533,0.279914,0.293528,0.31449,0.318986,0.34049,
	0.367845,0.375007,0.388779,0.39476,0.406577,0.418407,0.425026,0.43391,0.115489,0.153122,
	0.20137,0.214372,0.242651,0.259452,0.271232,0.283506,0.308844,0.316242,0.32176,0.346473,
	0.357072,0.371264,0.387054,0.415065,0.425979,0.105873,0.15049,0.178694,0.19109,0.214227,
	0.224868,0.254229,0.269422,0.297074,0.31641,0.340962,0.350929,0.357503,0.379616,0.406447,
	0.441294,0.460655,0.148205,0.184997,0.209429,0.230873,0.252504,0.276562,0.290241,0.299428,
	0.306483,0.318678,0.339948,0.359892,0.376597,0.394306,0.411601,0.429351,0.442051,0.178267,
	0.23121,0.248141,0.261049,0.274753,0.287126,0.301017,0.316553,0.331258,0.346578,0.356343,
	0.371992,0.391065,0.39716,0.405855,0.42154,0.433167,0.0818327,0.130023,0.181482,0.218445,
	0.227169,0.260478,0.276303,0.290187,0.302683,0.314346,0.317565,0.327318,0.350202,0.364067,
	0.378657,0.399229,0.41126,0.158959,0.194153,0.200482,0.239815,0.256202,0.270617,0.279928,
	0.291677,0.324638,0.335548,0.344593,0.353398,0.374008,0.389338,0.407072,0.429214,0.439384,
	0.13937,0.214241,0.249153,0.262912,0.300738,0.319474,0.327932,0.339264,0.345877,0.360758,
	0.377875,0.389604,0.395612,0.412283,0.431313,0.437804,0.459013,0.158461,0.201442,0.229945,
	0.256686,0.2872,0.299562,0.309594,0.321666,0.335328,0.340489,0.347109,0.369803,0.381844,
	0.400009,0.407782,0.423448,0.433999,0.115902,0.164519,0.192033,0.206392,0.23014,0.241564,
	0.25929,0.269414,0.289775,0.299306,0.314701,0.338119,0.347042,0.37126,0.380735,0.38967,
	0.410681,0.0954641,0.120935,0.171904,0.192329,0.208385,0.236661,0.264469,0.290409,0.305562,
	0.32002,0.329735,0.341951,0.362225,0.378627,0.399776,0.412042,0.426937,0.0765244,0.148824,
	0.170286,0.1989,0.225955,0.271435,0.282216,0.299209,0.314327,0.332349,0.349185,0.360365,
	0.372509,0.375316,0.381403,0.394748,0.39824,0.13168,0.155418,0.200041,0.221016,0.244805,
	0.26525,0.291607,0.315401,0.323612,0.342006,0.360835,0.368547,0.37463,0.388037,0.406282,
	0.41778,0.424659,0.128945,0.158617,0.207504,0.232391,0.241996,0.275163,0.291807,0.304589,
	0.328356,0.348878,0.354584,0.361452,0.366657,0.376261,0.388253,0.401706,0.42167,0.0847531,
	0.144881,0.163493,0.204349,0.24786,0.2631,0.274815,0.285477,0.305458,0.322563,0.328079,
	0.351235,0.368374,0.382791,0.397444,0.416454,0.430871,0.0840074,0.135321,0.162997,0.205367,
	0.241821,0.258843,0.272941,0.299455,0.311363,0.333512,0.347029,0.363076,0.379311,0.385606,
	0.402224,0.416354,0.42455,0.139989,0.177858,0.200324,0.249376,0.266311,0.283804,0.293187,
	0.299576,0.321155,0.336324,0.344766,0.356998,0.369519,0.379655,0.396388,0.404356,0.416574,
	0.147076,0.165801,0.177249,0.18273,0.214025,0.235448,0.259896,0.269111,0.296584,0.311088,
	0.317151,0.33647,0.36497,0.372285,0.394458,0.414546,0.4296,0.116688,0.156473,0.178371,
	0.209394,0.243135,0.258517,0.285988,0.311073,0.322853,0.34883,0.364372,0.385156,0.396268,
	0.403747,0.415583,0.43554,0.453863,0.162764,0.198535,0.227761,0.257608,0.273038,0.284108,
	0.29944,0.314218,0.337691,0.346333,0.356919,0.373788,0.384479,0.393329,0.408232,0.415878,
	0.443275,0.116384,0.174319,0.193311,0.221845,0.24301,0.265486,0.280428,0.3086,0.332731,
	0.353183,0.385539,0.400152,0.416846,0.430188,0.446937,0.459844,0.484937,0.138746,0.179503,
	0.21213,0.241658,0.249752,0.261302,0.276726,0.286747,0.295024,0.310876,0.322534,0.3324,
	0.349908,0.356612,0.382253,0.400623,0.419267,0.137542,0.188908,0.204822,0.23484,0.257798,
	0.289237,0.300299,0.329093,0.340438,0.346706,0.352856,0.379992,0.38789,0.409263,0.425078,
	0.437179,0.444264,0.149244,0.164752,0.192312,0.227408,0.248319,0.267239,0.28163,0.294676,
	0.316238,0.329648,0.339191,0.349276,0.369222,0.391152,0.408676,0.420909,0.426995,0.131742,
	0.177568,0.20558,0.233554,0.258632,0.271847,0.281661,0.296821,0.313477,0.332356,0.34915,
	0.352513,0.369843,0.376887,0.389902,0.404039,0.439821,0.129897,0.186526,0.21181,0.243386,
	0.267813,0.28169,0.296669,0.312774,0.326523,0.3337,0.350471,0.366826,0.385017,0.399471,
	0.419172,0.423575,0.438258,0.129045,0.148848,0.181683,0.213869,0.240414,0.259456,0.274177,
	0.292159,0.314131,0.325991,0.33858,0.357101,0.369324,0.376263,0.389903,0.418074,0.430231,
	0.140462,0.215015,0.244659,0.269943,0.298829,0.312447,0.336464,0.347849,0.358189,0.370711,
	0.386266,0.405851,0.418113,0.428324,0.437489,0.446205,0.459274,0.12756,0.157548,0.223411,
	0.235681,0.256191,0.271037,0.276026,0.298884,0.320632,0.349332,0.360572,0.379637,0.392041,
	0.408288,0.41877,0.433888,0.443498,0.116942,0.164221,0.229464,0.25412,0.285386,0.298987,
	0.316394,0.327737,0.341213,0.357346,0.370233,0.387133,0.398025,0.412334,0.421238,0.435437,
	0.441612,0.150679,0.180801,0.211745,0.228479,0.255588,0.272247,0.291648,0.304284,0.310504,
	0.328107,0.337472,0.36513,0.383238,0.389674,0.404231,0.41587,0.426158,0.13067,0.158706,
	0.184978,0.219099,0.234467,0.252476,0.270417,0.27962,0.296749,0.320483,0.329018,0.341914,
	0.356072,0.378279,0.384539,0.390046,0.417365,0.145659,0.182728,0.197889,0.220565,0.231055,
	0.247389,0.27591,0.300598,0.32187,0.329434,0.340978,0.347406,0.356534,0.370672,0.388466,
	0.420762,0.442789,0.131959,0.185812,0.202912,0.218312,0.226104,0.249487,0.270124,0.286387,
	0.297756,0.313003,0.316846,0.334897,0.344937,0.362273,0.399007,0.415054,0.423063,0.171461,
	0.20929,0.238393,0.258245,0.29075,0.310392,0.326796,0.35551,0.360472,0.3673,0.377338,
	0.383488,0.398097,0.417919,0.433751,0.447154,0.46545,0.129293,0.144404,0.17273,0.215114,
	0.248272,0.277652,0.299891,0.30916,0.319353,0.335283,0.358673,0.375009,0.398405,0.412348,
	0.418744,0.426281,0.442124,0.172305,0.205185,0.220035,0.247603,0.277901,0.296177,0.302768,
	0.312524,0.325525,0.355209,0.368689,0.382312,0.392242,0.412473,0.417307,0.428925,0.44487,
	0.107461,0.170044,0.192852,0.226081,0.262557,0.276851,0.2905,0.300934,0.324338,0.345697,
	0.353681,0.377779,0.390725,0.409954,0.420187,0.422235,0.441866,0.0976053,0.139708,0.175301,
	0.228986,0.251131,0.270622,0.292676,0.303091,0.326519,0.343995,0.356215,0.369919,0.389538,
	0.398167,0.420586,0.425948,0.443788,0.13188,0.145626,0.158717,0.191441,0.228119,0.245914,
	0.271726,0.285823,0.314104,0.319178,0.331914,0.3475,0.36598,0.386244,0.403841,0.421966,
	0.435905,0.112575,0.170152,0.1962,0.230681,0.250902,0.279758,0.296844,0.322884,0.345273,
	0.355072,0.36449,0.377843,0.398354,0.410278,0.428198,0.439031,0.453815,0.0827667,0.130232,
	0.178571,0.224195,0.250608,0.268633,0.284809,0.288505,0.306969,0.318577,0.333351,0.340429,
	0.353532,0.373966,0.387767,0.398494,0.412389,0.151294,0.208907,0.24045,0.269328,0.291335,
	0.310029,0.317076,0.327049,0.339157,0.357268,0.361326,0.375106,0.383193,0.392716,0.400556,
	0.413722,0.423666,0.0861021,0.130745,0.163834,0.194855,0.241684,0.251902,0.256774,0.282805,
	0.307534,0.313932,0.335699,0.351102,0.363955,0.377639,0.387762,0.409197,0.424532,0.149395,
	0.173346,0.210495,0.24222,0.255553,0.270058,0.288912,0.323015,0.348704,0.362597,0.379561,
	0.389629,0.399567,0.4079,0.41985,0.432498,0.441143,0.117028,0.161888,0.180571,0.194028,
	0.212216,0.244448,0.260669,0.280721,0.294021,0.297394,0.312806,0.325237,0.347229,0.350885,
	0.357319,0.370377,0.387034,0.128781,0.174014,0.222018,0.232126,0.247162,0.264472,0.285935,
	0.324743,0.338422,0.349808,0.367192,0.375781,0.385897,0.399466,0.414765,0.427602,0.438449,
	0.178487,0.204078,0.25002,0.265134,0.28147,0.286906,0.296593,0.311568,0.319806,0.345901,
	0.354039,0.373016,0.383799,0.393364,0.398804,0.410475,0.427772,0.0898087,0.151021,0.174144,
	0.201235,0.218836,0.232794,0.25022,0.26499,0.293729,0.302069,0.327005,0.337045,0.364822,
	0.386036,0.397253,0.416544,0.444209,0.148769,0.178147,0.201566,0.251541,0.277521,0.288767,
	0.298912,0.308631,0.323918,0.337714,0.358439,0.37182,0.382773,0.388008,0.42289,0.434889,
	0.451944,0.116565,0.169483,0.198001,0.222041,0.24346,0.264877,0.280762,0.293549,0.304554,
	0.319997,0.346035,0.350225,0.364877,0.386287,0.397278,0.400814,0.412414,0.140001,0.177271,
	0.2161,0.239569,0.270861,0.301758,0.320981,0.329298,0.340251,0.370782,0.374977,0.386382,
	0.397013,0.410634,0.417954,0.439699,0.465427,0.152417,0.175266,0.190577,0.20998,0.227925,
	0.258756,0.283242,0.287372,0.308132,0.318032,0.331809,0.347206,0.365245,0.395184,0.418917,
	0.429027,0.446135,0.160541,0.177337,0.208286,0.261503,0.278546,0.307053,0.320308,0.340532,
	0.351113,0.358767,0.371271,0.390767,0.399866,0.407185,0.415349,0.429543,0.448921,0.0981115,
	0.184175,0.220715,0.232539,0.250641,0.265616,0.282201,0.302805,0.30902,0.31595,0.328375,
	0.335658,0.348514,0.354119,0.368222,0.380623,0.390701,0.125477,0.188011,0.201934,0.244118,
	0.280797,0.294112,0.301579,0.309908,0.320634,0.34168,0.364519,0.383665,0.399118,0.403904,
	0.411475,0.425422,0.439594,0.112308,0.162315,0.191481,0.235512,0.242443,0.263878,0.291495,
	0.307304,0.324999,0.32957,0.34535,0.35428,0.384397,0.394822,0.405313,0.416844,0.423492,
	0.116479,0.153696,0.180251,0.19926,0.21744,0.252754,0.261896,0.290729,0.312526,0.333502,
	0.342694,0.35152,0.364659,0.382995,0.395576,0.406227,0.422965,0.0782063,0.148522,0.17906,
	0.199923,0.223234,0.242778,0.25918,0.294968,0.310624,0.323587,0.344621,0.360502,0.370311,
	0.386531,0.399036,0.406945,0.420366,0.145674,0.170536,0.19136,0.196193,0.21689,0.232604,
	0.245754,0.27559,0.287335,0.309016,0.315254,0.34633,0.358415,0.371508,0.381914,0.400031,
	0.411754,0.100028,0.13229,0.165956,0.193402,0.225949,0.244173,0.26185,0.276087,0.286655,
	0.299438,0.311754,0.339737,0.35866,0.383074,0.396801,0.410824,0.449413,0.146488,0.199293,
	0.217977,0.24622,0.263588,0.279784,0.314182,0.326571,0.346735,0.353918,0.3648,0.381911,
	0.401844,0.407745,0.421561,0.436449,0.455741,0.107174,0.152495,0.194754,0.23575,0.253101,
	0.262348,0.282684,0.300142,0.310859,0.32947,0.338408,0.355884,0.369004,0.378954,0.385853,
	0.397568,0.41235,0.109854,0.154938,0.208312,0.233663,0.252158,0.262422,0.280651,0.289883,
	0.301815,0.313918,0.32318,0.342187,0.350699,0.361938,0.380097,0.396628,0.418188,0.113658,
	0.133705,0.180173,0.21138,0.248999,0.258898,0.271047,0.280221,0.297073,0.320469,0.340618,
	0.351849,0.358379,0.386239,0.396668,0.40315,0.412565,0.093209,0.139833,0.179182,0.199063,
	0.23301,0.24386,0.269891,0.277405,0.302523,0.310307,0.333773,0.34459,0.354937,0.372694,
	0.394009,0.405079,0.424012,0.129614,0.152117,0.184388,0.214176,0.249647,0.25893,0.283146,
	0.29842,0.319103,0.339998,0.359316,0.382936,0.392815,0.409367,0.426785,0.436467,0.462834,
	0.120082,0.155997,0.201851,0.217403,0.235926,0.251093,0.275238,0.288446,0.307993,0.319555,
	0.339172,0.359364,0.378091,0.402328,0.415955,0.434145,0.449492,0.130052,0.154217,0.196823,
	0.225966,0.248047,0.264181,0.276009,0.284017,0.297434,0.31933,0.346547,0.353876,0.375991,
	0.391038,0.415219,0.42397,0.433356,0.134664,0.157415,0.183242,0.205096,0.228939,0.247606,
	0.262267,0.271826,0.287477,0.311757,0.325264,0.337299,0.350165,0.389189,0.398486,0.41303,
	0.437195,0.0880821,0.163076,0.185834,0.21319,0.237257,0.253184,0.268479,0.294161,0.32533,
	0.334532,0.36589,0.384427,0.395402,0.407877,0.426187,0.447139,0.457298,0.101,0.1477,
	0.185195,0.21434,0.248383,0.261938,0.300237,0.309411,0.32948,0.343705,0.360918,0.375282,
	0.38906,0.405813,0.415136,0.428248,0.451173,0.117487,0.172719,0.200242,0.21634,0.23322,
	0.245205,0.261798,0.282962,0.30402,0.318726,0.332854,0.34158,0.363101,0.372751,0.385748,
	0.399648,0.408507,0.146499,0.205623,0.223668,0.262026,0.284331,0.305608,0.318508,0.328317,
	0.337227,0.351808,0.36087,0.387653,0.398202,0.412079,0.429956,0.445785,0.461238,0.130668,
	0.15328,0.201685,0.225663,0.234401,0.255873,0.269483,0.290109,0.300037,0.316744,0.344256,
	0.372011,0.388302,0.416099,0.424778,0.453009,0.473605,0.0935834,0.157673,0.209764,0.235338,
	0.259823,0.277472,0.28478,0.317991,0.328678,0.340006,0.347178,0.364177,0.379776,0.38864,
	0.406828,0.423357,0.449758,0.102692,0.147045,0.182398,0.215204,0.225126,0.25326,0.266164,
	0.282182,0.296113,0.307553,0.324309,0.333559,0.350836,0.35871,0.369991,0.392848,0.403631,
	0.102232,0.142173,0.202042,0.218196,0.250876,0.274353,0.28575,0.325354,0.337215,0.353659,
	0.37059,0.38134,0.390472,0.396564,0.413435,0.424163,0.440423,0.0975372,0.166995,0.181601,
	0.207108,0.242737,0.257606,0.264116,0.288758,0.301149,0.323023,0.336523,0.350898,0.3636,
	0.390325,0.397034,0.406764,0.417337,0.106308,0.158234,0.196467,0.236524,0.258404,0.282326,
	0.290998,0.308428,0.323446,0.32786,0.340163,0.346021,0.356046,0.372913,0.393372,0.41828,
	0.438257,0.121812,0.162338,0.190419,0.212004,0.251146,0.268152,0.285244,0.301014,0.318567,
	0.339444,0.361142,0.373572,0.394346,0.399801,0.408563,0.426097,0.43796,0.145873,0.184225,
	0.200397,0.221424,0.236693,0.283515,0.313352,0.319793,0.33676,0.351621,0.366565,0.391208,
	0.406205,0.420038,0.446805,0.461617,0.467987,0.126951,0.147155,0.177399,0.205679,0.222768,
	0.2457,0.276073,0.292251,0.325521,0.347871,0.359908,0.365989,0.38691,0.393366,0.410677,
	0.422749,0.430299,0.0926255,0.165313,0.183415,0.221306,0.245296,0.279918,0.29288,0.302916,
	0.311796,0.318463,0.3419,0.352441,0.370232,0.378516,0.385896,0.401631,0.409055,0.132088,
	0.154348,0.185739,0.198794,0.213601,0.236907,0.250628,0.263525,0.280848,0.303458,0.320056,
	0.32872,0.337032,0.35947,0.369323,0.390588,0.397284,0.107452,0.178166,0.19498,0.218967,
	0.251285,0.274866,0.287295,0.304418,0.312874,0.328455,0.347097,0.369123,0.383139,0.404648,
	0.418848,0.439186,0.448115,0.160502,0.167729,0.208132,0.235969,0.255649,0.273104,0.282156,
	0.292743,0.306677,0.325346,0.333234,0.340617,0.354114,0.362994,0.374284,0.384236,0.397615,
	0.122517,0.173391,0.212969,0.226017,0.249399,0.257665,0.275866,0.295666,0.301476,0.305935,
	0.324604,0.338282,0.376935,0.397586,0.420077,0.429148,0.439342,0.146627,0.18519,0.215078,
	0.253996,0.269057,0.281222,0.297689,0.307957,0.330116,0.352681,0.367494,0.375439,0.39794,
	0.410285,0.424246,0.435729,0.453332,0.170265,0.194879,0.234531,0.270917,0.293073,0.307564,
	0.314507,0.333272,0.3497,0.352121,0.36727,0.384403,0.394541,0.412919,0.42767,0.443735,
	0.451425,0.141126,0.177858,0.234747,0.255151,0.274855,0.288148,0.300371,0.31379,0.32754,
	0.351761,0.366432,0.386421,0.407136,0.419659,0.431686,0.457918,0.470111,0.149923,0.174644,
	0.197092,0.236323,0.263329,0.281669,0.29105,0.297966,0.329311,0.355383,0.366161,0.385621,
	0.399445,0.412529,0.426503,0.431387,0.440149,0.0800145,0.138678,0.164506,0.234237,0.256929,
	0.265171,0.293889,0.300691,0.326141,0.337138,0.354166,0.362909,0.378904,0.395488,0.410782,
	0.427502,0.441895,0.134481,0.157854,0.177082,0.208822,0.231677,0.247232,0.269131,0.28039,
	0.306232,0.32259,0.338923,0.352059,0.367572,0.379844,0.383661,0.396363,0.401896,0.110722,
	0.134867,0.170162,0.183199,0.228039,0.245007,0.261427,0.27318,0.283768,0.304406,0.310975,
	0.329279,0.360268,0.374565,0.391752,0.411125,0.429402,0.101687,0.165464,0.208786,0.226307,
	0.243687,0.248975,0.286018,0.301468,0.315596,0.324623,0.338518,0.344203,0.360828,0.373918,
	0.396109,0.406549,0.414514,0.126005,0.153,0.205003,0.222491,0.251502,0.267004,0.280232,
	0.296427,0.322846,0.334958,0.34352,0.357828,0.363008,0.379448,0.394273,0.403915,0.434501,
	0.171276,0.21166,0.23394,0.275293,0.298253,0.310992,0.329059,0.340245,0.349497,0.362799,
	0.380313,0.39856,0.409754,0.422129,0.431921,0.440535,0.448178,0.095994,0.167235,0.204207,
	0.221051,0.247439,0.277378,0.291544,0.315528,0.327202,0.354091,0.368442,0.375527,0.386373,
	0.410704,0.423064,0.439788,0.446244,0.17521,0.207647,0.237518,0.245324,0.263182,0.273782,
	0.287366,0.302517,0.310782,0.324088,0.330511,0.357191,0.362749,0.386814,0.413168,0.4274,
	0.435504,0.144333,0.167683,0.233047,0.251958,0.280602,0.292576,0.305382,0.315318,0.332846,
	0.348909,0.361805,0.379305,0.390702,0.394901,0.407814,0.428001,0.452755,0.15149,0.16887,
	0.215958,0.224958,0.240594,0.266098,0.286118,0.295569,0.316947,0.330042,0.343859,0.351846,
	0.359638,0.372771,0.399002,0.427983,0.43527,0.154237,0.186336,0.201139,0.229845,0.267336,
	0.280731,0.286708,0.30371,0.331119,0.355971,0.378128,0.386137,0.397969,0.422344,0.43592,
	0.445377,0.476323,0.126878,0.169464,0.194318,0.215073,0.232171,0.267155,0.277833,0.298943,
	0.318473,0.326046,0.333748,0.345387,0.355628,0.370192,0.380475,0.402919,0.420005,0.105944,
	0.155164,0.183592,0.210032,0.251364,0.260218,0.293509,0.312432,0.319093,0.338167,0.359003,
	0.377474,0.384076,0.40111,0.426871,0.455342,0.468271,0.12853,0.1716,0.212701,0.231142,
	0.242722,0.259707,0.291283,0.309241,0.327542,0.349888,0.367441,0.386394,0.417207,0.426134,
	0.442031,0.472049,0.484847,0.116511,0.160902,0.185523,0.208517,0.227164,0.255335,0.28034,
	0.295692,0.308324,0.316198,0.33412,0.34389,0.364056,0.383989,0.397602,0.407816,0.420838,
	0.14102,0.21655,0.233815,0.256958,0.270011,0.281732,0.294086,0.306958,0.321959,0.330551,
	0.350648,0.364735,0.370392,0.380508,0.393217,0.407787,0.424272,0.111315,0.161844,0.212348,
	0.227761,0.235386,0.249042,0.258482,0.275538,0.295966,0.313748,0.327774,0.353991,0.368047,
	0.376177,0.404187,0.429931,0.43162,0.0921289,0.134095,0.171704,0.199914,0.23029,0.241072,
	0.253931,0.271792,0.300041,0.31155,0.326102,0.333513,0.361519,0.381523,0.398326,0.424972,
	0.435111,0.129205,0.161745,0.17365,0.20758,0.261526,0.27555,0.291034,0.29973,0.324773,
	0.340285,0.363313,0.376933,0.39589,0.40274,0.411361,0.427181,0.440307,0.133302,0.187836,
	0.203969,0.233746,0.254548,0.272062,0.283171,0.30043,0.319237,0.330325,0.35269,0.362463,
	0.372526,0.38829,0.417213,0.447368,0.460206,0.104081,0.149678,0.180593,0.2053,0.215205,
	0.238992,0.251275,0.265524,0.290347,0.295594,0.304997,0.318963,0.336727,0.353796,0.359937,
	0.377334,0.400111,0.108229,0.165908,0.19847,0.222532,0.264097,0.284268,0.295874,0.312794,
	0.331847,0.341866,0.373305,0.387603,0.399041,0.414709,0.434664,0.440666,0.448808,0.119767,
	0.16609,0.188209,0.230584,0.252668,0.267083,0.288471,0.31319,0.325652,0.343862,0.367071,
	0.379026,0.387402,0.39708,0.412145,0.419317,0.431462,0.167878,0.185395,0.211356,0.223693,
	0.237904,0.258839,0.276718,0.292701,0.301774,0.313579,0.332859,0.368396,0.380461,0.393463,
	0.403317,0.410669,0.42356,0.145044,0.170042,0.210483,0.244667,0.261349,0.274108,0.286501,
	0.307949,0.325912,0.357127,0.368403,0.381374,0.397555,0.411255,0.424069,0.442434,0.448239,
	0.0979732,0.172811,0.189148,0.22562,0.25198,0.269296,0.279077,0.297874,0.317108,0.337144,
	0.357638,0.36328,0.378136,0.392775,0.40892,0.415119,0.431807,0.149189,0.173689,0.205869,
	0.223509,0.232954,0.258752,0.271394,0.28447,0.305805,0.320351,0.342637,0.366588,0.388868,
	0.399163,0.404764,0.420848,0.444713,0.117016,0.152456,0.181553,0.228682,0.261749,0.280662,
	0.29234,0.309637,0.321458,0.337131,0.364598,0.378604,0.391999,0.405544,0.415174,0.425551,
	0.441398,0.123356,0.137056,0.178879,0.210323,0.23963,0.256064,0.282427,0.296938,0.317371,
	0.331644,0.335862,0.344934,0.361796,0.373921,0.382376,0.391557,0.404396,0.0887843,0.119549,
	0.171894,0.185293,0.224161,0.236917,0.251656,0.280974,0.287476,0.304284,0.31592,0.3332,
	0.361122,0.377494,0.402478,0.419428,0.428984,0.139629,0.172956,0.203317,0.252863,0.279087,
	0.292083,0.300012,0.325713,0.346144,0.356779,0.363138,0.377757,0.392265,0.399075,0.415369,
	0.425001,0.439237,0.127579,0.174025,0.205196,0.242894,0.269451,0.282529,0.293913,0.308667,
	0.317417,0.328651,0.351556,0.364834,0.383547,0.395859,0.407078,0.42123,0.443129,0.143805,
	0.181698,0.190338,0.20488,0.23774,0.260342,0.281216,0.290577,0.302948,0.330408,0.349241,
	0.356209,0.360445,0.380799,0.393577,0.406774,0.418652,0.151692,0.170018,0.193165,0.22353,
	0.254601,0.276232,0.298962,0.310689,0.321933,0.328378,0.342388,0.367126,0.380799,0.393165,
	0.418912,0.435348,0.447996,0.181263,0.207006,0.218692,0.236329,0.263523,0.270001,0.305276,
	0.31838,0.336841,0.343051,0.358149,0.382484,0.39198,0.405261,0.417049,0.427502,0.433514,
	0.104675,0.172936,0.218625,0.236659,0.255369,0.270988,0.281836,0.290162,0.307205,0.326767,
	0.339015,0.345348,0.364301,0.374655,0.391672,0.426711,0.449357,0.117019,0.137522,0.200169,
	0.221517,0.249068,0.25663,0.284616,0.295914,0.303893,0.313167,0.327407,0.354145,0.374738,
	0.380887,0.387684,0.401905,0.418495,0.0959244,0.175627,0.199141,0.223644,0.245533,0.273396,
	0.296328,0.308205,0.320982,0.329033,0.340035,0.355362,0.36351,0.383424,0.405369,0.419777,
	0.42758,0.127623,0.163129,0.180614,0.200729,0.214685,0.23047,0.261134,0.282415,0.303543,
	0.325651,0.342053,0.351917,0.358872,0.378091,0.389402,0.406642,0.41466,0.12696,0.221378,
	0.246872,0.262523,0.273393,0.303117,0.325572,0.352206,0.364785,0.371523,0.381304,0.388797,
	0.395909,0.412808,0.426143,0.43224,0.451873,0.120921,0.171943,0.203076,0.23002,0.256145,
	0.272795,0.293393,0.297375,0.314097,0.329657,0.34436,0.367647,0.384509,0.397236,0.411581,
	0.421332,0.432051,0.163033,0.20467,0.243339,0.275509,0.28176,0.295901,0.304209,0.321991,
	0.34065,0.359187,0.376925,0.385328,0.40825,0.416648,0.421853,0.438592,0.463028,0.118533,
	0.163438,0.196345,0.223052,0.229582,0.247638,0.266901,0.285173,0.299904,0.321675,0.337885,
	0.354903,0.368348,0.374718,0.383399,0.39845,0.417502,0.159311,0.215943,0.230677,0.25459,
	0.270915,0.289926,0.299207,0.310201,0.323841,0.334935,0.342722,0.363436,0.382905,0.39213,
	0.396179,0.40616,0.417289,0.12541,0.154593,0.173258,0.207287,0.237859,0.255962,0.273738,
	0.288365,0.314998,0.327935,0.334959,0.363727,0.371745,0.396359,0.40878,0.427952,0.446523
	};
	
	
	//for every test point:
	unsigned counter = 0;
	for ( size_t pts = 0 ; pts != test_points.size() ; ++pts )
	{
		//for every distance class:
		for ( size_t dist = 0 ; dist != classes_to_compute_knn.size() ; ++dist )
		{
			//std::cout << (*classes_to_compute_knn[dist])( test_points[pts]) << ",";
			//if ( counter % 10 == 9 )  std::cout << std::endl;
			BOOST_CHECK( fabs( (*classes_to_compute_knn[dist])( test_points[pts]) - result[counter] ) <= 5e-07 );
			++counter;
		}
	}		
}	




BOOST_AUTO_TEST_CASE(Distance_to_k_th_nearest_neighbor_periodic_domain_5d_kd_trees)
{
	//and now kd trees:

	std::vector< std::pair< double , double > > coordinates_of_grid(5);
	coordinates_of_grid[0] = std::pair< double,double >(0,1);
	coordinates_of_grid[1] = std::pair< double,double >(0,1);
	coordinates_of_grid[2] = std::pair< double,double >(0,1);
	coordinates_of_grid[3] = std::pair< double,double >(0,1);
	coordinates_of_grid[4] = std::pair< double,double >(0,1);

	
	//150 random points from R^5.
	std::vector< std::vector<double> > point_cloud = 
	{
	{0.0834491136,0.1779654133,0.4660556072,0.9012534379,0.5278197562},{0.4226206406,0.9520258962,0.758530349,0.8480201229,0.9766116729},{0.5977266638,0.1149142182,0.5232524716,0.1020889988,0.5119932909},
	{0.2711605153,0.8722113823,0.1650474358,0.7590726558,0.6219455686},{0.2334552624,0.1426430494,0.3078001516,0.5127559074,0.7537710918},{0.5213451183,0.048260008,0.1171943273,0.8425379661,0.4401132523},
	{0.4934588228,0.9310606949,0.6055611819,0.8583622978,0.1591549711},{0.6904272269,0.3707491688,0.9402776032,0.2538115475,0.4506233654},{0.6859491069,0.0036608037,0.10070653,0.9493759079,0.467630194},
	{0.8316467148,0.5517980717,0.0530057724,0.541103472,0.6380860347},{0.5978446342,0.0052102418,0.1664767258,0.1295384325,0.5667432889},{0.1400902977,0.3020685925,0.9834268689,0.4287117575,0.7183386472},
	{0.7368947901,0.085055118,0.8047293324,0.3464451593,0.5628888928},{0.7584376174,0.6761060369,0.2481096517,0.4506944369,0.2866351088},{0.0174080082,0.5826915498,0.2013254741,0.0690666109,0.6637215966},
	{0.8963396226,0.4346969973,0.5366920386,0.1864755638,0.793938949},{0.4703953473,0.8998906589,0.2146033144,0.7245412155,0.5346764219},{0.8069159831,0.0554156622,0.435896911,0.5484183766,0.5096576582},
	{0.798520593,0.0641915055,0.3239880204,0.9386689118,0.0557170159},{0.530292036,0.2768214454,0.1465934187,0.2411683465,0.6128252845},{0.3530370982,0.1426775174,0.6023191165,0.561994422,0.4964279565},
	{0.324631419,0.7056863343,0.7003368,0.390044505,0.7481189952},{0.6095938848,0.2573646542,0.7254240075,0.8075247274,0.1642136844},{0.393126196,0.5145787923,0.094276082,0.7294868501,0.9212118047},
	{0.5256770824,0.9961912902,0.6040280494,0.7194164298,0.6646895669},{0.5757870409,0.9228846987,0.5187003482,0.5765369246,0.0548731403},{0.01892799,0.7966554496,0.4025369862,0.9533785735,0.657591725},
	{0.6390500749,0.2819597092,0.0462486721,0.3264356456,0.0678714598},{0.9890820112,0.6545676773,0.9895538145,0.6786556507,0.6729139483},{0.3993956982,0.230341801,0.6619966628,0.0914732311,0.8457622584},
	{0.1578961639,0.395656131,0.7078005325,0.31785637,0.799688499},{0.1970842553,0.2578141042,0.6121606366,0.895383196,0.0099467547},{0.4440435269,0.8098267191,0.3113613564,0.9574479449,0.590720393},
	{0.9892045073,0.6384359652,0.113328327,0.0851192193,0.2706557543},{0.4174217242,0.0269506972,0.3091899455,0.5228551389,0.9364830989},{0.3051545946,0.8718159075,0.8864578567,0.8777936108,0.2914440557},
	{0.5721974459,0.255364249,0.2955540344,0.2383654471,0.4111704903},{0.5746684466,0.8488304173,0.3052570869,0.5827505253,0.0440383067},{0.8940656511,0.6254615316,0.7737742225,0.7862985716,0.5748524205},
	{0.1149175356,0.1231047418,0.0073562327,0.1489178708,0.1661646403},{0.8169913939,0.0589694376,0.6719363823,0.706263863,0.6720062813},{0.1250505799,0.6658676548,0.7596249031,0.2439054565,0.8981688344},
	{0.7740507799,0.9550943908,0.5599841925,0.8664209775,0.2164539213},{0.7468787492,0.3491825978,0.3749921625,0.5969773906,0.5457275284},{0.5082658448,0.3032837098,0.6396805949,0.0767936446,0.6775225908},
	{0.3115240003,0.3442898365,0.4645590561,0.6598785846,0.0875235305},{0.5962873036,0.5788370646,0.866793555,0.8045623263,0.8786653674},{0.2912724817,0.1452073315,0.1192471741,0.3550970335,0.6413987856},
	{0.6897958643,0.7324675771,0.2894716042,0.1806752072,0.5615636073},{0.0443081178,0.4395110959,0.3479123823,0.4973605119,0.0052013861},{0.2848336392,0.6381691515,0.9037866145,0.5798630451,0.8606929414},
	{0.9962989811,0.6965007521,0.7684573706,0.6817337216,0.0348774616},{0.4469122239,0.2762809929,0.7385759256,0.0328260458,0.9326410054},{0.5957374696,0.1654470081,0.3903057303,0.081382965,0.8578045396},
	{0.2540016014,0.6546578028,0.8285822051,0.4711142455,0.5774126945},{0.4869810222,0.3995124626,0.6691552743,0.7868390677,0.8419497332},{0.3616856213,0.6670947331,0.0732981151,0.7414509852,0.5219372443},
	{0.2986048725,0.2752836884,0.9118515207,0.3655172742,0.1339146204},{0.9688093858,0.9839632511,0.0486876701,0.6362382285,0.57918212},{0.6837973881,0.9770498336,0.9009042082,0.9714206185,0.6519823587},
	{0.7795881715,0.0317756881,0.8326280683,0.8814528899,0.8136011944},{0.5725898354,0.2234702245,0.5974976968,0.4569877312,0.213307021},{0.7687702996,0.4392909671,0.8704199663,0.921301692,0.8155898936},
	{0.8652025636,0.4180430397,0.9787512729,0.6994813962,0.6923105665},{0.0511887514,0.4375224516,0.9453371717,0.0047599326,0.0097769545},{0.724004366,0.6219634346,0.1438660068,0.7190580415,0.5052032052},
	{0.4562779004,0.3600603256,0.9501040571,0.7044810567,0.8574787069},{0.8937145416,0.164020116,0.058241178,0.106201407,0.7696755296},{0.4888281797,0.766221792,0.936441195,0.062458213,0.0697754715},
	{0.9425974668,0.8562304049,0.7708064395,0.2568905996,0.6197091083},{0.335584769,0.1298263362,0.4283242184,0.8620519785,0.7777455822},{0.8629977736,0.8917749887,0.3128077728,0.0634608201,0.926209972},
	{0.7719056602,0.861920055,0.0497314511,0.6510132572,0.5403517021},{0.3530974325,0.0193289635,0.221168828,0.3852439316,0.0185912373},{0.0121225882,0.0861171798,0.5727441111,0.1638555448,0.6454807117},
	{0.8268414487,0.6949828749,0.8489062793,0.0978258594,0.0601115541},{0.1035729912,0.560648279,0.9121483676,0.2378044515,0.6134148247},{0.7951907262,0.9540303422,0.166627543,0.1205986922,0.240749388},
	{0.2657136505,0.4355481011,0.1904063555,0.574914173,0.8218553839},{0.2939299098,0.2204667355,0.9821966698,0.679141239,0.5713087774},{0.4602067114,0.4914364009,0.4681642714,0.1970379374,0.7063334736},
	{0.0375891852,0.0569641995,0.7455471558,0.2937743538,0.4793669109},{0.3807092302,0.5626776498,0.0227227262,0.3776322922,0.3102924996},{0.1073450749,0.1354784286,0.6111545078,0.0798917336,0.883374803},
	{0.8246012561,0.9048718265,0.0842031094,0.922818,0.1667780017},{0.3447176681,0.8322002466,0.6765862678,0.1630718487,0.4693403409},{0.1308114366,0.1080321325,0.1855266048,0.1161769815,0.547545319},
	{0.7396875334,0.4003160442,0.9041407737,0.4025102544,0.3810156607},{0.6799464941,0.1193180687,0.6461838551,0.203010628,0.9633627844},{0.6320310263,0.8174625379,0.2596935623,0.5561422193,0.5005034981},
	{0.2930128623,0.9211458443,0.8391498474,0.1270901843,0.6022140398},{0.0262923602,0.2269908616,0.423874124,0.6669442728,0.7403559329},{0.6004114933,0.0432240716,0.6893737733,0.7959725389,0.3219275989},
	{0.8859459318,0.6765051291,0.2261610297,0.9232959126,0.3750975686},{0.291880592,0.6774714992,0.4957603388,0.8096024855,0.9660385114},{0.0824616721,0.606248806,0.1800767379,0.0973097328,0.6039455077},
	{0.9669107667,0.8577861548,0.7781087584,0.9469715115,0.6491487199},{0.7758012626,0.4867858763,0.9445895557,0.1278092652,0.2407444329},{0.1551746733,0.9757911062,0.7507881073,0.1180643593,0.1485485181},
	{0.3832001633,0.3472629497,0.5424197274,0.2347819048,0.5460624669},{0.0356127378,0.3928180397,0.3578011249,0.4599245598,0.5875645839},{0.0915482768,0.8792861651,0.5032183493,0.4352760403,0.4537544101},
	{0.0198585391,0.2713452806,0.2556442935,0.1045154808,0.9693496483},{0.2023309676,0.9850838145,0.0787245242,0.7884351299,0.4749226195},{0.903540008,0.8357364505,0.0291661206,0.6786283345,0.6844293238},
	{0.6532926024,0.9651191607,0.0164092269,0.5536394559,0.20937135},{0.4534949544,0.2433945769,0.7650616835,0.8789803525,0.1802423056},{0.0479354102,0.5243833917,0.1338649723,0.3890133824,0.7519316291},
	{0.9194105798,0.1130673101,0.6573031303,0.8761956429,0.8871361001},{0.4819012247,0.0394615314,0.2906131288,0.5298820648,0.139980295},{0.9022430899,0.9868422102,0.7773927888,0.054481644,0.1240634422},
	{0.7944725945,0.8888423275,0.152707455,0.4338914778,0.9359518671},{0.1997224065,0.1918049152,0.1200931163,0.5627134435,0.752902451},{0.3878603228,0.8341636786,0.6820615984,0.2949337328,0.9880458449},
	{0.7945821437,0.5599427966,0.0107656661,0.3007575846,0.0132312805},{0.4316681421,0.2866483754,0.7283170826,0.9506913973,0.781435441},{0.4010647747,0.8797198473,0.1713188484,0.8730948223,0.4066781153},
	{0.7693045647,0.315029734,0.0817150525,0.7694409494,0.0741932406},{0.7115833587,0.2275419647,0.9498028783,0.1129892636,0.8711559433},{0.1673052812,0.8056030611,0.5349432651,0.2845978811,0.3322721664},
	{0.7902489975,0.2670829475,0.3567338747,0.7903393579,0.8519300122},{0.9817763602,0.5290622578,0.8170695412,0.1410530538,0.174742341},{0.5457802389,0.5899946659,0.0248836742,0.2630379079,0.9055795665},
	{0.6633629678,0.3279535186,0.8295126562,0.308469286,0.6487535806},{0.6586975821,0.0229490069,0.0493111487,0.8900570092,0.2659982985},{0.170515669,0.2176624145,0.3272540492,0.2422923984,0.5008245346},
	{0.5240949828,0.7313395606,0.8073427537,0.143435156,0.0210109602},{0.6617908955,0.9260073218,0.9820222741,0.8356835872,0.5093474893},{0.5819300218,0.2114108466,0.2355436478,0.4126429521,0.8040479794},
	{0.8829782095,0.0962476581,0.0428924561,0.0863333354,0.4965009892},{0.2279743734,0.8029537937,0.1028554635,0.6660963602,0.6046901501},{0.431163084,0.645723691,0.9890625442,0.7445918641,0.948092517},
	{0.0918637037,0.7363437202,0.930700952,0.3582612039,0.4598332951},{0.2691765854,0.1872886773,0.9881781302,0.3337876112,0.5859420574},{0.9778241748,0.9903361972,0.2781383186,0.3943724416,0.2876809081},
	{0.3436939148,0.829814577,0.8227764482,0.3955216904,0.6047144763},{0.8864226528,0.8145404733,0.9444013245,0.4600515633,0.7838328138},{0.711984874,0.5752081149,0.2798226338,0.9386363029,0.9097414555},
	{0.177447621,0.5581884892,0.7912682011,0.6638984308,0.5882288776},{0.5586297519,0.3793644228,0.3251419647,0.4594679768,0.1841192399},{0.6058601916,0.8781331608,0.4892053111,0.4373561433,0.3042695038},
	{0.052442651,0.3735250605,0.2064282887,0.8883306035,0.2979985315},{0.4378396841,0.4180250121,0.9897906918,0.2218166715,0.3110161072},{0.5262017339,0.6203828678,0.8560292989,0.8867694151,0.0002818264},
	{0.6696474531,0.3079569174,0.3395057747,0.2152758248,0.7913365713},{0.2165289272,0.3357743523,0.6355468514,0.2182823461,0.9831528417},{0.3455453599,0.4667630245,0.9083769515,0.7256064662,0.292692615},
	{0.2951296815,0.8800065387,0.8613584242,0.5988326801,0.2341764281},{0.0876358571,0.6400671883,0.0284326931,0.5182432309,0.0632727058},{0.1320868398,0.512312592,0.4902214163,0.6152805286,0.3113408913}
    };	
	
	std::vector< Gudhi::Topological_inference_with_cubical_complexes::Distance_to_k_th_closest_point_periodic_k_d_tree* >
	classes_to_compute_knn;	
	for ( size_t i = 5 ; i != 90 ; i=i+5 )
	{
		classes_to_compute_knn.push_back(
		new Gudhi::Topological_inference_with_cubical_complexes::Distance_to_k_th_closest_point_periodic_k_d_tree
		( point_cloud, coordinates_of_grid , i ) );
	}  
	
	
	
	//150 random points to test the distance:
	std::vector< std::vector<double> > test_points =
	{
        {0.6156253421,0.0866876692,0.6546369994,0.7980181805,0.3036833229},{0.3805250507,0.5532490518,0.7087723003,0.8031594998,0.4108131663},{0.3181996481,0.8681068709,0.5779479195,0.1869301829,0.7758894945},
		{0.7863677139,0.3242641296,0.4091913456,0.5515330643,0.8118765294},{0.9259262064,0.4886787878,0.5708910702,0.64194575,0.8815862499},{0.7040850376,0.7995467752,0.5784777766,0.8222728742,0.62802323},
		{0.3415478489,0.1666721171,0.2690072777,0.1288172284,0.7955359623},{0.5544780709,0.6510210049,0.7063087481,0.5599074461,0.2960705815},{0.5480834248,0.7081139782,0.119110408,0.2898809859,0.5940711163},
		{0.403030508,0.4703414962,0.9957842538,0.7432314595,0.2526669232},{0.3263766619,0.6720029712,0.0591752841,0.3556350819,0.4690573232},{0.2064622687,0.5762362943,0.6930124161,0.8952994619,0.0152675807},
		{0.1605475396,0.0204712148,0.1160777733,0.247559764,0.8003176553},{0.6224387698,0.455551483,0.8776516232,0.3776043199,0.2195397872},{0.9540606579,0.1067592113,0.2777793263,0.1421314587,0.5134446463},
		{0.0406719297,0.6811534585,0.8866318022,0.4287882582,0.5006276881},{0.870287016,0.2856612143,0.1356323606,0.8231363078,0.335178002},{0.7992206207,0.2248386785,0.652452287,0.4542188642,0.1205765211},
		{0.5416745862,0.4066060479,0.8934159388,0.9497496316,0.7745553942},{0.5069284546,0.7684237019,0.5580311327,0.9252927683,0.7101417952},{0.0541560685,0.3102081812,0.5665760315,0.4848779158,0.047492543},
		{0.7745181869,0.5853341455,0.5212711582,0.7912951233,0.8805664999},{0.2830229441,0.032924467,0.8854985246,0.2631558124,0.6652764447},{0.6711386403,0.7294949514,0.078212844,0.8218062653,0.4102030019},
		{0.1846682311,0.2522296603,0.0682438477,0.3138512317,0.5468119136},{0.1068816297,0.8906914974,0.3227836362,0.9410896939,0.7084375294},{0.6661021973,0.2722590868,0.9509407876,0.3421589313,0.0829368159},
		{0.2404863781,0.9042210197,0.214130345,0.7180377392,0.5869330307},{0.4962969481,0.7940989283,0.9015325028,0.8254819594,0.1567669387},{0.7736787216,0.1744941329,0.021256563,0.6872016888,0.2290441715},
		{0.176827118,0.524485471,0.7884605827,0.9387880943,0.8023405797},{0.9890728234,0.6826795805,0.712908695,0.592580826,0.7100058955},{0.8278585384,0.3474175616,0.7225242155,0.4581460671,0.041690419},
		{0.3977338392,0.7302346188,0.6874800553,0.5466448646,0.0409471621},{0.2029349247,0.9303811735,0.3612236425,0.4019850001,0.7212867679},{0.0979896556,0.5700703752,0.8520184332,0.7930803259,0.1817028737},
		{0.4875848652,0.2855914799,0.1294624861,0.6426335031,0.2070947611},{0.9098975097,0.2684065762,0.4937620726,0.2055649131,0.5361320146},{0.4482507212,0.4339038245,0.4895149507,0.9056274116,0.6298295192},
		{0.3551141147,0.7516602229,0.0647157019,0.6934965344,0.9270651848},{0.9786594613,0.4435110099,0.3609527228,0.8899996865,0.3251161473},{0.153831613,0.4260658224,0.334146712,0.7129775598,0.853001711},
		{0.9990624285,0.3813134674,0.3835411465,0.098076452,0.4698095596},{0.9865754729,0.0481281504,0.9224064483,0.7275167136,0.2618650477},{0.0151539787,0.7301042816,0.1246206658,0.0785468051,0.6520250721},
		{0.1044145077,0.2527292348,0.3791852833,0.8313468685,0.6523172415},{0.677102146,0.5618522577,0.0822143261,0.3621671796,0.9020784714},{0.9044152494,0.487927774,0.4113524007,0.5217605885,0.4028691752},
		{0.9884912123,0.1933964493,0.8975249475,0.9582569813,0.0643633602},{0.6265765063,0.6546486341,0.6538947937,0.7916958395,0.2581015776},{0.1925964975,0.8758495869,0.7928991041,0.8404584907,0.0501623263},
		{0.5288344054,0.0940175443,0.3880097063,0.5240216907,0.2688257839},{0.151839297,0.4815489671,0.1200810564,0.5240271261,0.4553921684},{0.6643557919,0.7364383612,0.6186252583,0.9226146308,0.033159988},
		{0.7722630778,0.9259055518,0.8249877302,0.7974182731,0.7244633394},{0.1482382049,0.1489458592,0.760812395,0.4457491287,0.1411920183},{0.8833106682,0.1271617629,0.7506208313,0.9937048177,0.6780842175},
		{0.0076348747,0.578255268,0.6369468854,0.1451479641,0.5575621051},{0.4688744638,0.067973051,0.9584744396,0.4920711045,0.5704484717},{0.7850542171,0.0124715585,0.3579944114,0.5550846129,0.9165706104},
		{0.783577254,0.4116507513,0.7414916439,0.8980796521,0.4050635439},{0.1215615212,0.7875910646,0.141919357,0.5981509653,0.5387970002},{0.7539457558,0.4804520637,0.4587839125,0.2325882488,0.6723778283},
		{0.5750644468,0.3006371737,0.1820672811,0.0022302219,0.9543300846},{0.1849272884,0.3658491941,0.5859630774,0.2220380711,0.1284696721},{0.6724957384,0.4000277666,0.9780439236,0.0553804729,0.3112557554},
		{0.9480736756,0.5819517735,0.4126889601,0.2862954626,0.9022503761},{0.3177717375,0.2309198242,0.1540829227,0.7880887291,0.8507864},{0.7116010217,0.6213103663,0.2930113161,0.3727468587,0.9516895316},
		{0.3819570281,0.8180342214,0.2295541067,0.7412023214,0.8717434809},{0.7886661366,0.2495642588,0.048790921,0.1196943321,0.3640789606},{0.4985892794,0.0430643624,0.1038529542,0.9094649146,0.3767474603},
		{0.4307315028,0.4978282363,0.0282152942,0.2946561221,0.888703462},{0.6723954086,0.784596333,0.0863898948,0.7806103502,0.2525719581},{0.5714831352,0.4371523738,0.4938304564,0.8903684383,0.567739527},
		{0.3035034607,0.9636900444,0.9626759095,0.9843311792,0.2204016687},{0.4772450773,0.2324201055,0.7633715817,0.3358677777,0.9510003324},{0.3538140105,0.9728750708,0.5895036547,0.9765105813,0.6960510097},
		{0.857701164,0.8197515106,0.1978863701,0.0605267917,0.5568545174},{0.1860871471,0.5490134645,0.7524937373,0.3706550032,0.4462244327},{0.4302063505,0.2829570582,0.6203799653,0.6574771008,0.9612034687},
		{0.8402685211,0.9488183272,0.3951937449,0.683283203,0.6873453115},{0.5700790631,0.1601365875,0.7559559469,0.1201780995,0.2192969187},{0.5646118473,0.1211143576,0.7025673755,0.6807703031,0.0820197742},
		{0.7155782501,0.3898127335,0.8637648029,0.2135735601,0.130055618},{0.5476660505,0.0127128572,0.0949698025,0.1445218846,0.094173088},{0.7984791892,0.6017809003,0.4891143411,0.4428788116,0.5716396598},
		{0.2536687476,0.3673323095,0.6265329041,0.6520072613,0.8490074677},{0.7154158063,0.3787585783,0.1292777357,0.3725556275,0.3312357822},{0.4494496752,0.0517784276,0.2138768274,0.3533489655,0.6084807601},
		{0.2920235356,0.1409710939,0.1904358806,0.4535624995,0.9371193126},{0.5526888377,0.5882126573,0.827937352,0.858154837,0.7465898003},{0.5247145414,0.6173376355,0.8801819091,0.2626326918,0.0044825338},
		{0.3793205258,0.2487422219,0.5722060839,0.229396201,0.0556458414},{0.3603424882,0.3157115481,0.4432009482,0.318166235,0.0832488397},{0.1691172069,0.4944495703,0.1648998465,0.7261980609,0.8623710414},
		{0.3742564206,0.7017267803,0.7860315635,0.5872677292,0.7537247543},{0.7428908623,0.0345434777,0.9549517429,0.8074096772,0.849237503},{0.7985906939,0.0946825973,0.3014169701,0.8058131211,0.1560831016},
		{0.2603009364,0.1148938027,0.8615724845,0.4697891294,0.9368353917},{0.5338982816,0.2307557308,0.6977123844,0.2302658497,0.1019402801},{0.9767520463,0.718872142,0.6688464943,0.4160435961,0.9022799341},
		{0.0166666615,0.4389940493,0.7579744349,0.5370429766,0.1495361354},{0.3979732869,0.4755268127,0.5795001904,0.9908350941,0.1669989298},{0.8949154501,0.2682687736,0.6063345135,0.2275951863,0.017111354},
		{0.4043473324,0.3614203415,0.5158000104,0.1085771453,0.8999612131},{0.977441187,0.1501589406,0.8518329735,0.3219707806,0.8540179215},{0.4475535215,0.0862195999,0.8618425971,0.9331469764,0.4519766094},
		{0.4441091421,0.9066383392,0.0553346425,0.5744697268,0.9810071255},{0.9267698682,0.9928793737,0.5727517998,0.1823302133,0.5895151296},{0.7495447183,0.6460728368,0.470270294,0.6661383023,0.0136451058},
		{0.6532815052,0.325947958,0.675558131,0.2983020572,0.5749601708},{0.59615145,0.4472527723,0.581177111,0.5676356531,0.77461539},{0.7588686878,0.4715215371,0.7138359281,0.6454586445,0.6092199853},
		{0.6003141189,0.7735408389,0.7334262682,0.5612569719,0.5469313441},{0.8622437229,0.4840239387,0.3084369798,0.2880281513,0.3305151053},{0.7834043312,0.2515344019,0.8712158261,0.3956733951,0.6839187194},
		{0.7728605138,0.4649363579,0.1062444423,0.2089061183,0.2042179685},{0.2607685262,0.5400420581,0.7663626072,0.5831070796,0.1841536507},{0.4371812299,0.9649007791,0.4027777456,0.0520443097,0.6338951529},
		{0.5583380917,0.7672120824,0.3233393638,0.3555802056,0.5658726953},{0.5354220215,0.7534636182,0.9579864799,0.6235753896,0.9860881409},{0.9109844903,0.9209281902,0.2032730938,0.0754138096,0.4649585146},
		{0.20800173,0.8517249846,0.9164571092,0.2328442468,0.1823428017},{0.2348916477,0.4451530215,0.5246549759,0.7947810211,0.6921255561},{0.6184620375,0.1530659061,0.9847622805,0.7646336413,0.6839535253},
		{0.4769351536,0.1649140997,0.5587477549,0.4871179559,0.1689623857},{0.0984163808,0.6076200083,0.966948214,0.1425744521,0.2599090473},{0.6092678127,0.14414063,0.6357082536,0.4362932767,0.9267203519},
		{0.2324069291,0.3263996779,0.682590036,0.591015572,0.0428221659},{0.2328919342,0.7355384333,0.8440155704,0.7060465526,0.3706862081},{0.8768993784,0.8748575826,0.4019506322,0.9472513744,0.7443623699},
		{0.0835958309,0.2113429978,0.4997533069,0.3181976748,0.5457607873},{0.4077742845,0.1392330623,0.8563316758,0.1835489655,0.016605665},{0.9483211499,0.9293324859,0.939470334,0.891758071,0.5508040495},
		{0.1962504159,0.1287407137,0.7581467074,0.7718736557,0.280741273},{0.4436504522,0.8687367935,0.5447200509,0.848678587,0.4843223449},{0.5696217748,0.9954605456,0.6056541163,0.2167235308,0.562121178},
		{0.9663550276,0.1291832027,0.860948727,0.2084195951,0.1163226201},{0.1134767649,0.3865185466,0.9239129042,0.7176613885,0.3527642083},{0.3190164822,0.7672213989,0.2879216047,0.6133680425,0.4904669598},
		{0.2613077373,0.1525334122,0.048044354,0.678694438,0.5231821293},{0.1799443669,0.7424147304,0.896189518,0.6574717208,0.9911259802},{0.043235773,0.0527958591,0.0389118574,0.0312990851,0.419792705},
		{0.1913606562,0.6254830866,0.5446361115,0.6739698576,0.0345029861},{0.9619858956,0.8867313419,0.8675114349,0.7746892928,0.0028468464},{0.4198780232,0.6695927838,0.4730970562,0.8541963291,0.0656680651},
		{0.4705536903,0.2660579658,0.9417342963,0.9129116719,0.1820742511},{0.2058038623,0.0234297065,0.9110244056,0.4768905714,0.1417649831},{0.2174516551,0.0631806597,0.2748815217,0.5283225405,0.3984943791}
	};	
	
	
	std::vector<double> result = 
	{
	0.0848295,0.167682,0.190143,0.224555,0.260605,0.276135,0.299343,0.308722,0.325974,0.332576,
	0.355922,0.363254,0.382531,0.395273,0.411649,0.428977,0.445709,0.162351,0.213094,0.233111,
	0.273276,0.291066,0.299465,0.307219,0.316869,0.329778,0.349699,0.363402,0.369242,0.379291,
	0.393963,0.399038,0.416935,0.437603,0.129412,0.175697,0.205706,0.222322,0.237218,0.257294,
	0.281272,0.290291,0.30679,0.32004,0.334267,0.356391,0.375863,0.399707,0.420824,0.42818,
	0.437744,0.123879,0.183859,0.223911,0.250866,0.272731,0.284377,0.29067,0.316663,0.330722,
	0.346474,0.360469,0.373744,0.384353,0.390481,0.416428,0.424699,0.435283,0.17581,0.203497,
	0.220094,0.234968,0.271108,0.301817,0.316558,0.327927,0.34733,0.360667,0.371183,0.385407,
	0.396205,0.410337,0.442214,0.456112,0.460307,0.148153,0.185379,0.206594,0.246517,0.264187,
	0.281161,0.293867,0.312138,0.321497,0.331995,0.347012,0.355599,0.36766,0.375588,0.405274,
	0.415481,0.440765,0.116477,0.145421,0.171826,0.210955,0.23101,0.248492,0.271623,0.302032,
	0.308633,0.314454,0.327844,0.339908,0.353454,0.367661,0.385579,0.404055,0.420144,0.168283,
	0.212246,0.237069,0.250381,0.264192,0.266797,0.2842,0.291756,0.305796,0.324195,0.346134,
	0.365502,0.377756,0.394817,0.413001,0.428632,0.442215,0.146687,0.187936,0.211868,0.235742,
	0.247644,0.265255,0.272803,0.282406,0.296366,0.310016,0.326382,0.334602,0.347457,0.370799,
	0.38372,0.394413,0.399228,0.130958,0.180165,0.212301,0.226185,0.251889,0.278232,0.292023,
	0.31947,0.334075,0.342909,0.362938,0.37578,0.377176,0.386597,0.407513,0.414678,0.438542,
	0.118381,0.18156,0.209455,0.233214,0.244773,0.256898,0.269794,0.29571,0.301167,0.306994,
	0.322214,0.330463,0.341341,0.35721,0.385369,0.406063,0.417629,0.131054,0.17019,0.193571,
	0.209064,0.233001,0.26,0.275177,0.305442,0.315445,0.343536,0.362786,0.381563,0.398439,
	0.428502,0.440018,0.444011,0.461674,0.116072,0.15786,0.197647,0.238552,0.25556,0.265478,
	0.28922,0.296169,0.30705,0.316529,0.32188,0.333712,0.347416,0.36699,0.375887,0.383202,
	0.40527,0.0918235,0.157325,0.209278,0.256533,0.279914,0.293528,0.31449,0.318986,0.34049,
	0.367845,0.375007,0.388779,0.39476,0.406577,0.418407,0.425026,0.43391,0.115489,0.153122,
	0.20137,0.214372,0.242651,0.259452,0.271232,0.283506,0.308844,0.316242,0.32176,0.346473,
	0.357072,0.371264,0.387054,0.415065,0.425979,0.105873,0.15049,0.178694,0.19109,0.214227,
	0.224868,0.254229,0.269422,0.297074,0.31641,0.340962,0.350929,0.357503,0.379616,0.406447,
	0.441294,0.460655,0.148205,0.184997,0.209429,0.230873,0.252504,0.276562,0.290241,0.299428,
	0.306483,0.318678,0.339948,0.359892,0.376597,0.394306,0.411601,0.429351,0.442051,0.178267,
	0.23121,0.248141,0.261049,0.274753,0.287126,0.301017,0.316553,0.331258,0.346578,0.356343,
	0.371992,0.391065,0.39716,0.405855,0.42154,0.433167,0.0818327,0.130023,0.181482,0.218445,
	0.227169,0.260478,0.276303,0.290187,0.302683,0.314346,0.317565,0.327318,0.350202,0.364067,
	0.378657,0.399229,0.41126,0.158959,0.194153,0.200482,0.239815,0.256202,0.270617,0.279928,
	0.291677,0.324638,0.335548,0.344593,0.353398,0.374008,0.389338,0.407072,0.429214,0.439384,
	0.13937,0.214241,0.249153,0.262912,0.300738,0.319474,0.327932,0.339264,0.345877,0.360758,
	0.377875,0.389604,0.395612,0.412283,0.431313,0.437804,0.459013,0.158461,0.201442,0.229945,
	0.256686,0.2872,0.299562,0.309594,0.321666,0.335328,0.340489,0.347109,0.369803,0.381844,
	0.400009,0.407782,0.423448,0.433999,0.115902,0.164519,0.192033,0.206392,0.23014,0.241564,
	0.25929,0.269414,0.289775,0.299306,0.314701,0.338119,0.347042,0.37126,0.380735,0.38967,
	0.410681,0.0954641,0.120935,0.171904,0.192329,0.208385,0.236661,0.264469,0.290409,0.305562,
	0.32002,0.329735,0.341951,0.362225,0.378627,0.399776,0.412042,0.426937,0.0765244,0.148824,
	0.170286,0.1989,0.225955,0.271435,0.282216,0.299209,0.314327,0.332349,0.349185,0.360365,
	0.372509,0.375316,0.381403,0.394748,0.39824,0.13168,0.155418,0.200041,0.221016,0.244805,
	0.26525,0.291607,0.315401,0.323612,0.342006,0.360835,0.368547,0.37463,0.388037,0.406282,
	0.41778,0.424659,0.128945,0.158617,0.207504,0.232391,0.241996,0.275163,0.291807,0.304589,
	0.328356,0.348878,0.354584,0.361452,0.366657,0.376261,0.388253,0.401706,0.42167,0.0847531,
	0.144881,0.163493,0.204349,0.24786,0.2631,0.274815,0.285477,0.305458,0.322563,0.328079,
	0.351235,0.368374,0.382791,0.397444,0.416454,0.430871,0.0840074,0.135321,0.162997,0.205367,
	0.241821,0.258843,0.272941,0.299455,0.311363,0.333512,0.347029,0.363076,0.379311,0.385606,
	0.402224,0.416354,0.42455,0.139989,0.177858,0.200324,0.249376,0.266311,0.283804,0.293187,
	0.299576,0.321155,0.336324,0.344766,0.356998,0.369519,0.379655,0.396388,0.404356,0.416574,
	0.147076,0.165801,0.177249,0.18273,0.214025,0.235448,0.259896,0.269111,0.296584,0.311088,
	0.317151,0.33647,0.36497,0.372285,0.394458,0.414546,0.4296,0.116688,0.156473,0.178371,
	0.209394,0.243135,0.258517,0.285988,0.311073,0.322853,0.34883,0.364372,0.385156,0.396268,
	0.403747,0.415583,0.43554,0.453863,0.162764,0.198535,0.227761,0.257608,0.273038,0.284108,
	0.29944,0.314218,0.337691,0.346333,0.356919,0.373788,0.384479,0.393329,0.408232,0.415878,
	0.443275,0.116384,0.174319,0.193311,0.221845,0.24301,0.265486,0.280428,0.3086,0.332731,
	0.353183,0.385539,0.400152,0.416846,0.430188,0.446937,0.459844,0.484937,0.138746,0.179503,
	0.21213,0.241658,0.249752,0.261302,0.276726,0.286747,0.295024,0.310876,0.322534,0.3324,
	0.349908,0.356612,0.382253,0.400623,0.419267,0.137542,0.188908,0.204822,0.23484,0.257798,
	0.289237,0.300299,0.329093,0.340438,0.346706,0.352856,0.379992,0.38789,0.409263,0.425078,
	0.437179,0.444264,0.149244,0.164752,0.192312,0.227408,0.248319,0.267239,0.28163,0.294676,
	0.316238,0.329648,0.339191,0.349276,0.369222,0.391152,0.408676,0.420909,0.426995,0.131742,
	0.177568,0.20558,0.233554,0.258632,0.271847,0.281661,0.296821,0.313477,0.332356,0.34915,
	0.352513,0.369843,0.376887,0.389902,0.404039,0.439821,0.129897,0.186526,0.21181,0.243386,
	0.267813,0.28169,0.296669,0.312774,0.326523,0.3337,0.350471,0.366826,0.385017,0.399471,
	0.419172,0.423575,0.438258,0.129045,0.148848,0.181683,0.213869,0.240414,0.259456,0.274177,
	0.292159,0.314131,0.325991,0.33858,0.357101,0.369324,0.376263,0.389903,0.418074,0.430231,
	0.140462,0.215015,0.244659,0.269943,0.298829,0.312447,0.336464,0.347849,0.358189,0.370711,
	0.386266,0.405851,0.418113,0.428324,0.437489,0.446205,0.459274,0.12756,0.157548,0.223411,
	0.235681,0.256191,0.271037,0.276026,0.298884,0.320632,0.349332,0.360572,0.379637,0.392041,
	0.408288,0.41877,0.433888,0.443498,0.116942,0.164221,0.229464,0.25412,0.285386,0.298987,
	0.316394,0.327737,0.341213,0.357346,0.370233,0.387133,0.398025,0.412334,0.421238,0.435437,
	0.441612,0.150679,0.180801,0.211745,0.228479,0.255588,0.272247,0.291648,0.304284,0.310504,
	0.328107,0.337472,0.36513,0.383238,0.389674,0.404231,0.41587,0.426158,0.13067,0.158706,
	0.184978,0.219099,0.234467,0.252476,0.270417,0.27962,0.296749,0.320483,0.329018,0.341914,
	0.356072,0.378279,0.384539,0.390046,0.417365,0.145659,0.182728,0.197889,0.220565,0.231055,
	0.247389,0.27591,0.300598,0.32187,0.329434,0.340978,0.347406,0.356534,0.370672,0.388466,
	0.420762,0.442789,0.131959,0.185812,0.202912,0.218312,0.226104,0.249487,0.270124,0.286387,
	0.297756,0.313003,0.316846,0.334897,0.344937,0.362273,0.399007,0.415054,0.423063,0.171461,
	0.20929,0.238393,0.258245,0.29075,0.310392,0.326796,0.35551,0.360472,0.3673,0.377338,
	0.383488,0.398097,0.417919,0.433751,0.447154,0.46545,0.129293,0.144404,0.17273,0.215114,
	0.248272,0.277652,0.299891,0.30916,0.319353,0.335283,0.358673,0.375009,0.398405,0.412348,
	0.418744,0.426281,0.442124,0.172305,0.205185,0.220035,0.247603,0.277901,0.296177,0.302768,
	0.312524,0.325525,0.355209,0.368689,0.382312,0.392242,0.412473,0.417307,0.428925,0.44487,
	0.107461,0.170044,0.192852,0.226081,0.262557,0.276851,0.2905,0.300934,0.324338,0.345697,
	0.353681,0.377779,0.390725,0.409954,0.420187,0.422235,0.441866,0.0976053,0.139708,0.175301,
	0.228986,0.251131,0.270622,0.292676,0.303091,0.326519,0.343995,0.356215,0.369919,0.389538,
	0.398167,0.420586,0.425948,0.443788,0.13188,0.145626,0.158717,0.191441,0.228119,0.245914,
	0.271726,0.285823,0.314104,0.319178,0.331914,0.3475,0.36598,0.386244,0.403841,0.421966,
	0.435905,0.112575,0.170152,0.1962,0.230681,0.250902,0.279758,0.296844,0.322884,0.345273,
	0.355072,0.36449,0.377843,0.398354,0.410278,0.428198,0.439031,0.453815,0.0827667,0.130232,
	0.178571,0.224195,0.250608,0.268633,0.284809,0.288505,0.306969,0.318577,0.333351,0.340429,
	0.353532,0.373966,0.387767,0.398494,0.412389,0.151294,0.208907,0.24045,0.269328,0.291335,
	0.310029,0.317076,0.327049,0.339157,0.357268,0.361326,0.375106,0.383193,0.392716,0.400556,
	0.413722,0.423666,0.0861021,0.130745,0.163834,0.194855,0.241684,0.251902,0.256774,0.282805,
	0.307534,0.313932,0.335699,0.351102,0.363955,0.377639,0.387762,0.409197,0.424532,0.149395,
	0.173346,0.210495,0.24222,0.255553,0.270058,0.288912,0.323015,0.348704,0.362597,0.379561,
	0.389629,0.399567,0.4079,0.41985,0.432498,0.441143,0.117028,0.161888,0.180571,0.194028,
	0.212216,0.244448,0.260669,0.280721,0.294021,0.297394,0.312806,0.325237,0.347229,0.350885,
	0.357319,0.370377,0.387034,0.128781,0.174014,0.222018,0.232126,0.247162,0.264472,0.285935,
	0.324743,0.338422,0.349808,0.367192,0.375781,0.385897,0.399466,0.414765,0.427602,0.438449,
	0.178487,0.204078,0.25002,0.265134,0.28147,0.286906,0.296593,0.311568,0.319806,0.345901,
	0.354039,0.373016,0.383799,0.393364,0.398804,0.410475,0.427772,0.0898087,0.151021,0.174144,
	0.201235,0.218836,0.232794,0.25022,0.26499,0.293729,0.302069,0.327005,0.337045,0.364822,
	0.386036,0.397253,0.416544,0.444209,0.148769,0.178147,0.201566,0.251541,0.277521,0.288767,
	0.298912,0.308631,0.323918,0.337714,0.358439,0.37182,0.382773,0.388008,0.42289,0.434889,
	0.451944,0.116565,0.169483,0.198001,0.222041,0.24346,0.264877,0.280762,0.293549,0.304554,
	0.319997,0.346035,0.350225,0.364877,0.386287,0.397278,0.400814,0.412414,0.140001,0.177271,
	0.2161,0.239569,0.270861,0.301758,0.320981,0.329298,0.340251,0.370782,0.374977,0.386382,
	0.397013,0.410634,0.417954,0.439699,0.465427,0.152417,0.175266,0.190577,0.20998,0.227925,
	0.258756,0.283242,0.287372,0.308132,0.318032,0.331809,0.347206,0.365245,0.395184,0.418917,
	0.429027,0.446135,0.160541,0.177337,0.208286,0.261503,0.278546,0.307053,0.320308,0.340532,
	0.351113,0.358767,0.371271,0.390767,0.399866,0.407185,0.415349,0.429543,0.448921,0.0981115,
	0.184175,0.220715,0.232539,0.250641,0.265616,0.282201,0.302805,0.30902,0.31595,0.328375,
	0.335658,0.348514,0.354119,0.368222,0.380623,0.390701,0.125477,0.188011,0.201934,0.244118,
	0.280797,0.294112,0.301579,0.309908,0.320634,0.34168,0.364519,0.383665,0.399118,0.403904,
	0.411475,0.425422,0.439594,0.112308,0.162315,0.191481,0.235512,0.242443,0.263878,0.291495,
	0.307304,0.324999,0.32957,0.34535,0.35428,0.384397,0.394822,0.405313,0.416844,0.423492,
	0.116479,0.153696,0.180251,0.19926,0.21744,0.252754,0.261896,0.290729,0.312526,0.333502,
	0.342694,0.35152,0.364659,0.382995,0.395576,0.406227,0.422965,0.0782063,0.148522,0.17906,
	0.199923,0.223234,0.242778,0.25918,0.294968,0.310624,0.323587,0.344621,0.360502,0.370311,
	0.386531,0.399036,0.406945,0.420366,0.145674,0.170536,0.19136,0.196193,0.21689,0.232604,
	0.245754,0.27559,0.287335,0.309016,0.315254,0.34633,0.358415,0.371508,0.381914,0.400031,
	0.411754,0.100028,0.13229,0.165956,0.193402,0.225949,0.244173,0.26185,0.276087,0.286655,
	0.299438,0.311754,0.339737,0.35866,0.383074,0.396801,0.410824,0.449413,0.146488,0.199293,
	0.217977,0.24622,0.263588,0.279784,0.314182,0.326571,0.346735,0.353918,0.3648,0.381911,
	0.401844,0.407745,0.421561,0.436449,0.455741,0.107174,0.152495,0.194754,0.23575,0.253101,
	0.262348,0.282684,0.300142,0.310859,0.32947,0.338408,0.355884,0.369004,0.378954,0.385853,
	0.397568,0.41235,0.109854,0.154938,0.208312,0.233663,0.252158,0.262422,0.280651,0.289883,
	0.301815,0.313918,0.32318,0.342187,0.350699,0.361938,0.380097,0.396628,0.418188,0.113658,
	0.133705,0.180173,0.21138,0.248999,0.258898,0.271047,0.280221,0.297073,0.320469,0.340618,
	0.351849,0.358379,0.386239,0.396668,0.40315,0.412565,0.093209,0.139833,0.179182,0.199063,
	0.23301,0.24386,0.269891,0.277405,0.302523,0.310307,0.333773,0.34459,0.354937,0.372694,
	0.394009,0.405079,0.424012,0.129614,0.152117,0.184388,0.214176,0.249647,0.25893,0.283146,
	0.29842,0.319103,0.339998,0.359316,0.382936,0.392815,0.409367,0.426785,0.436467,0.462834,
	0.120082,0.155997,0.201851,0.217403,0.235926,0.251093,0.275238,0.288446,0.307993,0.319555,
	0.339172,0.359364,0.378091,0.402328,0.415955,0.434145,0.449492,0.130052,0.154217,0.196823,
	0.225966,0.248047,0.264181,0.276009,0.284017,0.297434,0.31933,0.346547,0.353876,0.375991,
	0.391038,0.415219,0.42397,0.433356,0.134664,0.157415,0.183242,0.205096,0.228939,0.247606,
	0.262267,0.271826,0.287477,0.311757,0.325264,0.337299,0.350165,0.389189,0.398486,0.41303,
	0.437195,0.0880821,0.163076,0.185834,0.21319,0.237257,0.253184,0.268479,0.294161,0.32533,
	0.334532,0.36589,0.384427,0.395402,0.407877,0.426187,0.447139,0.457298,0.101,0.1477,
	0.185195,0.21434,0.248383,0.261938,0.300237,0.309411,0.32948,0.343705,0.360918,0.375282,
	0.38906,0.405813,0.415136,0.428248,0.451173,0.117487,0.172719,0.200242,0.21634,0.23322,
	0.245205,0.261798,0.282962,0.30402,0.318726,0.332854,0.34158,0.363101,0.372751,0.385748,
	0.399648,0.408507,0.146499,0.205623,0.223668,0.262026,0.284331,0.305608,0.318508,0.328317,
	0.337227,0.351808,0.36087,0.387653,0.398202,0.412079,0.429956,0.445785,0.461238,0.130668,
	0.15328,0.201685,0.225663,0.234401,0.255873,0.269483,0.290109,0.300037,0.316744,0.344256,
	0.372011,0.388302,0.416099,0.424778,0.453009,0.473605,0.0935834,0.157673,0.209764,0.235338,
	0.259823,0.277472,0.28478,0.317991,0.328678,0.340006,0.347178,0.364177,0.379776,0.38864,
	0.406828,0.423357,0.449758,0.102692,0.147045,0.182398,0.215204,0.225126,0.25326,0.266164,
	0.282182,0.296113,0.307553,0.324309,0.333559,0.350836,0.35871,0.369991,0.392848,0.403631,
	0.102232,0.142173,0.202042,0.218196,0.250876,0.274353,0.28575,0.325354,0.337215,0.353659,
	0.37059,0.38134,0.390472,0.396564,0.413435,0.424163,0.440423,0.0975372,0.166995,0.181601,
	0.207108,0.242737,0.257606,0.264116,0.288758,0.301149,0.323023,0.336523,0.350898,0.3636,
	0.390325,0.397034,0.406764,0.417337,0.106308,0.158234,0.196467,0.236524,0.258404,0.282326,
	0.290998,0.308428,0.323446,0.32786,0.340163,0.346021,0.356046,0.372913,0.393372,0.41828,
	0.438257,0.121812,0.162338,0.190419,0.212004,0.251146,0.268152,0.285244,0.301014,0.318567,
	0.339444,0.361142,0.373572,0.394346,0.399801,0.408563,0.426097,0.43796,0.145873,0.184225,
	0.200397,0.221424,0.236693,0.283515,0.313352,0.319793,0.33676,0.351621,0.366565,0.391208,
	0.406205,0.420038,0.446805,0.461617,0.467987,0.126951,0.147155,0.177399,0.205679,0.222768,
	0.2457,0.276073,0.292251,0.325521,0.347871,0.359908,0.365989,0.38691,0.393366,0.410677,
	0.422749,0.430299,0.0926255,0.165313,0.183415,0.221306,0.245296,0.279918,0.29288,0.302916,
	0.311796,0.318463,0.3419,0.352441,0.370232,0.378516,0.385896,0.401631,0.409055,0.132088,
	0.154348,0.185739,0.198794,0.213601,0.236907,0.250628,0.263525,0.280848,0.303458,0.320056,
	0.32872,0.337032,0.35947,0.369323,0.390588,0.397284,0.107452,0.178166,0.19498,0.218967,
	0.251285,0.274866,0.287295,0.304418,0.312874,0.328455,0.347097,0.369123,0.383139,0.404648,
	0.418848,0.439186,0.448115,0.160502,0.167729,0.208132,0.235969,0.255649,0.273104,0.282156,
	0.292743,0.306677,0.325346,0.333234,0.340617,0.354114,0.362994,0.374284,0.384236,0.397615,
	0.122517,0.173391,0.212969,0.226017,0.249399,0.257665,0.275866,0.295666,0.301476,0.305935,
	0.324604,0.338282,0.376935,0.397586,0.420077,0.429148,0.439342,0.146627,0.18519,0.215078,
	0.253996,0.269057,0.281222,0.297689,0.307957,0.330116,0.352681,0.367494,0.375439,0.39794,
	0.410285,0.424246,0.435729,0.453332,0.170265,0.194879,0.234531,0.270917,0.293073,0.307564,
	0.314507,0.333272,0.3497,0.352121,0.36727,0.384403,0.394541,0.412919,0.42767,0.443735,
	0.451425,0.141126,0.177858,0.234747,0.255151,0.274855,0.288148,0.300371,0.31379,0.32754,
	0.351761,0.366432,0.386421,0.407136,0.419659,0.431686,0.457918,0.470111,0.149923,0.174644,
	0.197092,0.236323,0.263329,0.281669,0.29105,0.297966,0.329311,0.355383,0.366161,0.385621,
	0.399445,0.412529,0.426503,0.431387,0.440149,0.0800145,0.138678,0.164506,0.234237,0.256929,
	0.265171,0.293889,0.300691,0.326141,0.337138,0.354166,0.362909,0.378904,0.395488,0.410782,
	0.427502,0.441895,0.134481,0.157854,0.177082,0.208822,0.231677,0.247232,0.269131,0.28039,
	0.306232,0.32259,0.338923,0.352059,0.367572,0.379844,0.383661,0.396363,0.401896,0.110722,
	0.134867,0.170162,0.183199,0.228039,0.245007,0.261427,0.27318,0.283768,0.304406,0.310975,
	0.329279,0.360268,0.374565,0.391752,0.411125,0.429402,0.101687,0.165464,0.208786,0.226307,
	0.243687,0.248975,0.286018,0.301468,0.315596,0.324623,0.338518,0.344203,0.360828,0.373918,
	0.396109,0.406549,0.414514,0.126005,0.153,0.205003,0.222491,0.251502,0.267004,0.280232,
	0.296427,0.322846,0.334958,0.34352,0.357828,0.363008,0.379448,0.394273,0.403915,0.434501,
	0.171276,0.21166,0.23394,0.275293,0.298253,0.310992,0.329059,0.340245,0.349497,0.362799,
	0.380313,0.39856,0.409754,0.422129,0.431921,0.440535,0.448178,0.095994,0.167235,0.204207,
	0.221051,0.247439,0.277378,0.291544,0.315528,0.327202,0.354091,0.368442,0.375527,0.386373,
	0.410704,0.423064,0.439788,0.446244,0.17521,0.207647,0.237518,0.245324,0.263182,0.273782,
	0.287366,0.302517,0.310782,0.324088,0.330511,0.357191,0.362749,0.386814,0.413168,0.4274,
	0.435504,0.144333,0.167683,0.233047,0.251958,0.280602,0.292576,0.305382,0.315318,0.332846,
	0.348909,0.361805,0.379305,0.390702,0.394901,0.407814,0.428001,0.452755,0.15149,0.16887,
	0.215958,0.224958,0.240594,0.266098,0.286118,0.295569,0.316947,0.330042,0.343859,0.351846,
	0.359638,0.372771,0.399002,0.427983,0.43527,0.154237,0.186336,0.201139,0.229845,0.267336,
	0.280731,0.286708,0.30371,0.331119,0.355971,0.378128,0.386137,0.397969,0.422344,0.43592,
	0.445377,0.476323,0.126878,0.169464,0.194318,0.215073,0.232171,0.267155,0.277833,0.298943,
	0.318473,0.326046,0.333748,0.345387,0.355628,0.370192,0.380475,0.402919,0.420005,0.105944,
	0.155164,0.183592,0.210032,0.251364,0.260218,0.293509,0.312432,0.319093,0.338167,0.359003,
	0.377474,0.384076,0.40111,0.426871,0.455342,0.468271,0.12853,0.1716,0.212701,0.231142,
	0.242722,0.259707,0.291283,0.309241,0.327542,0.349888,0.367441,0.386394,0.417207,0.426134,
	0.442031,0.472049,0.484847,0.116511,0.160902,0.185523,0.208517,0.227164,0.255335,0.28034,
	0.295692,0.308324,0.316198,0.33412,0.34389,0.364056,0.383989,0.397602,0.407816,0.420838,
	0.14102,0.21655,0.233815,0.256958,0.270011,0.281732,0.294086,0.306958,0.321959,0.330551,
	0.350648,0.364735,0.370392,0.380508,0.393217,0.407787,0.424272,0.111315,0.161844,0.212348,
	0.227761,0.235386,0.249042,0.258482,0.275538,0.295966,0.313748,0.327774,0.353991,0.368047,
	0.376177,0.404187,0.429931,0.43162,0.0921289,0.134095,0.171704,0.199914,0.23029,0.241072,
	0.253931,0.271792,0.300041,0.31155,0.326102,0.333513,0.361519,0.381523,0.398326,0.424972,
	0.435111,0.129205,0.161745,0.17365,0.20758,0.261526,0.27555,0.291034,0.29973,0.324773,
	0.340285,0.363313,0.376933,0.39589,0.40274,0.411361,0.427181,0.440307,0.133302,0.187836,
	0.203969,0.233746,0.254548,0.272062,0.283171,0.30043,0.319237,0.330325,0.35269,0.362463,
	0.372526,0.38829,0.417213,0.447368,0.460206,0.104081,0.149678,0.180593,0.2053,0.215205,
	0.238992,0.251275,0.265524,0.290347,0.295594,0.304997,0.318963,0.336727,0.353796,0.359937,
	0.377334,0.400111,0.108229,0.165908,0.19847,0.222532,0.264097,0.284268,0.295874,0.312794,
	0.331847,0.341866,0.373305,0.387603,0.399041,0.414709,0.434664,0.440666,0.448808,0.119767,
	0.16609,0.188209,0.230584,0.252668,0.267083,0.288471,0.31319,0.325652,0.343862,0.367071,
	0.379026,0.387402,0.39708,0.412145,0.419317,0.431462,0.167878,0.185395,0.211356,0.223693,
	0.237904,0.258839,0.276718,0.292701,0.301774,0.313579,0.332859,0.368396,0.380461,0.393463,
	0.403317,0.410669,0.42356,0.145044,0.170042,0.210483,0.244667,0.261349,0.274108,0.286501,
	0.307949,0.325912,0.357127,0.368403,0.381374,0.397555,0.411255,0.424069,0.442434,0.448239,
	0.0979732,0.172811,0.189148,0.22562,0.25198,0.269296,0.279077,0.297874,0.317108,0.337144,
	0.357638,0.36328,0.378136,0.392775,0.40892,0.415119,0.431807,0.149189,0.173689,0.205869,
	0.223509,0.232954,0.258752,0.271394,0.28447,0.305805,0.320351,0.342637,0.366588,0.388868,
	0.399163,0.404764,0.420848,0.444713,0.117016,0.152456,0.181553,0.228682,0.261749,0.280662,
	0.29234,0.309637,0.321458,0.337131,0.364598,0.378604,0.391999,0.405544,0.415174,0.425551,
	0.441398,0.123356,0.137056,0.178879,0.210323,0.23963,0.256064,0.282427,0.296938,0.317371,
	0.331644,0.335862,0.344934,0.361796,0.373921,0.382376,0.391557,0.404396,0.0887843,0.119549,
	0.171894,0.185293,0.224161,0.236917,0.251656,0.280974,0.287476,0.304284,0.31592,0.3332,
	0.361122,0.377494,0.402478,0.419428,0.428984,0.139629,0.172956,0.203317,0.252863,0.279087,
	0.292083,0.300012,0.325713,0.346144,0.356779,0.363138,0.377757,0.392265,0.399075,0.415369,
	0.425001,0.439237,0.127579,0.174025,0.205196,0.242894,0.269451,0.282529,0.293913,0.308667,
	0.317417,0.328651,0.351556,0.364834,0.383547,0.395859,0.407078,0.42123,0.443129,0.143805,
	0.181698,0.190338,0.20488,0.23774,0.260342,0.281216,0.290577,0.302948,0.330408,0.349241,
	0.356209,0.360445,0.380799,0.393577,0.406774,0.418652,0.151692,0.170018,0.193165,0.22353,
	0.254601,0.276232,0.298962,0.310689,0.321933,0.328378,0.342388,0.367126,0.380799,0.393165,
	0.418912,0.435348,0.447996,0.181263,0.207006,0.218692,0.236329,0.263523,0.270001,0.305276,
	0.31838,0.336841,0.343051,0.358149,0.382484,0.39198,0.405261,0.417049,0.427502,0.433514,
	0.104675,0.172936,0.218625,0.236659,0.255369,0.270988,0.281836,0.290162,0.307205,0.326767,
	0.339015,0.345348,0.364301,0.374655,0.391672,0.426711,0.449357,0.117019,0.137522,0.200169,
	0.221517,0.249068,0.25663,0.284616,0.295914,0.303893,0.313167,0.327407,0.354145,0.374738,
	0.380887,0.387684,0.401905,0.418495,0.0959244,0.175627,0.199141,0.223644,0.245533,0.273396,
	0.296328,0.308205,0.320982,0.329033,0.340035,0.355362,0.36351,0.383424,0.405369,0.419777,
	0.42758,0.127623,0.163129,0.180614,0.200729,0.214685,0.23047,0.261134,0.282415,0.303543,
	0.325651,0.342053,0.351917,0.358872,0.378091,0.389402,0.406642,0.41466,0.12696,0.221378,
	0.246872,0.262523,0.273393,0.303117,0.325572,0.352206,0.364785,0.371523,0.381304,0.388797,
	0.395909,0.412808,0.426143,0.43224,0.451873,0.120921,0.171943,0.203076,0.23002,0.256145,
	0.272795,0.293393,0.297375,0.314097,0.329657,0.34436,0.367647,0.384509,0.397236,0.411581,
	0.421332,0.432051,0.163033,0.20467,0.243339,0.275509,0.28176,0.295901,0.304209,0.321991,
	0.34065,0.359187,0.376925,0.385328,0.40825,0.416648,0.421853,0.438592,0.463028,0.118533,
	0.163438,0.196345,0.223052,0.229582,0.247638,0.266901,0.285173,0.299904,0.321675,0.337885,
	0.354903,0.368348,0.374718,0.383399,0.39845,0.417502,0.159311,0.215943,0.230677,0.25459,
	0.270915,0.289926,0.299207,0.310201,0.323841,0.334935,0.342722,0.363436,0.382905,0.39213,
	0.396179,0.40616,0.417289,0.12541,0.154593,0.173258,0.207287,0.237859,0.255962,0.273738,
	0.288365,0.314998,0.327935,0.334959,0.363727,0.371745,0.396359,0.40878,0.427952,0.446523
	};
	
	
	//for every test point:
	unsigned counter = 0;
	for ( size_t pts = 0 ; pts != test_points.size() ; ++pts )
	{
		//for every distance class:
		for ( size_t dist = 0 ; dist != classes_to_compute_knn.size() ; ++dist )
		{
			//std::cout << (*classes_to_compute_knn[dist])( test_points[pts]) << ",";
			//if ( counter % 10 == 9 )  std::cout << std::endl;
			BOOST_CHECK( fabs( (*classes_to_compute_knn[dist])( test_points[pts]) - result[counter] ) <= 5e-06 );
			++counter;
		}
	}		
}	




/*

 //to test Clement's cde, remoporary test, to be remove d later
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

	
	//Gudhi::Topological_inference_with_cubical_complexes::Distance_to_k_th_closest_point_periodic_k_d_tree
	//kd_tree( point_cloud, coordinates_of_grid , 50 );	
	//Gudhi::Topological_inference_with_cubical_complexes::Distance_to_k_th_closest_point<periodic_Euclidean_distance_squared>
	//brute_force( point_cloud, period_eu , 50 );
	
	
	std::vector< Gudhi::Topological_inference_with_cubical_complexes::Distance_to_k_th_closest_point_periodic_k_d_tree* >
	kd_tree;
	std::vector< Gudhi::Topological_inference_with_cubical_complexes::Distance_to_k_th_closest_point<periodic_Euclidean_distance_squared>* >
	brute_force;
			
	for ( size_t i = 5 ; i != 90 ; i=i+5 )
	{
		kd_tree.push_back(
		new Gudhi::Topological_inference_with_cubical_complexes::Distance_to_k_th_closest_point_periodic_k_d_tree( point_cloud, coordinates_of_grid , i ) );
		
		brute_force.push_back(
		new Gudhi::Topological_inference_with_cubical_complexes::Distance_to_k_th_closest_point<periodic_Euclidean_distance_squared>( point_cloud, period_eu , i ) );
	}  
	
	std::cerr << "Here \n";
	    
	
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
		
	//for every test point:
	unsigned counter = 0;
	for ( size_t pts = 0 ; pts != test_points.size() ; ++pts )
	{		
		//if ( fabs(brute_force( test_points[pts]) - kd_tree( test_points[pts])) >= 5e-05 )
		//{
		//	std::cout << "test point : " << test_points[pts][0] << " " << test_points[pts][1] << " " << test_points[pts][2] << std::endl;
		//	std::cerr << "brute_force( test_points[pts]) : " << brute_force( test_points[pts]) << std::endl;
		//	std::cerr << "k_d_tree( test_points[pts]) : " << kd_tree( test_points[pts]) << std::endl;
		//	getchar();
		//}
		for ( size_t dist = 0 ; dist != kd_tree.size() ; ++dist )
		{			
			if ( fabs( (*kd_tree[dist])( test_points[pts]) - (*brute_force[dist])( test_points[pts]) ) >= 5e-05 )
			{
				std::cout << (*kd_tree[dist])( test_points[pts]) << " " << (*brute_force[dist])( test_points[pts]) << std::endl;
				getchar();
			}
			++counter;
		}	
	}		
}	
*/


