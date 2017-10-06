/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Pawel Dlotko
 *
 *    Copyright (C) 2017 Swansea University
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

#include <gudhi/reader_utils.h>
#include <gudhi/Off_reader.h>

#include <gudhi/Topological_inference.h>
#include <gudhi/Bitmap_cubical_complex.h>
#include <gudhi/Persistent_cohomology.h>

// standard stuff
#include <iostream>
#include <sstream>
#include <vector>
#include <sstream>

#include <gudhi/Points_off_io.h>

using Point = std::vector<double>;
using Points_off_reader = Gudhi::Points_off_reader<Point>;


int main()
{
	//This example show how to create a function, which is distance to k-th nearest neighbor and use it for topological inference.
	
	//Below we have a sampling of S^1 that we will be the point cloud used here:
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
	
	//First we need to decide the rectagle that will encapsulate our cubical grid. In this case,
	//we are picking up a square [-1,1]^2
	std::vector< std::pair< double,double > > coorfinates_of_grid(2);
	coorfinates_of_grid[0] = coorfinates_of_grid[1] = std::pair<double,double>( -2.0 , 2.0 );
	
	//Later we need to decide on the resolution. Here our two dimensional grid will have 
	//100 maximal cubes resolution in each direction (giving 100*100=10000 grid points)
	std::vector< unsigned > resolution_of_a_grid(2);	
	resolution_of_a_grid[0] = resolution_of_a_grid[1] = 100;
	
	
	
	//A few typedefs to make the code easy to read:
	
	
	//********************************************************************************************
	
/*	
	//Use this set of typedefs if you want to use distance to k-th nearest neighbor as your function:
	
	//Some distance functions require additional parameters. For instance, distance to the k-th
	//nearest neighbor require the parameter k, which is known here as number_of_nearest_neighbors
	unsigned number_of_nearest_neighbors = 5;	
	//distance function typedefs
    typedef Gudhi::Topological_inference_with_cubical_complexes::Euclidan_distance_squared local_distance;    
    typedef Gudhi::Topological_inference_with_cubical_complexes::Distance_to_k_th_closest_point<local_distance> actual_distance;    
    local_distance eu;
    actual_distance f( point_cloud ,eu ,  number_of_nearest_neighbors );
*/  
    
    
    //********************************************************************************************    
    
    
 /*   
    //Use this set of typedefs if you want to use kernels_centerd_in_point_cloud class as your function.
    //Class kernels_centerd_in_point_cloud is a class that for every point x under considertation, iterate 
    //through the initial point cloud, and for every point p in this point cloud, compute a kernel (which is a template
    //parameter of the kernels_centerd_in_point_cloud of x and p. The overal value of a kernel is the accumulated value
    //for all points. 
    //There is a number of kernels you can choose over here. Please look for them at:
    //functions_for_topological_inference/functions_for_topological_inference.h
    //Sample examples tested here are: Euclidan_distance_squared, Manhattan_distance and Max_norm_distance
    //They can be used both in periodic and non-periodic wersion regardless of the (periodic or not periodic)
    //cubical complex that is going to be created for topological_inference.
    
    //distance function typedefs
    typedef Gudhi::Topological_inference_with_cubical_complexes::Euclidan_distance_squared local_distance;
    typedef Gudhi::Topological_inference_with_cubical_complexes::kernels_centerd_in_point_cloud<local_distance> actual_distance;
    local_distance eu;
    actual_distance f( eu , point_cloud );
 */ 
    
      //********************************************************************************************    
    
    /*
    //to use Sum_of_distances_from_points, please use those templates:   
    typedef Gudhi::Topological_inference_with_cubical_complexes::Sum_of_distances_from_points actual_distance;
    actual_distance f( point_cloud );
	*/
    
   
    //Use the code below to constrct periodic version of Manhattan_distance distance
    //on the periodic grid [-1,1]^2
    //define the periodic grid    
    std::vector< std::pair< double , double > > coordinates_of_grid(2);
    coordinates_of_grid[0] = coordinates_of_grid[1] = std::pair<double,double>( -1.0,1.0 );
    
    //Periodic version of Manhattan_distance
    Gudhi::Topological_inference_with_cubical_complexes::Manhattan_distance manhattan;
    typedef Gudhi::Topological_inference_with_cubical_complexes::periodic_domain_distance< Gudhi::Topological_inference_with_cubical_complexes::Manhattan_distance > local_distance;   
    local_distance periodic_sum( coordinates_of_grid,manhattan );
    
    //now we have distance to the k-th closest points 
    unsigned number_of_nearest_neighbors = 5;
    typedef Gudhi::Topological_inference_with_cubical_complexes::Distance_to_k_th_closest_point
    < local_distance  > actual_distance;
    
    actual_distance f( point_cloud , periodic_sum , number_of_nearest_neighbors );
    
   
  
    
    
  
	//topological inference typedefs
    typedef Gudhi::Cubical_complex::Bitmap_cubical_complex_base<double> Bitmap_cubical_complex_base;
    typedef Gudhi::Cubical_complex::Bitmap_cubical_complex<Bitmap_cubical_complex_base> Bitmap_cubical_complex;
    typedef Gudhi::Topological_inference_with_cubical_complexes::Topological_inference< Bitmap_cubical_complex , double ,   actual_distance > topological_inference;
  
    //computations of persistence typedefs
    typedef Gudhi::persistent_cohomology::Field_Zp Field_Zp;
    typedef Gudhi::persistent_cohomology::Persistent_cohomology<topological_inference, Field_Zp> Persistent_cohomology;

	//now we can create our topological inference object
    topological_inference b( coorfinates_of_grid , resolution_of_a_grid , f );
    
    //you can vizualize some of the 2d functions with gnuplot. Proper output is obtained by using this format.
    //later call: load 'to_view' matrix with image
    b.write_to_file_with_newlines_at_the_ends_of_structure("to_view");
    

    //And compute its persistence diagram
    Persistent_cohomology pcoh(b);
    pcoh.init_coefficients(2); 
    
    //we want only the intervals of persistence >= 0.1
    pcoh.compute_persistent_cohomology(0.1);
    
    //and print it out to the screen: 
    std::vector< std::tuple<size_t, size_t, int> > intervals = pcoh.get_persistent_pairs();   
    for ( size_t i = 0 ; i != intervals.size() ; ++i )
    {
		std::cout << "Dimension : " << b.get_dimension_of_a_cell(std::get<0>(intervals[i])) << " birth: " << b.filtration(std::get<0>(intervals[i])) 
		<< " death: " << b.filtration(std::get<1>(intervals[i])) << std::endl;
	}
	
	return 0;
}
