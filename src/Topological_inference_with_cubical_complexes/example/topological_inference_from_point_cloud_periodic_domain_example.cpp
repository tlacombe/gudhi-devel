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
	//This example show how to create a function, which is distance to k-th nearest neighbor on a periodic domain.
	
	//Below we have a sof a x axis from -1 to 1. The sampling is embeded in R^2:
	std::vector< std::vector< double > > point_cloud =
	{ 
	{-0.037521644,0},{-0.0825336012,0},{0.7911890061,0},{-0.9246191038,0},{0.4524056804,0},{0.0775374882,0},
	{0.9467886547,0},{0.8274683105,0},{0.2227936597,0},{-0.9703822481,0},{-0.5710002943,0},{0.4003993953,0},
	{0.2283410849,0},{0.8515946642,0},{-0.5128273047,0},{0.4887397168,0},{-0.811148053,0},{-0.8842573576,0},
	{-0.2220111005,0},{-0.5831138957,0},{0.8204358765,0},{-0.5207901569,0},{-0.5540787992,0},{-0.4000139898,0},
	{0.7256436953,0},{0.2093629842,0},{-0.3356358311,0},{-0.1762533984,0},{0.9359343066,0},{0.8788811448,0},
	{-0.866219101,0},{-0.529315426,0},{0.5705546658,0},{-0.7242854806,0},{0.5188311604,0},{0.729877892,0},
	{0.0464174165,0},{0.5329394881,0},{0.2439277298,0},{-0.6725864992,0},{-0.3394736839,0},{0.0889316835,0},
	{-0.4414633382,0},{0.5508343112,0},{0.6158041162,0},{-0.4060751721,0},{0.6979859914,0},{0.7339886054,0},
	{0.6153524141,0},{0.4652338373,0},{-0.4786614878,0},{0.260527797,0},{-0.2397670946,0},{-0.8819530639,0},
	{0.8997337194,0},{0.375566321,0},{-0.2439842103,0},{0.3730698815,0},{0.5866884915,0},{0.9692299939,0},
	{-0.1300423318,0},{0.5369295818,0},{-0.3530211342,0},{0.8829268264,0},{-0.6039790777,0},{0.4238893455,0},
	{-0.3912836928,0},{0.6745120543,0},{0.4432661091,0},{0.8482065122,0},{-0.5663325754,0},{0.7848620885,0},
	{-0.7805975657,0},{-0.966094452,0},{-0.1505859527,0},{0.2879623915,0},{-0.6182084312,0},{0.7471420313,0},
	{-0.1730242288,0},{-0.6181439795,0},{-0.878390268,0},{0.524438641,0},{-0.3151322999,0},{-0.5649360898,0},
	{0.9157870538,0},{0.5519235665,0},{0.5554750767,0},{0.6858953899,0},{0.9699861645,0},{-0.8798311735,0},
	{0.0980144818,0},{-0.8763509346,0},{-0.7061346853,0},{0.6577362157,0},{0.5012289365,0},{-0.4615081623,0},
	{-0.0274199871,0},{0.9074162804,0},{-0.6656602919,0},{-0.1186002814,0}
	};	
	
	//First we need to decide the rectagle that will encapsulate our cubical grid. In this case,
	//we are picking up a square [-2,2]^2
	std::vector< std::pair< double,double > > coorfinates_of_grid(2);
	coorfinates_of_grid[0] = coorfinates_of_grid[1] = std::pair<double,double>( -1.0 , 1.0 );
	
	//Later we need to decide on the resolution. Here our two dimensional grid will have 
	//100 maximal cubes resolution in each direction (giving 100*100=10000 grid points)
	std::vector< unsigned > resolution_of_a_grid(2);	
	resolution_of_a_grid[0] = resolution_of_a_grid[1] = 100;
	
	//Then we need to decide in which direction the periodic boundary conditions are to be imposed.
	//In this case, we impose them in both directions
	std::vector<bool> directions_in_which_periodic_b_cond_are_to_be_imposed(2);
	directions_in_which_periodic_b_cond_are_to_be_imposed[0] = directions_in_which_periodic_b_cond_are_to_be_imposed[1] = true;
	
	//Lastly, some distance functions require additional parameters. For instance, distance to the k-th
	//nearest neighbor require the parameter k, which is known here as number_of_nearest_neighbors
	unsigned number_of_nearest_neighbors = 5;
	
	//A few typedefs to make the code easy to read:
	//distance function typedefs
    Gudhi::Topological_inference_with_cubical_complexes::Euclidan_distance_squared eu;
    Gudhi::Topological_inference_with_cubical_complexes::Distance_to_k_th_closest_point<Gudhi::Topological_inference_with_cubical_complexes::Euclidan_distance_squared> 
    f( point_cloud ,eu ,  number_of_nearest_neighbors );
  
	//topological inference typedefs
	typedef Gudhi::Cubical_complex::Bitmap_cubical_complex_periodic_boundary_conditions_base<double> Periodic_bitmap_cubical_complex_base;
    typedef Gudhi::Cubical_complex::Bitmap_cubical_complex<Periodic_bitmap_cubical_complex_base> Periodic_bitmap_cubical_complex;
    typedef Gudhi::Topological_inference_with_cubical_complexes::Topological_inference< Periodic_bitmap_cubical_complex , double ,   
    Gudhi::Topological_inference_with_cubical_complexes::Distance_to_k_th_closest_point<Gudhi::Topological_inference_with_cubical_complexes::Euclidan_distance_squared> > topological_inference;
  
  
    //computations of persistence typedefs
    typedef Gudhi::persistent_cohomology::Field_Zp Field_Zp;
    typedef Gudhi::persistent_cohomology::Persistent_cohomology<topological_inference, Field_Zp> Persistent_cohomology;

	//now we can create our topological inference object
    topological_inference b( coorfinates_of_grid , resolution_of_a_grid , f , directions_in_which_periodic_b_cond_are_to_be_imposed );
        

    //And compute its persistence diagram. The second parameter (set to true) will allow us to see the
    //infinite generator in dimension 2.
    Persistent_cohomology pcoh(b,true);
    pcoh.init_coefficients(2); 
    
    //we want only the intervals of persistence >= 0.1
    pcoh.compute_persistent_cohomology(0.1);
    
    
    //and print it out to the screen. Note that one of the infinite generator in dimension 1 
    //is born very close to zero, which is not the case for the other one. This is because 
    //of the point cloud which is sampled from x axis, and which glue the two sides of periodic 
    //domanin very fast.
    std::vector< std::tuple<size_t, size_t, int> > intervals = pcoh.get_persistent_pairs();   
    for ( size_t i = 0 ; i != intervals.size() ; ++i )
    {
		std::cout << "Dimension : " << b.get_dimension_of_a_cell(std::get<0>(intervals[i])) << " birth: " << b.filtration(std::get<0>(intervals[i])) 
		<< " death: " << b.filtration(std::get<1>(intervals[i])) << std::endl;
	}
	
	return 0;
}
