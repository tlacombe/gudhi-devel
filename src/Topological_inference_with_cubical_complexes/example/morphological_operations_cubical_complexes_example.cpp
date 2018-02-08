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
#include <gudhi/Topological_inference.h>
#include <gudhi/Bitmap_cubical_complex.h>
#include <gudhi/Persistent_cohomology.h>

// standard stuff
#include <iostream>
#include <sstream>
#include <vector>
#include <sstream>

#include <gudhi/Points_off_io.h> 
#include <gudhi/Morphological_operations_cubical_complex.h>


int main()
{
	//This file is a demonstration how to use a class Morphological_operations_cubical_complex.
	//In most typical case, the class take as an input a filtered cubical complex and a predicate. It then
	//pick up the top dimensional cubes which satisfy the predicate - they will constitute
	//a 'set'. The cubes that do not satisfy the predicate constitute set's 'complement'.
	//Then a distance function (in terms of a minimal number of hops from a cube belonging to 
	//a set to a cube in set's complement is computed and re-efined on this cubical complex.
	//When it comes to define neighborhood of a cube, at the moment two types of neighborhood are
	//implemented; 'full_face' (neighs of a maximal cube cube C consist of all maximal cubes D such 
	//that C intersect with D is a top dimensional face. The other type is 'all', where neighberhood 
	//of C consist of all cubes having nonempty intersection with C. 
	//Later on, we can compute persistent homology of the cubical complex with this filtration.
	
	//There are various ways we can have the initial complex defined. In this demonstration we present the 
	//following options:
	//(A) We start from a sublevel set of a function and treshold it.
	//(B) we manually pick up the cells we want to be in the set.
	//(C) We use topological_inference class to define the right function.
	//When presenting those cases we will use different types of neighborhoods, and different options as
	//of having periodic or non-periodic grid.
	//Please consilt the details below:

{	
	//CASE A:
	//We start by defining a function on a cubical complex by hand. In this case, this is a distance from a cicle in a plane,
	//in a 9 by 9 grid:
	//4	3	3	3	3	3	3	3	4
	//3	2	1	1	1	1	1	2	3
	//3	1	0	0	0	0	0	1	3
	//3	1	0	1	1	1	0	1	3
	//3	1	0	1	2	1	0	1	3
	//3	1	0	1	1	1	0	1	3
	//3	1	0	0	0	0	0	1	3
	//3	2	1	1	1	1	1	2	3
	//4	3	3	3	3	3	3	3	4


    std::cout << std::endl << std::endl << std::endl << "Case A "  << std::endl << std::endl << std::endl;
	std::vector< unsigned > sizes(2);
	sizes[0] = sizes[1] = 9;
	std::vector<double> data = {
								4,3,3,3,3,3,3,3,4,
								3,2,1,1,1,1,1,2,3,
								3,1,0,0,0,0,0,1,3,
								3,1,0,1,1,1,0,1,3,
								3,1,0,1,2,1,0,1,3,
								3,1,0,1,1,1,0,1,3,
								3,1,0,0,0,0,0,1,3,
								3,2,1,1,1,1,1,2,3,
								4,3,3,3,3,3,3,3,4 
							};
	
	//In case you would like to use periodic cubical complexes, set those variables and change the 
	//typedef in the next paragraph from the ones that use Bitmap_cubical_complex_base into the ones
	//that use Bitmap_cubical_complex_periodic_boundary_conditions_base 
	std::vector<bool> directions_in_which_periodic_boundary_conditions_are_to_be_imposed(2);
	directions_in_which_periodic_boundary_conditions_are_to_be_imposed[0] =
	directions_in_which_periodic_boundary_conditions_are_to_be_imposed[1] = false;
	
	//construct a cubical complex based on the data. Here is the code for standarc cubical complex:
    typedef Gudhi::Cubical_complex::Bitmap_cubical_complex_base<double> Bitmap_cubical_complex_base;
    typedef Gudhi::Cubical_complex::Bitmap_cubical_complex<Bitmap_cubical_complex_base> Bitmap_cubical_complex;
    //and here for the periodic domain. If you want to use the periodic domain, please comment out the
    //two lines above, and uncomment the four lines below:
    //typedef Gudhi::Cubical_complex::Bitmap_cubical_complex_periodic_boundary_conditions_base<double> Bitmap_cubical_complex_base_periodic;
    //typedef Gudhi::Cubical_complex::Bitmap_cubical_complex<Bitmap_cubical_complex_base_periodic> Bitmap_cubical_complex;
    
    
    Bitmap_cubical_complex b( sizes , data , directions_in_which_periodic_boundary_conditions_are_to_be_imposed);
    
    //now we construct an object of Morphological_operations_cubical_complex type:
		
    //This is the preictor, which defines set elements as those below certain value (given as a parameter).
    //Experiment with changing the 0.5 cutoff value. 
    typedef Gudhi::Topological_inference_with_cubical_complexes::Filtration_below_certain_value<double> Predictor_type;
    Predictor_type pred(0.5);        
    
    //This is a Morphological_operations_cubical_complex build based on Bitmap_cubical_complex and Predictor_type
    typedef Gudhi::Topological_inference_with_cubical_complexes::Morphological_operations_cubical_complex<Bitmap_cubical_complex,Predictor_type> MOCC;
       
	MOCC mor( &b , pred );	
	size_t index = 0;
	std::cout << "Here is the complex after running predicate : " << std::endl;
	for ( auto it = b.top_dimensional_cells_iterator_begin() ; it != b.top_dimensional_cells_iterator_end() ; ++it )
	{		
		std::cout << b.get_cell_data( *it ) << " ";		
		if ( index % 9 == 8 )std::cout << std::endl;
		++index;
	}		
	
	//now we can run the dilation of the set. When experimenting, please change the type of neighborhood
	//from 'all' to 'full_face' and observe the difference in the dylated complex.
	mor.dilation( 1 , Gudhi::Topological_inference_with_cubical_complexes::considered_neighberhoods::all );
	
	//display the values on top dimensional cubes of the complex. Note that I am adding newlines just for display
	//pusposes:
	index = 0;
	std::cout << "Here is the complex after running dilation procedure : " << std::endl;
	for ( auto it = b.top_dimensional_cells_iterator_begin() ; it != b.top_dimensional_cells_iterator_end() ; ++it )
	{		
		std::cout << b.get_cell_data( *it ) << " ";		
		if ( index % 9 == 8 )std::cout << std::endl;
		++index;
	}
	
	//and compute persistence of this complex:
	typedef Gudhi::persistent_cohomology::Field_Zp Field_Zp;
    typedef Gudhi::persistent_cohomology::Persistent_cohomology<Bitmap_cubical_complex, Field_Zp> Persistent_cohomology;
  
    Persistent_cohomology pcoh(b,true);
    pcoh.init_coefficients(2);  // initializes the coefficient field for homology
    pcoh.compute_persistent_cohomology(0);
        
    std::vector< std::tuple<size_t, size_t, int> > intervals = pcoh.get_persistent_pairs();       

	std::cout << "And here are the persistence intervals : \n";
    for ( size_t i = 0 ; i != intervals.size() ; ++i )
    {
		std::cout << "dimenion : " << b.get_dimension_of_a_cell(std::get<0>(intervals[i])) << " " << b.filtration(std::get<0>(intervals[i])) << " " << b.filtration(std::get<1>(intervals[i])) << std::endl;
	}
}

    //************************************************************************************************************
	//************************************************************************************************************
	//************************************************************************************************************
	//************************************************************************************************************
	//************************************************************************************************************
	
{	
	//CASE B:
	std::cout << std::endl << std::endl << std::endl << "Case B "  << std::endl << std::endl << std::endl;
	//In this case we will start from 2 dimensional 5 by 5 complex. 	
	std::vector< unsigned > sizes1(2);
	sizes1[0] = sizes1[1] = 5;
	
	//we may want to construct a cubical complex with periodic boundary conditions.
	//If you set those to true, remember to change the typedefs of cubical complexes from
	//Bitmap_cubical_complex_base to Bitmap_cubical_complex_periodic_boundary_conditions_base
	//in the next paragraph of code.
	std::vector<bool> directions_in_which_periodic_boundary_conditions_are_to_be_imposed1(2);
	directions_in_which_periodic_boundary_conditions_are_to_be_imposed1[0] =
	directions_in_which_periodic_boundary_conditions_are_to_be_imposed1[1] = true;
	
	//construct a cubical complex based on the data. Here is the code for standarc cubical complex:
    typedef Gudhi::Cubical_complex::Bitmap_cubical_complex_base<double> Bitmap_cubical_complex_base;
    typedef Gudhi::Cubical_complex::Bitmap_cubical_complex<Bitmap_cubical_complex_base> Bitmap_cubical_complex;
    //and here for the periodic domain. If you want to use the periodic domain, please comment out the
    //two lines above, and uncomment the four lines below:
    //typedef Gudhi::Cubical_complex::Bitmap_cubical_complex_periodic_boundary_conditions_base<double> Bitmap_cubical_complex_base_periodic;
    //typedef Gudhi::Cubical_complex::Bitmap_cubical_complex<Bitmap_cubical_complex_base_periodic> Bitmap_cubical_complex;
        
    
    //Now we need to pick up the cells we want to put into a set. They will form an annulus 
    //centered in the center of bitmap. To do so, we can start from taking the positions of
    //the maximal cubes we want to take. Here is the filtration of the top dimensional cubes 
    //in the complex:
    
    //0 0 0 0 0
    //0 1 1 1 0
    //0 1 0 1 0
    //0 1 1 1 0
    //0 0 0 0 0
    std::vector< std::vector<unsigned> > maximal_cubes_we_select = 
    {
		{1,1},
		{2,1},
		{3,1},
		{1,2},
		{3,2},
		{1,3},
		{2,3},
		{3,3}
	};
    
    //Now we need to translate maximal_cubes_we_select into real positions in bitmap:
   
		
	//This is a Morphological_operations_cubical_complex build based on Bitmap_cubical_complex and Predictor_type
	typedef Gudhi::Topological_inference_with_cubical_complexes::Filtration_below_certain_value<double> Predictor_type;
    typedef Gudhi::Topological_inference_with_cubical_complexes::Morphological_operations_cubical_complex<Bitmap_cubical_complex,Predictor_type> MOCC;
       
	MOCC mor1( maximal_cubes_we_select , sizes1 , directions_in_which_periodic_boundary_conditions_are_to_be_imposed1 );			
	
	//Let us print out the complex:
	Bitmap_cubical_complex* cmplx = mor1.get_complex();
	size_t index = 0;
	std::cout << "Here are the values on top dimensional cubes: \n";
	for ( auto it = cmplx->top_dimensional_cells_iterator_begin() ; it != cmplx->top_dimensional_cells_iterator_end() ; ++it )
	{		
		std::cout << cmplx->get_cell_data( *it ) << " ";		
		if ( index % 5 == 4 )std::cout << std::endl;
		++index;
	}
	
	
	//now we can run the dilation of the set. When experimenting, please change the type of neighborhood
	//from 'all' to 'full_face' and observe the difference in the dylated complex.
	mor1.dilation( 1 , Gudhi::Topological_inference_with_cubical_complexes::considered_neighberhoods::all );
	
	//display the values on top dimensional cubes of the complex. Note that I am adding newlines just for display
	//pusposes:
	cmplx = mor1.get_complex();
	index = 0;
	std::cout << "Here are the values after dilation: \n";
	for ( auto it = cmplx->top_dimensional_cells_iterator_begin() ; it != cmplx->top_dimensional_cells_iterator_end() ; ++it )
	{		
		std::cout << cmplx->get_cell_data( *it ) << " ";		
		if ( index % 5 == 4 )std::cout << std::endl;
		++index;
	}
	
	//and compute persistence of this complex:
	typedef Gudhi::persistent_cohomology::Field_Zp Field_Zp;
    typedef Gudhi::persistent_cohomology::Persistent_cohomology<Bitmap_cubical_complex, Field_Zp> Persistent_cohomology;
  
    Persistent_cohomology pcoh1(*cmplx,true);
    pcoh1.init_coefficients(2);  // initializes the coefficient field for homology
    pcoh1.compute_persistent_cohomology(0);
        
    std::vector< std::tuple<size_t, size_t, int> > intervals1 = pcoh1.get_persistent_pairs();       

	std::cout << "And here are the persistence intervals : \n";
    for ( size_t i = 0 ; i != intervals1.size() ; ++i )
    {
		std::cout << "dimenion : " << cmplx->get_dimension_of_a_cell(std::get<0>(intervals1[i])) << ", birth : " 
		<< cmplx->filtration(std::get<0>(intervals1[i])) << ", death: " << cmplx->filtration(std::get<1>(intervals1[i])) 
		<< std::endl;
	}
}

	//************************************************************************************************************
	//************************************************************************************************************
	//************************************************************************************************************
	//************************************************************************************************************
	//************************************************************************************************************
	
{	
	//CASE C:
	std::cout << std::endl << std::endl << std::endl << "Case C "  << std::endl << std::endl << std::endl;	
	//first we will start from having those points sampled from unit circle.
	std::vector< std::vector<double> > point_cloud = 
	{
	{0,1},{0.0998334166,0.9950041653},{0.1986693308,0.9800665778},{0.2955202067,0.9553364891},{0.3894183423,0.921060994},
	{0.4794255386,0.8775825619},{0.5646424734,0.8253356149},{0.6442176872,0.7648421873},{0.7173560909,0.6967067093},
	{0.7833269096,0.6216099683},{0.8414709848,0.5403023059},{0.8912073601,0.4535961214},{0.932039086,0.3623577545},
	{0.9635581854,0.2674988286},{0.98544973,0.1699671429},{0.9974949866,0.0707372017},{0.999573603,-0.0291995223},
	{0.9916648105,-0.1288444943},{0.9738476309,-0.2272020947},{0.9463000877,-0.3232895669},{0.9092974268,-0.4161468365},
	{0.8632093666,-0.5048461046},{0.8084964038,-0.5885011173},{0.7457052122,-0.6662760213},{0.6754631806,-0.7373937155},
	{0.5984721441,-0.8011436155},{0.5155013718,-0.8568887534},{0.4273798802,-0.904072142},{0.3349881502,-0.9422223407},
	{0.2392493292,-0.9709581651},{0.1411200081,-0.9899924966},{0.0415806624,-0.9991351503},{-0.0583741434,-0.9982947758},
	{-0.1577456941,-0.9874797699},{-0.255541102,-0.9667981926},{-0.3507832277,-0.9364566873},{-0.4425204433,-0.8967584163},
	{-0.5298361409,-0.8481000317},{-0.6118578909,-0.7909677119},{-0.6877661592,-0.7259323042},{-0.7568024953,-0.6536436209},
	{-0.8182771111,-0.5748239465},{-0.8715757724,-0.4902608213},{-0.9161659367,-0.4007991721},{-0.9516020739,-0.30733287},
	{-0.9775301177,-0.2107957994},{-0.9936910036,-0.1121525269},{-0.9999232576,-0.0123886635},{-0.9961646088,0.0874989834},
	{-0.9824526126,0.1865123694},{-0.9589242747,0.2836621855},{-0.9258146823,0.3779777427},{-0.8834546557,0.4685166713},
	{-0.8322674422,0.5543743362},{-0.7727644876,0.6346928759},{-0.7055403256,0.7086697743},{-0.6312666379,0.7755658785},
	{-0.5506855426,0.8347127848},{-0.4646021794,0.8855195169},{-0.3738766648,0.9274784307},{-0.2794154982,0.9601702867},
	{-0.1821625043,0.9832684384},{-0.0830894028,0.996542097},{0.0168139005,0.9998586364}
	};
	
	
	//now we build a topological inference object based on the point_cloud:
	//now creat an object of topoogical inference based on it:
	std::vector< std::pair< double,double > > coorfinates_of_grid(2);
	coorfinates_of_grid[0] = coorfinates_of_grid[1] = std::pair<double,double>( -2,2 );
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
	mor.dilation(1 , Gudhi::Topological_inference_with_cubical_complexes::considered_neighberhoods::all );

    // Compute the persistence diagram of the complex
    Persistent_cohomology pcoh(b);
    pcoh.init_coefficients(2);  // initializes the coefficient field for homology
    pcoh.compute_persistent_cohomology(0.2);
    
 
    std::cout << "Here are the persistence intervals obtained as a result of dilation: \n";
    std::vector< std::tuple<size_t, size_t, int> > intervals = pcoh.get_persistent_pairs();   
    for ( size_t i = 0 ; i != intervals.size() ; ++i )
    {
		std::cout << "dimenion : " << b.get_dimension_of_a_cell(std::get<0>(intervals[i])) << ", birth : " 
		<< b.filtration(std::get<0>(intervals[i])) << ", death: " << b.filtration(std::get<1>(intervals[i])) 
		<< std::endl;		
	}
}
	
	
	return 0;	
}
