#include <gudhi/SparseMsMatrix.h>
#include <gudhi/Fake_simplex_tree.h>
#include <gudhi/Simplex_tree.h>

#include <gudhi/Rips_complex.h>
#include <gudhi/distance_functions.h>

#include <gudhi/Points_off_io.h>
#include <gudhi/pick_n_random_points.h>

#include <CGAL/Epick_d.h>
#include <CGAL/point_generators_d.h>
#include <CGAL/Random.h>

#include <boost/program_options.hpp>

#include <cmath>
#include <chrono>

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <limits>  // for std::numeric_limits

using Point = CGAL::Epick_d< CGAL::Dynamic_dimension_tag>::Point_d;
using Filtration_value = Fake_simplex_tree::Filtration_value;
using Rips_complex = Gudhi::rips_complex::Rips_complex<Filtration_value>;

using Fake_simplex_tree = Gudhi::Fake_simplex_tree ;
using Simplex_tree = Gudhi::Simplex_tree<>;
using Vector_of_points = std::vector<Point>;

void program_options(int argc, char * const argv[]
                     , std::size_t & number_of_points
                     , double & begin_thresold
                     , double & steps
                     , double & end_thresold
                     , int 	  & repetetions
                     , char   & manifold
                     , int 	  & dimension
                     , std::string & in_file_name
                     , std::string & out_file_name
                     ) {
  namespace po = boost::program_options;

  po::options_description visible("Allowed options", 100);
  visible.add_options()
    ("help,h", "produce help message")
    ("number,n", po::value<std::size_t>(&number_of_points)->default_value(0),
       "Number of generated point_vector.")
    
    ("begin_thresold,b", po::value<double>(&begin_thresold)->default_value(0),
       "Initial threshold for rips complex.")
  	("steps,s", po::value<double>(&steps)->default_value(0.1),
   		"Steps of the threshold")
    ("end_thresold,e", po::value<double>(&end_thresold)->default_value(1),
  		"Final threshold for rips complex.")
      
    ("repetetions,r", po::value<int>(&repetetions)->default_value(1),
    	"Num of repetetions of the experiments.")
  	("manifold,m", po::value<char>(&manifold)->default_value('s'),
   		"Steps of the threshold")
       
    ("dimensions,D", po::value<int>(&dimension)->default_value(2),
     "Dimension of the manifold.")

    ("input_file_name,i", po::value<std::string>(&in_file_name),
     "The input file.")
    ("out_file_name,o", po::value<std::string>(&out_file_name),
     "The output file.");

  po::options_description all;
  all.add(visible);

  po::variables_map vm;
  po::store(po::command_line_parser(argc, argv).
            options(all).run(), vm);
  po::notify(vm);

  if (vm.count("help")) {
    std::cout << std::endl;
    std::cout << "Computes rips complexes of different threshold values, from 'begin_thresold' to 'end_thresold', with priodic steps of 'steps' from a n random uniform point_vector on a selected manifold, . \n";
    std::cout << "Strongly collapses all the rips complexes and output the results in out_file. \n";
    std::cout << "The experiments are repeted 'repete' num of times for each threshold value. \n";
    std::cout << "type -m for manifold options, 's' for uni sphere, 'b' for unit ball, 'f' for file. \n";
    std::cout << "type -i 'filename' for Input file option for exported point sample. \n";
    std::cout << std::endl << std::endl;

    std::cout << "Usage: " << argv[0] << " [options]" << std::endl << std::endl;
    std::cout << visible << std::endl;
    std::abort();
  }
}

void generate_points_sphere(Vector_of_points& W, int nbP, int dim, int radius)
{
	CGAL::Random_points_on_sphere_d<Point> rp(dim, radius);
	for (int i = 0; i < nbP; i++)
		W.push_back(*rp++);
}
void generate_points_ball(Vector_of_points& W, int nbP, int dim, int radius)
{
	CGAL::Random_points_in_ball_d<Point> rp(dim, radius); 
	for (int i = 0; i < nbP; i++)
		W.push_back(*rp++);
}

void generate_points_cube(Vector_of_points& W, int nbP, int dim)
{
	CGAL::Random_points_in_cube_d<Point> rp(dim, 6);
	for (int i = 0; i < nbP; i++)
		W.push_back(*rp++);
}

void add_point_vectors(Vector_of_points& V, Vector_of_points& U, int nbP, int dim) // Adds two point vectors of the same size (nbP), by Modifying the first one, V = V+W.
{
	for (int i = 0; i < nbP; i++)
	{
		std::vector<double> coords;
		for(int j =0; j< dim; j++)
			coords.push_back(V[i][j]+U[i][j]); 
		V[i] = Point(coords);
	}
}


int main(int argc, char * const argv[])
{
	
	std::string out_file_name = "default";
	std::string in_file_name = "default";
	std::size_t number_of_points;
  	
  	int 	dimension;
  	double 	begin_thresold;
  	double 	end_thresold;
  	double 	steps;
  	int    	repetetions;
	char   	manifold;

	Vector_of_points *point_vector;
	Vector_of_points file_all_points;
	
	int dime = 20000; // pseudo variable... of no use
	std::string manifold_full = "sphere";

	int  originalDim 		= 0;
	int  collDim 			= 0;
	// long originalNumSimp 	= 0;

	long colNumSimp 		= 0;
	long originalMaxNumSimp = 0;
	long colMaxNumSimp 		= 0;

	int  colNumVert 		= 0;
	long originalNumMxSimp 	= 0;
	long colNumMxSimp 		= 0;
	double collapseTime 	= 0;
	
	int radius  = 1;
	long originalNumSimpFrmSmplxtr 	= 0;

	program_options(argc, argv, number_of_points, begin_thresold, steps, end_thresold, repetetions, manifold, dimension, in_file_name, out_file_name);
    std::ofstream myfileDetl (out_file_name); 

    if(manifold == 'f' || manifold =='F')
	{
	  	Gudhi::Points_off_reader<Point> off_reader(in_file_name);
	  	if (!off_reader.is_valid()) 
	  	{
	      	std::cerr << "Sparse matrix for Rips complex - Unable to read file " << in_file_name << "\n";
	      	exit(-1);  // ----- >>
	    }
	  	
	  	file_all_points = Vector_of_points(off_reader.get_point_cloud());

	  	std::cout << "Successfully read " << file_all_points.size() << " point_vector.\n";
	  	std::cout << "Ambient dimension is " << file_all_points[0].dimension() << ".\n";
	}
   
   	std::cout << "The current input values to run the program is: "<< std::endl;
	std::cout << "number_of_points, begin_thresold, steps, end_thresold, repetetions, manifold, dimension, in_file_name, out_file_name" << std::endl;
    std::cout << number_of_points << ", " << begin_thresold << ", " << steps << ", " << end_thresold << ", " << repetetions << ", " << manifold << ", " << dimension<< ", " << in_file_name << ", " << out_file_name << std::endl;

	myfileDetl << "Thresold, InitialNumVertices, AfterCollNumVert, InitialDimension, CollDimension, InitialNumSimp, CollapsedNumSimplices, TimeInMS" << std::endl;

	Fake_simplex_tree 	* stree;      
	Fake_simplex_tree 	* coll_tree; 	
	SparseMsMatrix 		* mat;

	Simplex_tree 		*simplexTreeToCnt;

	double threshold =  begin_thresold;

	while(threshold <= end_thresold)
	{
		std::cout <<  "Process started for threshold value: " << threshold << std::endl;
		
		for(int i = 0; i < repetetions; i++)
		{
			std::cout << "Current iteration/repetetion is: " << (i+1) << " and threshold: " << threshold << std::endl;
			
			point_vector = new Vector_of_points();

			if(manifold == 's' || manifold == 'S'){
				generate_points_sphere(*point_vector, number_of_points, dimension, radius);
				std::cout << number_of_points << " points successfully chosen randomly from "<< dimension <<"-sphere of radius " << radius << std::endl;
		    }
			else if(manifold == 'b' || manifold == 'B'){
				generate_points_ball(*point_vector, number_of_points, dimension, radius); 
				std::cout << number_of_points << " points successfully chosen randomly from "<< dimension <<"-ball of radius " << radius << std::endl;
			
			}
			else if(manifold == 'f' || manifold =='F')
			{
			  	// Subsampling from all points for each iterations
  				Gudhi::subsampling::pick_n_random_points(file_all_points, number_of_points, std::back_inserter(*point_vector));
  				std::cout << number_of_points << " points succesfully chosen randomly from "<< dimension <<"input points " << std::endl;
			}
			else
			{
				std::cerr << "Wrong choice for input manifold..." <<std::endl;	
				exit(-1); 
			}

			
			stree 		= new Fake_simplex_tree();
			coll_tree 	= new Fake_simplex_tree();

			simplexTreeToCnt = new Simplex_tree();
			
			Rips_complex rips_complex_from_points(*point_vector, threshold, Gudhi::Euclidean_distance());
			
			// std::cout << "Formation of the toplex map began" << std::endl;
			// auto toplex_form_began = std::chrono::high_resolution_clock::now();
			rips_complex_from_points.create_complex(*stree, dime);
			// auto toplex_form_end= std::chrono::high_resolution_clock::now();
			// auto toplexFormTime = std::chrono::duration<double, std::milli>(toplex_form_end- toplex_form_began).count();
			
			// std::cout << "Formation of the simplex tree began" << std::endl;
			// auto simplx_tree_form_began = std::chrono::high_resolution_clock::now();
			rips_complex_from_points.create_complex(*simplexTreeToCnt,dime);
			// auto simplx_tree_form_end = std::chrono::high_resolution_clock::now();
			// auto simplexTreeFormTime = std::chrono::duration<double, std::milli>(simplx_tree_form_end- simplx_tree_form_began).count();
			
			// std::cout << " Formation eneded and Counting by simplex tree began" << std::endl;
			// auto simplx_tree_count_began = std::chrono::high_resolution_clock::now();
			
			// auto simplx_tree_count_end = std::chrono::high_resolution_clock::now();
			// std::cout << "Counting by simplex tree ended" << std::endl;
			// auto stree_formed  = std::chrono::high_resolution_clock::now();
			// auto simplexTreeCountTime = std::chrono::duration<double, std::milli>(simplx_tree_count_end- simplx_tree_count_began).count();
			
			std::cout << "Toplex map created... Next stop matrix formation" << std::endl;

			mat = new SparseMsMatrix(*stree);

			auto matrix_formed  = std::chrono::high_resolution_clock::now();
			std::cout << "Matrix formed ... Next action COLLAPSE!!" << std::endl;

			mat->strong_collapse();
			auto collapse_done = std::chrono::high_resolution_clock::now();
			std::cout << "Collapse done !" << std::endl;
			collapseTime += std::chrono::duration<double, std::milli>(collapse_done- matrix_formed).count();
			std::cout << "Time of Collapse And Simplex tree formation : ~" << collapseTime/(i+1) << " ms\n" << std::endl;

			coll_tree = new Fake_simplex_tree(mat->collapsed_tree());
	
			originalDim 		+= mat->initial_dimension();
			collDim 			+= mat->collapsed_dimension();
			
			colNumVert			+= coll_tree->num_vertices();

			originalNumMxSimp 	+= stree->num_simplices();
			colNumMxSimp 		+= coll_tree->num_simplices();

			originalMaxNumSimp 	+= mat->max_num_inital_simplices();
			colMaxNumSimp	   	+= mat->max_num_collapsed_simplices();
			
			// auto toplex_count_began = std::chrono::high_resolution_clock::now();
			// originalNumSimp 	+= ((stree->filtration_simplex_range().size()) - 1);
			// auto toplex_count_end = std::chrono::high_resolution_clock::now();
			// auto toplexCountTime = std::chrono::duration<double, std::milli>(toplex_count_end- toplex_count_began).count();

			originalNumSimpFrmSmplxtr 	+= simplexTreeToCnt->num_simplices();
			colNumSimp					+= ((coll_tree->filtration_simplex_range().size()) - 1);

			delete stree;
			delete coll_tree;
			delete mat;
			delete point_vector;
			delete simplexTreeToCnt;

			// std::cout << "Toplex formation time: " << toplexFormTime <<" Toplex count time: " << toplexCountTime << " Simplex tree formation time: " << simplexTreeFormTime << " Simplex tree count time: "
			// 		 << simplexTreeCountTime << std::endl;	

		}

		originalDim 		= originalDim/repetetions;
		collDim 			= collDim/repetetions;
		colNumVert			= colNumVert/repetetions;

		originalNumMxSimp 	= originalNumMxSimp/repetetions ;
		colNumMxSimp 		= colNumMxSimp/repetetions ;

		originalMaxNumSimp 	= originalMaxNumSimp/repetetions;
		colMaxNumSimp	   	= colMaxNumSimp/repetetions;

		// originalNumSimp 	= originalNumSimp/repetetions ;
		colNumSimp			= colNumSimp/repetetions;
		collapseTime   		= collapseTime/repetetions;

		originalNumSimpFrmSmplxtr = originalNumSimpFrmSmplxtr/repetetions;

		std::cout << "Rips complex is of dimension " << originalDim << " with " << originalNumMxSimp << "  maximal simplices and " << originalNumSimpFrmSmplxtr  << " simplices and " << number_of_points << " vertices." << std::endl;
		std::cout << "Collapsed Rips complex is of dimension " << collDim << " with " <<  colNumMxSimp << " maximal simplices and " <<  colNumSimp <<  " simplices and " << colNumVert << " vertices." << std::endl;

		std::cout << "** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** " << std::endl;

		myfileDetl << threshold << "," << number_of_points << "," << colNumVert << "," << originalDim << "," << collDim << "," << originalNumSimpFrmSmplxtr  <<  "," << colNumSimp << "," << collapseTime << std::endl;

		threshold = threshold+steps;

		originalDim 		= 0;
		collDim 			= 0;
		colNumVert			= 0;
		originalNumMxSimp 	= 0;
		colNumMxSimp 		= 0;
		originalMaxNumSimp 	= 0;
		colMaxNumSimp	   	= 0;
		// originalNumSimp 	= 0;
		colNumSimp			= 0;
		collapseTime   		= 0;
		originalNumSimpFrmSmplxtr = 0;
	}

	myfileDetl.close();
	return 0;
}
