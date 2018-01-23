#include <gudhi/SparseMsMatrix.h>
#include <gudhi/Fake_simplex_tree.h>

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
    std::cout << "The experiments are repeted 'repete' num of times. \n";
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
  	
  	int dimension;
  	double begin_thresold;
  	double end_thresold;
  	double steps;
  	int    repetetions;
	char   manifold;

	Vector_of_points point_vector;
	Vector_of_points file_all_points;
	
	int dime = 10; // pseudo variable... of no use
	std::string manifold_full = "sphere";

	int originalDim = 0;
	int collDim = 0;
	long originalNumSimp = 0;
	long colNumSimp = 0;
	long originalMaxNumSimp = 0;
	long colMaxNumSimp = 0;

	int  colNumVert = 0;
	long originalNumMxSimp =0;
	long colNumMxSimp = 0;
	double collapseTime = 0;
	
	// int dimBall	= 2;
	int radius  = 1;
	


	program_options(argc, argv, number_of_points, begin_thresold, steps, end_thresold, repetetions, manifold, dimension, in_file_name, out_file_name);
    std::ofstream myfileDetl (out_file_name); 

    if(manifold == 'f' || manifold =='F')
	{
	  	Gudhi::Points_off_reader<Point> off_reader(in_file_name);
	  	if (!off_reader.is_valid()) {
	      	std::cerr << "Sparse matrix for Rips complex - Unable to read file " << in_file_name << "\n";
	      	exit(-1);  // ----- >>
	    	}
	  	file_all_points = Vector_of_points(off_reader.get_point_cloud());

	  	std::cout << "Successfully read " << file_all_points.size() << " point_vector.\n";
	  	std::cout << "Ambient dimension is " << file_all_points[0].dimension() << ".\n";
	}
   	//add_point_vectors(point_vector, noisePoints, number_of_points);
	
	// std::cout << "Number of point_vector : " << number_of_points << " and " << "Threshold value : " << threshold << std::endl;
	// Rips_complex rips_complex_from_points();
	
	// std::ofstream myfileOrig ("./output/rips1_original.xls");
	// std::ofstream myfileCol  ("./output/rips1_collapsed.xls");
	

	// myfileOrig << "Threshold, Dimension, NumOfSimplices" << std::endl;
	// myfileOrig << "Threshold, Dimension, NumOfSimplices"<< std::endl;
	
	myfileDetl << "Thresold, InitialNumVertices, AfterCollNumVert, InitialDimension, CollDimension, InitialNumSimp, CollapsedNumSimplices, TimeInMS" << std::endl;


	double threshold =  begin_thresold;
	Fake_simplex_tree stree;

	while(threshold <= end_thresold)
	{

		for(int i = 0; i < repetetions; i++){
			
			if(manifold == 's' || manifold == 'S'){
				generate_points_sphere(point_vector,number_of_points,dimension, radius);
				myfileDetl << number_of_points << " point_vector chosen randomly from "<< dimension <<"-sphere of radius " << radius << std::endl;
		    }
			else if(manifold == 'b' || manifold == 'B'){
				generate_points_ball(point_vector,number_of_points,dimension,radius); 
				myfileDetl << number_of_points << " point_vector chosen randomly from "<< dimension <<"-ball of radius " << radius << std::endl;
			
			}
			else if(manifold == 'f' || manifold =='F')
			{
			  	// Subsampling from all points for each iterations
  				Gudhi::subsampling::pick_n_random_points(file_all_points, number_of_points, std::back_inserter(point_vector));
  				myfileDetl << number_of_points << " point_vector chosen randomly from "<< dimension <<"input points " << std::endl;
			}
			
			Rips_complex rips_complex_from_points(point_vector, threshold, Gudhi::Euclidean_distance());
			rips_complex_from_points.create_complex(stree, dime);

			// auto stree_formed  = std::chrono::high_resolution_clock::now();
			std::cout << "Simplex tree created ... Next stop matrix formation" << std::endl;

			SparseMsMatrix mat(stree);
			auto matrix_formed  = std::chrono::high_resolution_clock::now();
			std::cout << "Matrix formed ... Next action COLLAPSE!!" << std::endl;

			Fake_simplex_tree coll_tree = mat.collapsed_tree();
			auto collapse_done = std::chrono::high_resolution_clock::now();

			std::cout << "Collapse done !" << std::endl;
			collapseTime += std::chrono::duration<double, std::milli>(collapse_done- matrix_formed).count();
			// std::cout << "Time for formation of Matrix : " << (matrix_formed - stree_formed)/CLOCKS_PER_SEC << " seconds" << std::endl;
			std::cout << "Time for Collapse : " << collapseTime << " ms\n" << std::endl;
			
			originalDim += mat.initial_dimension();
			collDim 	+= mat.collapsed_dimension();
			
			colNumVert	+= coll_tree.num_vertices();

			originalNumMxSimp 	+= stree.num_simplices();
			colNumMxSimp 		+= coll_tree.num_simplices();

			originalMaxNumSimp += mat.max_num_inital_simplices();
			colMaxNumSimp	   += mat.max_num_collapsed_simplices();

			originalNumSimp += (stree.filtration_simplex_range().size() - 1);
			colNumSimp		+= (coll_tree.filtration_simplex_range().size() - 1);
		}

		std::cout << "Rips complex is of dimension " << originalDim << " with " << originalNumMxSimp/repetetions << " maximal simplices, " << originalNumSimp/repetetions << "/" << originalMaxNumSimp/repetetions << " simplices and " << number_of_points << " vertices." << std::endl;
		std::cout << "Collapsed Rips complex is of dimension " << collDim << " with " <<  colNumMxSimp/repetetions << " maximal simplices, " <<  colNumSimp/repetetions << "/" << colMaxNumSimp/repetetions << " simplices and " << colNumVert/repetetions << " vertices." << std::endl;

		std::cout << "** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** " << std::endl;


		// myfileOrig << threshold << "," << originalDim << "," << originalNumSimp << std::endl;
		// myfileOrig << threshold << "," << collDim << "," << colNumMxSimp << std::endl;
		myfileDetl << threshold << "," << number_of_points << "," << colNumVert/repetetions << "," << originalDim/repetetions << "," << collDim/repetetions << "," << originalNumSimp/repetetions  <<  "," << colNumSimp/repetetions << "," << collapseTime/repetetions << std::endl;

		threshold = threshold+steps;
	}
	// myfileOrig.close();
	// myfileCol.close();
	myfileDetl.close();

	return 0;
}
