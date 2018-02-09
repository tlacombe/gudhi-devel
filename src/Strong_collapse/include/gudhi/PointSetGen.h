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

class PointSetGen 
{
  public:
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
    
    /*
      Generates point sets on <nbSpheres> spheres wedged at origin, sphere can have different radii from <init_radius> with steps of <multiplier_step>
      Number of points on the sphere can also be different from <init_nbP> for the smallest sphere and then multiplied by <multiplier_step>
    */
    void generate_points_wedged_sphere(Vector_of_points& W, int init_nbP, int dim, int init_radius, int multiplier_step, int nbSpheres)
    {
      int radius = init_radius;
      int nbP = init_nbP;
      std::vector<double> translation;
      for(int d = 0; d< dim; d++){
        translation.push_back(0);
      }

      for(int s = 0; s < nbSpheres; s++)
      {
        CGAL::Random_points_on_sphere_d<Point> rp(dim, radius); 
        for (int i = 0; i < nbP; i++)
        {
          W.push_back(add_point(*rp++, translation, dim)); 
        }
        nbP = nbP*multiplier_step;
        radius = radius*multiplier_step;
        translation.at(dim-1) = (radius - init_radius);
        
      }
      
    }

    void generate_points_concentric_sphere(Vector_of_points& W, int init_nbP, int dim, int init_radius, int multiplier_step, int nbSpheres)
    {
      int radius = init_radius;
      int nbP = init_nbP;
     
      for(int s = 0; s < nbSpheres; s++)
      {
        CGAL::Random_points_on_sphere_d<Point> rp(dim, radius); 
        for (int i = 0; i < nbP; i++)
        {
          W.push_back(*rp++); 
        }
        nbP = nbP*multiplier_step;
        radius = radius*multiplier_step;
      }
      
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
    		V[i] = add_point(V[i], U[i], dim);
    	}
    }

    //returns x = x+y;
    Point add_point(const Point & x, const Point & y, int dim)
    {
      std::vector<double> coords;
      for(int i =0; i< dim; i++)
        coords.push_back(x[i]+y[i]); 
      return Point(coords);
    }
    void print_point(const Point & x, int dim)
    { std::cout<< "(";
      for(int i = 0; i < dim; i++){
        std::cout<< x[i] << ", " ;
      }
      std::cout<< ")" << std::endl;  
    }
    
    PointSetGen(){}

    ~PointSetGen(){}
};