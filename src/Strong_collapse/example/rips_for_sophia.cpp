#include <gudhi/FormatTower.h>
#include <gudhi/Rips_complex.h>
#include <gudhi/distance_functions.h>

#include <CGAL/Epick_d.h>
#include <CGAL/point_generators_d.h>
#include <CGAL/Random.h>

#include <gudhi/PointSetGen.h>

#include <cmath>
#include <iostream>
#include <string>
#include <vector>
#include <limits>
#include <algorithm>

using Fake_simplex_tree   = Gudhi::Fake_simplex_tree ;
using Filtration_value    = Fake_simplex_tree::Filtration_value;
using Rips_complex        = Gudhi::rips_complex::Rips_complex<Filtration_value>;

using Vector_of_points = std::vector<Point>;

int main() {
	
    Vector_of_points points;
    Map map_empty;

    std::size_t init_num_pts = 50;
    int maxSplxDim = 10;
    int sphereDim = 2;          // ( actual dimension is sphereDim -1)
    int radius    = 1;
    int multiplier_steps = 2;
    int num_spheres = 4;
    int tot_num_points = init_num_pts*(pow(multiplier_steps, num_spheres) - 1);

    Fake_simplex_tree *subComplex ;   
    SparseMsMatrix * mat_coll       = new SparseMsMatrix(); 
    SparseMsMatrix * mat_prev_coll  = new SparseMsMatrix(tot_num_points, 100*tot_num_points); 

    FormatTower towerFormater(tot_num_points) ;
    FormatTower origTowerFormater(tot_num_points);
    PointSetGen point_generator;
    // point_generator.generate_points_sphere(points,tot_num_points,sphereDim,radius); // 2 for S2
   
    // point_generator.generate_points_wedged_sphere(points, init_num_pts, sphereDim, radius, multiplier_steps, num_spheres);
    point_generator.generate_points_concentric_sphere(points, init_num_pts, sphereDim, radius, multiplier_steps, num_spheres);
    std::cout << "Point Set Generated: "  <<std::endl;
    // for(int i = 0; i <tot_num_points; i++ )
    //     point_generator.print_point(points.at(i), sphereDim);

    int numSimComp = 15;
    double threshold_vals = 0.3;
    double incremental_step = 0.01;
     
    for(int i = 0 ; i < numSimComp ; ++i)
    {
	 	subComplex = new Fake_simplex_tree();

        Rips_complex rips_complex_from_points(points, threshold_vals, Gudhi::Euclidean_distance());
        rips_complex_from_points.create_complex(*subComplex, maxSplxDim);
        std::cout << "Rips complex for parameter " << threshold_vals << " created." <<std::endl;
        
        mat_coll = new SparseMsMatrix(*subComplex);

        if( i == numSimComp-1)
        {
            origTowerFormater.one_step_tower(* new SparseMsMatrix(tot_num_points, 100*tot_num_points), *mat_coll, map_empty, "OriginalTowerRips.txt");
            std::cout << "Original Tower computed" << std::endl;
        }

        mat_coll->strong_collapse();
        std::cout << "Subcomplex #" << (i+1) << " Collapsed" << std::endl;
        Map redmap = mat_coll->reduction_map();
         
        

        std::cout << "Uncollapsed Rips complex is of dimension " << mat_coll->initial_dimension() << " with " << subComplex->num_simplices() << " maximal simplices " << std::endl;
        std::cout << "Collapsed Rips complex is of dimension " << mat_coll->collapsed_dimension() << " with " <<  mat_coll->number_max_simplices() << " maximal simplices" << std::endl;

        std::cout << "** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** " << std::endl;
        
        towerFormater.one_step_tower(*mat_prev_coll, *mat_coll, redmap, "CollapsedTowerRips.txt");
        std::cout << "Tower updated for subcomplex #" << i+1 << std::endl; 
        delete mat_prev_coll;
        delete subComplex;
        mat_prev_coll = mat_coll; 
        threshold_vals = threshold_vals+incremental_step;
    }

    delete mat_prev_coll;
    return 0;
}
