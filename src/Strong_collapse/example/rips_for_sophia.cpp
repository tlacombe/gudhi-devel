#include <gudhi/FormatTower.h>
#include <gudhi/Rips_complex.h>
#include <gudhi/distance_functions.h>

#include <CGAL/Epick_d.h>
#include <CGAL/point_generators_d.h>
#include <CGAL/Random.h>
#include <cmath>

#include <iostream>
#include <string>
#include <vector>
#include <limits>

using Point = CGAL::Epick_d< CGAL::Dimension_tag<10> >::Point_d;
using Fake_simplex_tree     = Gudhi::Fake_simplex_tree ;

using Filtration_value = Fake_simplex_tree::Filtration_value;
using Rips_complex = Gudhi::rips_complex::Rips_complex<Filtration_value>;

using Vector_of_points = std::vector<Point>;

void generate_points_sphere(Vector_of_points& W, int nbP, int dim) {
	CGAL::Random_points_on_sphere_d<Point> rp(dim, 1);
	for (int i = 0; i < nbP; i++)
		W.push_back(*rp++);
}

int main() {
	
    Vector_of_points points;
    Map map_empty;
	Fake_simplex_tree subComplex;
    Fake_simplex_tree *st_prev      = new Fake_simplex_tree() ;
    Fake_simplex_tree *st_prev_coll = new Fake_simplex_tree() ;
    
    

    std::size_t num_pts = 100;
    int maxSplxDim = 10;
    int sphereDim = 2;
	
    FormatTower towerFormater;
    FormatTower origTowerFormater;
    generate_points_sphere(points,num_pts,sphereDim); // 2 for S2
	
    int numSimComp = 15;
    // double threshold_vals[3] = {0.1,0.2,0.4};// {0.25,0.5,0.75,1,1.1};
    double threshold_vals = 0.6;
    double incremental_step = 0.01;
     
    for(int i = 0 ; i < numSimComp ; ++i)
    {
	 	
        Rips_complex rips_complex_from_points(points, threshold_vals, Gudhi::Euclidean_distance());
        rips_complex_from_points.create_complex(subComplex, maxSplxDim);
        std::cout << "Simplex tree for parameter " << threshold_vals << " created." <<std::endl;
        
        if( i == numSimComp-1)
        {
            origTowerFormater.one_step_tower(*st_prev, subComplex, map_empty, "OriginalTowerRips.txt");
            std::cout << "Original Tower computed" << std::endl;
            // delete st_prev;
            // st_prev = new Fake_simplex_tree(subComplex);
        }
        
        SparseMsMatrix mat(subComplex);
        mat.strong_collapse();
	 	Fake_simplex_tree* subCpx_coll = new Fake_simplex_tree(mat.collapsed_tree());
        std::cout << "Subcomplex #" << (i+1) << " Collapsed" << std::endl;
	 	Map redmap = mat.reduction_map();

        std::cout << "Rips complex is of dimension " << mat.initial_dimension() << " with " << subComplex.num_simplices() << " maximal simplices " << std::endl;
        std::cout << "Collapsed Rips complex is of dimension " << mat.collapsed_dimension() << " with " <<  subCpx_coll->num_simplices() << " maximal simplices" << std::endl;

        std::cout << "** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** " << std::endl;
        
        towerFormater.one_step_tower(*st_prev_coll, *subCpx_coll, redmap, "CollapsedTowerRips.txt");
        std::cout << "Tower updated for subcomplex #" << i << std::endl; 
        delete st_prev_coll;
        st_prev_coll = subCpx_coll; 
        threshold_vals = threshold_vals+incremental_step;
    }

    delete st_prev_coll;
    //delete st_prev;
    return 0;
}
