#include <gudhi/FormatTower.h>

#include <gudhi/Rips_complex.h>
#include <gudhi/distance_functions.h>

#include <gudhi/PointSetGen.h>

using Fake_simplex_tree   = Gudhi::Fake_simplex_tree ;
using Filtration_value    = Fake_simplex_tree::Filtration_value;
using Rips_complex        = Gudhi::rips_complex::Rips_complex<Filtration_value>;

using Vector_of_points = std::vector<Point>;

int main(int argc, char * const argv[]) {
	
    PointSetGen point_generator;
    std::string out_file_name = "default";
    std::string in_file_name = "default";
    std::size_t number_of_points;
    
    int     dimension;
    double  begin_thresold;
    double  end_thresold;
    double  steps;
    int     repetetions = 1;
    char    manifold;

    Vector_of_points *point_vector;
    Vector_of_points file_all_points;
    
    int dime = 20000; // pseudo variable... of no use
    std::string manifold_full = "sphere";
    
    double radius  = 1;
    double r_min = 0.8;
    double r_max = 1;

    point_generator.program_options(argc, argv, number_of_points, begin_thresold, steps, end_thresold, repetetions, manifold, dimension, in_file_name, out_file_name);
    
    if(manifold == 'f' || manifold =='F') {
        Gudhi::Points_off_reader<Point> off_reader(in_file_name);
        if (!off_reader.is_valid()) {
            std::cerr << "Sparse matrix for Rips complex - Unable to read file " << in_file_name << "\n";
            exit(-1);  // ----- >>
        }

        file_all_points = Vector_of_points(off_reader.get_point_cloud());
        std::cout << "Successfully read " << file_all_points.size() << " point_vector.\n";
        std::cout << "Ambient dimension is " << file_all_points[0].dimension() << ".\n";
    }
   
    std::cout << "The current input values to run the program is: "<< std::endl;
    std::cout << "number_of_points, begin_thresold, steps, end_thresold, repetetions, manifold, dimension, in_file_name, out_file_name" << std::endl;
    std::cout << number_of_points << ", " << begin_thresold << ", " << steps << ", " << end_thresold << ", " << repetetions << ", " << manifold << ", " << dimension << ", " << in_file_name << ", " << out_file_name << std::endl;

    Map map_empty;

    std::string origFile = "./PersistenceOutput/original_tower_rips" ;
    std::string collFile = "./PersistenceOutput/collapsed_tower_rips" ;

    double totCollapseTime = 0.0;
    double totPrintTime = 0.0;

    Fake_simplex_tree * subComplex ;   

    SparseMsMatrix * mat_coll       = new SparseMsMatrix(); 
    SparseMsMatrix * mat_prev       = new SparseMsMatrix(number_of_points, 100*number_of_points); 
    SparseMsMatrix * mat_prev_coll  = new SparseMsMatrix(number_of_points, 100*number_of_points); 

    FormatTower towerFormater(number_of_points) ;
    FormatTower origTowerFormater(number_of_points);
   
    point_vector = new Vector_of_points();

    if(manifold == 's' || manifold == 'S'){
        point_generator.generate_points_sphere(*point_vector, number_of_points, dimension, radius);
        origFile = origFile+"_"+"_sphere_"+out_file_name+".txt";
        collFile = collFile+"_"+"_sphere_"+out_file_name+".txt";
        std::cout << number_of_points << " points successfully chosen randomly from "<< dimension <<"-sphere of radius " << radius << std::endl;
    }
    else if(manifold == 'b' || manifold == 'B'){
        point_generator.generate_points_ball(*point_vector, number_of_points, dimension, radius); 
        origFile = origFile+"_"+"_ball_"+out_file_name+".txt";
        collFile = collFile+"_"+"_ball_"+out_file_name+".txt";
        std::cout << number_of_points << " points successfully chosen randomly from "<< dimension <<"-ball of radius " << radius << std::endl;
    
    }
    else if( (manifold == 'a' || manifold == 'A')&& dimension == 2){
        point_generator.generate_points_2annulus(*point_vector, number_of_points, r_min, r_max); 
        origFile = origFile+"_"+"_annulus_"+out_file_name+".txt";
        collFile = collFile+"_"+"_annulus_"+out_file_name+".txt";
        std::cout << number_of_points << " points successfully chosen randomly from "<< 2 <<"-annulus of radii (" << r_min << ',' << r_max << ") " << std::endl;

    }
    else if( (manifold == 'a' || manifold == 'A') && dimension == 3){
        point_generator.generate_points_spherical_shell(*point_vector, number_of_points, r_min, r_max); 
        origFile = origFile+"_"+"_annulus_"+out_file_name+".txt";
        collFile = collFile+"_"+"_annulus_"+out_file_name+".txt";
        std::cout << number_of_points << " points successfully chosen randomly from spherical shell of radii (" << r_min << ',' << r_max << ") " << std::endl;

    }
    
    else if(manifold == 'f' || manifold =='F') {
        // Subsampling from all points for each iterations
        Gudhi::subsampling::pick_n_random_points(file_all_points, number_of_points, std::back_inserter(*point_vector));
        origFile = origFile+"_"+ in_file_name +"-"+ out_file_name+ ".txt";
        collFile = collFile+"_"+ in_file_name +"-"+ out_file_name+ ".txt";
        std::cout << number_of_points << " points succesfully chosen randomly from "<< dimension <<"input points " << std::endl;
    }
    else {
        std::cerr << "Wrong choice for input manifold..." <<std::endl;  
        exit(-1); 
    }
    std::cout << "Point Set Generated."  <<std::endl;
   
    for(int i = 0; i < number_of_points; i++ )
        point_generator.print_point(point_vector->at(i));
        
    double threshold =  begin_thresold;
    int i = 0;

    while(threshold <= end_thresold)
    {
	 	subComplex = new Fake_simplex_tree();
        Rips_complex rips_complex_from_points(*point_vector, threshold, Gudhi::Euclidean_distance());
        rips_complex_from_points.create_complex(*subComplex, dime);
        std::cout << "Rips complex for parameter " << threshold << " created." <<std::endl;
        
        mat_coll = new SparseMsMatrix(*subComplex);
        
        origTowerFormater.print_tower_for_two_cmplxs(*mat_prev, *mat_coll, map_empty, threshold, origFile);
        std::cout << "Original Tower updated" << std::endl;
        delete mat_prev;
        mat_prev = new SparseMsMatrix(*subComplex);

        totCollapseTime += mat_coll->strong_collapse();
        std::cout << "Subcomplex #" << (i+1) << " Collapsed" << std::endl;
        
        Map redmap = mat_coll->reduction_map();     
        
        std::cout << "Uncollapsed Rips complex is of dimension " << mat_coll->initial_dimension() << " with " << subComplex->num_simplices() << " maximal simplices " << std::endl;
        std::cout << "Collapsed Rips complex is of dimension " << mat_coll->collapsed_dimension() << " with " <<  mat_coll->number_max_simplices() << " maximal simplices" << std::endl;

        std::cout << "** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** " << std::endl;
        
        totPrintTime += towerFormater.print_tower_for_two_cmplxs(*mat_prev_coll, *mat_coll, redmap, threshold, collFile);
        std::cout << "Tower updated for subcomplex #" << i+1 << std::endl; 
        
        delete mat_prev_coll;
        delete subComplex;
        mat_prev_coll = mat_coll; 
        threshold = threshold+steps;
        i++;
    }
    
    std::ofstream myfile (collFile, std::ios::app);
    if(myfile.is_open()){
        myfile << "# Total time taken in all collapses is: " << totCollapseTime << " ms" << std::endl;
        myfile  << "# Total time taken to print the tower format: " <<   totPrintTime  <<" ms" << std::endl;
        myfile.close();
    }
    else {
        std::cerr << "Unable to open file";
        exit(-1) ;
    }

    std::cout << "Total time taken in all collapses is: " << totCollapseTime << " ms." <<std::endl;
    std::cout << "Total time taken to print the tower format: " <<   totPrintTime << " ms." <<std::endl;
    delete mat_prev_coll;
    delete mat_prev;
    return 0;
}
