#include <gudhi/FormatTower.h>

#include <gudhi/Rips_complex.h>
#include <gudhi/distance_functions.h>
#include <gudhi/reader_utils.h>

#include <gudhi/PointSetGen.h>
//#include "gudhi/Simplex_tree.h"

using Fake_simplex_tree   = Gudhi::Fake_simplex_tree ;
using Simplex_tree        = Gudhi::Simplex_tree<>;
using Filtration_value    = Fake_simplex_tree::Filtration_value;
using Rips_complex        = Gudhi::rips_complex::Rips_complex<Filtration_value>;
using Distance_matrix     = std::vector<std::vector<Filtration_value>>;

using Vector_of_points = std::vector<Point>;


template<class edge_list, class edge_filt_list >
size_t extract_sub_rips_complex(double threshold, Fake_simplex_tree & currentComplex, edge_list & edges, edge_filt_list & edge_f, size_t end_idx) {
    size_t idx = find_index(edge_f, threshold, 0, end_idx);
    for( idx; idx <= edge_f.size_t(); idx++)
        currentComplex.filtration_simplex_range();
    return idx;    
}

template<class edge_filt >
double find_index(const edge_filt & edge_f, double threshold, size_t begin_idx, size_t end_idx) {
size_t mid_idx = (end_idx - begin_idx)/2;
if( mid_idx <= edge_f.size()){
    if(edge_f.at(mid_idx -1) < threshold and edge_f.at(mid_idx) >= threshold)
        return mid_idx;
    if( edge_f.at(mid_idx) < threshold) 
       find_index(edge_f, threshold, mid_idx, end_idx);
    else
       find_index(edge_f, threshold, begin_idx, mid_idx);
    }
else 
 return -1;    
}

int main(int argc, char * const argv[]) {
	
    PointSetGen point_generator;
    std::string out_file_name   = "default";
    std::string in_file_name    = "default";
    std::size_t number_of_points;
    typedef int Vertex_handle;
    typedef std::vector< std::pair< Vertex_handle, Vertex_handle > > edge_list;
    typedef std::vector< Filtration_value > edge_filt_list;
    typedef typename boost::adjacency_list < boost::vecS, boost::vecS, boost::undirectedS
      , boost::property < vertex_filtration_t, Filtration_value >
      , boost::property < edge_filtration_t, Filtration_value >> OneSkeletonGraph;
    
    edge_list edges;
    edge_filt_list edge_f;

    int     dimension;
    double  begin_thresold;
    double  end_thresold;
    double  steps;
    int     repetetions = 1;
    char    manifold;

    Vector_of_points * point_vector;
    Vector_of_points file_all_points;
    
    int dim_max = 20000; // pseudo variable... of no use
    std::string manifold_full = "sphere";
    
    double radius  = 1;
    double r_min = 0.6;
    double r_max = 1;

    point_generator.program_options(argc, argv, number_of_points, begin_thresold, steps, end_thresold, repetetions, manifold, dimension, in_file_name, out_file_name);
    
    std::cout << "The current input values to run the program is: "<< std::endl;
    std::cout << "number_of_points, begin_thresold, steps, end_thresold, repetetions, manifold, dimension, in_file_name, out_file_name" << std::endl;
    std::cout << number_of_points << ", " << begin_thresold << ", " << steps << ", " << end_thresold << ", " << repetetions << ", " << manifold << ", " << dimension << ", " << in_file_name << ", " << out_file_name << std::endl;
    
    if(manifold == 'f' || manifold =='F') {
        Gudhi::Points_off_reader<Point> off_reader(in_file_name);
        if (!off_reader.is_valid()) {
            std::cerr << "Unable to read file " << in_file_name << "\n";
            exit(-1);  // ----- >>
        }

        file_all_points = Vector_of_points(off_reader.get_point_cloud());
        dimension = file_all_points[0].dimension() ;
        std::cout << "Successfully read " << file_all_points.size() << " point_vector.\n";
        std::cout << "Ambient dimension is " << dimension << ".\n";
    }
   
    Map map_empty;

    std::string origFile ("./PersistenceOutput/original_tower_rips" );
    std::string collFile  ("./PersistenceOutput/collapsed_tower_rips") ;
    
    std::string otherStats ("./PersistenceOutput/maximal_simplx_cnt");
    otherStats = otherStats+"_"+ out_file_name+ ".txt";

    double currentCollapseTime = 0.0;
    double totCollapseTime = 0.0;
    double maxCollapseTime = 0.0;
    double totPrintTime = 0.0;
    double maxCreationAndcollapseTime = 0.0;
   
    point_vector = new Vector_of_points();
    Distance_matrix distances;

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
    
    else if(manifold == 'f' || manifold =='f') {
        // Subsampling from all points for each iterations
        Gudhi::subsampling::pick_n_random_points(file_all_points, number_of_points, std::back_inserter(*point_vector));
        origFile = origFile+"_"+ out_file_name+ ".txt";
        collFile = collFile+"_"+ out_file_name+ ".txt";
        std::cout << number_of_points << " points succesfully chosen randomly of dimension "<< dimension << " ." << std::endl;
    }
    else if (manifold == 'm'){
        std::string csv_file_name(in_file_name);
        distances = read_lower_triangular_matrix_from_csv_file<Filtration_value>(csv_file_name);
        number_of_points = distances.size();
        std::cout << "Read the distance matrix succesfully, of size: " << number_of_points << std::endl;
        origFile = origFile+"_"+ out_file_name+ ".txt";
        collFile = collFile+"_"+ out_file_name+ ".txt";
    }
    else {
        std::cerr << "Wrong parameters for input manifold..." <<std::endl;  
        exit(-1); 
    }
    std::cout << "Point Set Generated."  <<std::endl;
    std::ofstream statsfile (otherStats, std::ios::app);
    if(statsfile.is_open()){
        statsfile << " #number_of_points, begin_thresold, steps, end_thresold, repetetions, manifold, dimension, in_file_name, out_file_name" << std::endl;
        statsfile << " Original_maximal_simplex, Original_complex_dimension, Collapsed_maximal_simplex, Collapsed_complex_dimension" << std::endl;
    }
    else {
        std::cerr << "Unable to open stats file";
        exit(-1) ;
    }
    // for(int i = 0; i < number_of_points; i++ )
    //     point_generator.print_point(point_vector->at(i));

    Fake_simplex_tree * subComplex ; 
    SparseMsMatrix * mat_coll       = new SparseMsMatrix(); 
    SparseMsMatrix * mat_prev_coll  = new SparseMsMatrix(number_of_points, 100*number_of_points); 

    FormatTower towerFormater(number_of_points) ;

    auto threshold =  begin_thresold;
    int iterations = (end_thresold - begin_thresold)/steps;
    int i = 1;
    std::cout << "Total number of iterations to be run are: " << iterations << std::endl;
    
    while(threshold <= end_thresold) {
	 	auto begin1 = std::chrono::high_resolution_clock::now();
        subComplex = new Fake_simplex_tree();
    
        if(manifold == 'm'){
            Rips_complex rips_complex_from_file(distances, threshold);
            rips_complex_from_file.create_complex(*subComplex, dim_max);
        }
        else{
            Rips_complex rips_complex_from_points(*point_vector, threshold, Gudhi::Euclidean_distance());
            rips_complex_from_points.create_complex(*subComplex, dim_max);
        }
        auto inter1 = std::chrono::high_resolution_clock::now();
        mat_coll = new SparseMsMatrix(*subComplex);
        
        auto end1 = std::chrono::high_resolution_clock::now();
        auto matrixForm = std::chrono::duration<double, std::milli>(end1- inter1).count();
        auto creationTime = std::chrono::duration<double, std::milli>(end1- begin1).count();

        std::cout << " SparseMsMatrix using ToplexMap is created in time: (" << creationTime << ", " << matrixForm <<") ms.\n";
        
        currentCollapseTime = mat_coll->strong_collapse();
        std::cout << " Subcomplex #" << i << " Collapsed" << std::endl;
        
        if( maxCreationAndcollapseTime < (currentCollapseTime + creationTime))
        	maxCreationAndcollapseTime = currentCollapseTime + creationTime;
        
        totCollapseTime += currentCollapseTime ;
        if( maxCollapseTime < currentCollapseTime)
            maxCollapseTime = currentCollapseTime;
        
        
        Map redmap = mat_coll->reduction_map();     
        std::cout << "Uncollapsed Rips complex is of dimension " << mat_coll->initial_dimension() << " with " << subComplex->num_simplices() << " maximal simplices " << std::endl;
        std::cout << "Collapsed Rips complex is of dimension " << mat_coll->collapsed_dimension() << " with " <<  mat_coll->number_max_simplices() << " maximal simplices" << std::endl;

        if(statsfile.is_open()){
            statsfile << subComplex->num_simplices() << "," << mat_coll->initial_dimension() << "," << mat_coll->number_max_simplices() << "," << mat_coll->collapsed_dimension() << std::endl;
        }
        else {
            std::cerr << "Unable to open stats file";
            exit(-1) ;
        }
        
        totPrintTime += towerFormater.print_tower_for_two_cmplxs(*mat_prev_coll, *mat_coll, redmap, threshold, collFile);
        std::cout << "Tower updated for subcomplex #" << i << std::endl; 
        std::cout << "** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** " << std::endl;
        delete subComplex;
        mat_prev_coll = mat_coll; 
        threshold = threshold+steps;
        i++;
    }

    statsfile.close();
    
    std::ofstream myfile (collFile, std::ios::app);
    if(myfile.is_open()){
        myfile << "# The input parameter to run the experiment are: "<< std::endl;
        myfile << "# number_of_points, begin_thresold, steps, end_thresold, repetetions, manifold, dimension, in_file_name, out_file_name" << std::endl;
        myfile << "# "<< number_of_points << ", " << begin_thresold << ", " << steps << ", " << end_thresold << ", " << repetetions << ", " << manifold << ", " << dimension << ", " << in_file_name << ", " << out_file_name << std::endl;
        myfile << "# Maximum time taken of all collapses is: " << maxCollapseTime << " ms" << std::endl;
        myfile << "# Total time taken in all collapses is: " << totCollapseTime << " ms" << std::endl;
        myfile  << "# Total time taken to print the tower format: " <<   totPrintTime  <<" ms" << std::endl;
        myfile.close();
    }
    else {
        std::cerr << "Unable to open file";
        exit(-1) ;
    }
    
    std::cout << "Maximum of (Creation + Collapse) TIME : " << maxCreationAndcollapseTime << " ms.\n";
    std::cout << "# Maximum time taken of all collapses is: " << maxCollapseTime << " ms" << std::endl;
    std::cout << "Total time taken in all collapses is: " << totCollapseTime << " ms." <<std::endl;
    std::cout << "Total time taken to print the tower format: " <<   totPrintTime << " ms." <<std::endl;
    std::cout << "** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** " << std::endl;
    delete mat_prev_coll;
    //delete mat_prev;
    return 0;
}
