#include <gudhi/TowerAssembler_FlagComplex.h>
#include <gudhi/Rips_complex.h>
#include <gudhi/Simplex_tree.h>
#include <gudhi/Persistent_cohomology.h>
#include <gudhi/Rips_edge_list.h>
#include <gudhi/distance_functions.h>
#include <gudhi/reader_utils.h>
#include <gudhi/PointSetGen.h>

// Types definition
using Vector_of_points         = std::vector<Point>;
using Vector_of_SM_pointers    = std::vector<FlagComplexSpMatrix*>;

using Simplex_tree = Gudhi::Simplex_tree<Gudhi::Simplex_tree_options_fast_persistence>;
using Filtration_value = double;
using Rips_complex = Gudhi::rips_complex::Rips_complex<Filtration_value>;
using Rips_edge_list      = Gudhi::rips_edge_list::Rips_edge_list<Filtration_value>;
using Field_Zp = Gudhi::persistent_cohomology::Field_Zp;
using Persistent_cohomology = Gudhi::persistent_cohomology::Persistent_cohomology<Simplex_tree, Field_Zp>;
using Distance_matrix = std::vector<std::vector<Filtration_value>>;
class extract_sub_one_skeleton
{
    public:
        template<class Filtered_sorted_edge_list,  class Fil_vector >
        extract_sub_one_skeleton(double threshold, Filtered_sorted_edge_list & current_edge_t, Filtered_sorted_edge_list & edge_t, Fil_vector & edge_filt ) {
            std::cout << "The number of the remaining edges are: " << (edge_t.size()) << std::endl;
            std::cout << "Current number of the edges are: " << (current_edge_t.size()) << std::endl;
            //std::cout << "The longest edge of the current edges is: "<< std::get<0>(*(current_edge_t.end()-1)) << " and threshold is: " << threshold << std::endl;
            auto end_it = std::upper_bound(edge_filt.begin(), edge_filt.end(), threshold); // find_index(edge_t, threshold, 0, end_idx);
            size_t end_idx = std::distance(edge_filt.begin(), end_it);

            for( size_t idx = 0; idx < end_idx ; idx++) {
               current_edge_t.push_back(*edge_t.begin()); //{std::get<0>(),std::get<1>(edge_t.at(0)), std::get<2>(edge_t.at(0))})
               edge_filt.erase(edge_filt.begin());
               edge_t.erase(edge_t.begin());
            }
        }
};


int main(int argc, char * const argv[]) {
	
    // auto the_begin = std::chrono::high_resolution_clock::now();
    PointSetGen point_generator;
    std::string out_file_name   = "default";
    std::string in_file_name    = "default";
    std::size_t number_of_points;
    
    typedef size_t Vertex_handle;
    typedef std::vector< std::tuple<Filtration_value, Vertex_handle, Vertex_handle > > Filtered_sorted_edge_list;
    std::vector<Filtration_value> * edge_filt = new std::vector<Filtration_value>() ;

    int     dimension;
    double  begin_thresold;
    double  end_threshold;
    double  min_dist;
    double  max_dist;
    double  steps;
    int     repetetions = 1;
    char    manifold;

    Vector_of_points * point_vector;
    Vector_of_points file_all_points;

    std::string manifold_full = "sphere";
    
    double radius  = 1;
    double r_min  = 0.6;
    double r_max = 1;
    int dim_max  = 3;

    point_generator.program_options(argc, argv, number_of_points, begin_thresold, steps, end_threshold, repetetions, manifold, dimension, in_file_name, out_file_name);
    
    std::cout << "The current input values to run the program is: "<< std::endl;
    std::cout << "number_of_points, begin_thresold, steps, end_threshold, repetetions, manifold, dimension, in_file_name, out_file_name" << std::endl;
    std::cout << number_of_points << ", " << begin_thresold << ", " << steps << ", " << end_threshold << ", " << repetetions << ", " << manifold << ", " << dimension << ", " << in_file_name << ", " << out_file_name << std::endl;
    
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
    std::string filediag_aft ("./PersistenceOutput/sparse_persistence_diags.txt") ;
    std::string filediag_bfr ("./PersistenceOutput/original_persistence_diags.txt") ;
    
    
    std::string otherStats ("./PersistenceOutput/maximal_simplx_cnt");
    otherStats = otherStats+"_"+ out_file_name+ ".txt";


    double totAssembleTime = 0.0;
    double currentCreationTime = 0.0;
    double maxCreationTime     = 0.0;

   
    point_vector = new Vector_of_points();
    Distance_matrix distances, sparse_distances;

    if(manifold == 's' || manifold == 'S'){
        point_generator.generate_points_sphere(*point_vector, number_of_points, dimension, radius);
        origFile = origFile+"_sphere_"+out_file_name+".txt";
        collFile = collFile+"_sphere_"+out_file_name+".txt";
        std::cout << number_of_points << " points successfully chosen randomly from "<< dimension <<"-sphere of radius " << radius << std::endl;
    }
    else if(manifold == 'b' || manifold == 'B'){
        point_generator.generate_points_ball(*point_vector, number_of_points, dimension, radius); 
        origFile = origFile+"_ball_"+out_file_name+".txt";
        collFile = collFile+"_ball_"+out_file_name+".txt";
        std::cout << number_of_points << " points successfully chosen randomly from "<< dimension <<"-ball of radius " << radius << std::endl;
    
    }
    else if( (manifold == 'a' || manifold == 'A')&& dimension == 2){
        point_generator.generate_points_2annulus(*point_vector, number_of_points, r_min, r_max); 
        origFile = origFile+"_annulus_"+out_file_name+".txt";
        collFile = collFile+"_annulus_"+out_file_name+".txt";
        std::cout << number_of_points << " points successfully chosen randomly from "<< 2 <<"-annulus of radii (" << r_min << ',' << r_max << ") " << std::endl;
    }
    else if( (manifold == 'a' || manifold == 'A') && dimension == 3){
        point_generator.generate_points_spherical_shell(*point_vector, number_of_points, r_min, r_max); 
        origFile = origFile+"_annulus_"+out_file_name+".txt";
        collFile = collFile+"_annulus_"+out_file_name+".txt";
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
    //Preparing the statsfile to record the reduction in num of maximal simplices and the dimension of the complex.
    std::ofstream statsfile (otherStats, std::ios::app);
    if(statsfile.is_open()){
        statsfile << " #number_of_points, begin_thresold, steps, end_threshold, repetetions, manifold, dimension, in_file_name, out_file_name" << std::endl;
        statsfile << " Original_maximal_simplex, Original_complex_dimension, Collapsed_maximal_simplex, Collapsed_complex_dimension" << std::endl;
    }
    else {
        std::cerr << "Unable to open stats file";
        exit(-1) ;
    }

    // for(int i = 0; i < number_of_points; i++ )
    //     point_generator.print_point(point_vector->at(i));

    Filtered_sorted_edge_list * sub_skeleton  = new Filtered_sorted_edge_list();  
    Filtered_sorted_edge_list * edge_t = new Filtered_sorted_edge_list();
    TowerAssembler_FlagComplex twr_assembler(number_of_points) ;
    
    std::cout << "Computing the one-skeleton for threshold: " << end_threshold << std::endl; 
    
    auto begin_full_cmplx = std::chrono::high_resolution_clock::now();
    if(manifold == 'm'){ //Input is a distance 'm'atrix
        //Creating the edge list
        Rips_edge_list Rips_edge_list_from_file(distances, end_threshold);
        Rips_edge_list_from_file.create_edges(*edge_t);
        std::cout<< "Sorted edge list computed" << std::endl;

        //Creating the Rips Complex
        //Rips_complex rips_complex_from_file(distances, end_threshold);
        //rips_complex_from_file.create_complex(*subComplex, dim_max);
        //std::cout<< "Rips complex computed" << std::endl;
    }
    else{ //Point cloud input
         //Creating the edge list
        Rips_edge_list Rips_edge_list_from_points(*point_vector, end_threshold, Gudhi::Euclidean_distance());
        Rips_edge_list_from_points.create_edges(*edge_t);
        std::cout<< "Sorted edge list computed" << std::endl;
        //Creating the Rips Complex
        // Rips_complex rips_complex_from_points(*point_vector, end_threshold, Gudhi::Euclidean_distance());
        // rips_complex_from_points.create_complex(*subComplex, dim_max);
        // std::cout<< "Rips complex computed" << std::endl;
    }

    //An additional vector <edge_filt> to perform binary search to find the index of given threshold
    edge_filt->clear();
    for(auto edIt = edge_t->begin(); edIt != edge_t->end(); edIt++) {
        edge_filt->push_back(std::get<0>(*edIt));
    }
    min_dist = edge_filt->at(0);
    max_dist = edge_filt->at(edge_filt->size()-1);

    if(begin_thresold < min_dist){
        std::cout<< "Begin threshold re-set to the minimum filteration value, " << min_dist << "." <<std::endl;
    }
    int iterations = (end_threshold - min_dist)/steps;
    std::cout << "Total number of iterations to be run are: " << iterations << std::endl;

    auto end_full_cmplx = std::chrono::high_resolution_clock::now();
    currentCreationTime = std::chrono::duration<double, std::milli>(end_full_cmplx - begin_full_cmplx).count();
    maxCreationTime = currentCreationTime;
   
    auto threshold =  min_dist;  

    FlagComplexSpMatrix * mat_coll       = new FlagComplexSpMatrix(); 
    FlagComplexSpMatrix * mat_prev_coll  = new FlagComplexSpMatrix(number_of_points); 

    std::cout << "Going for collapse and tower assembly" << std::endl;

    int i = 1;
    Map * redmap;
    while(threshold <= end_threshold) {
        extract_sub_one_skeleton(threshold, *sub_skeleton, *edge_t, *edge_filt);
       
        mat_coll = new FlagComplexSpMatrix(number_of_points, *sub_skeleton);
        mat_coll->strong_collapse();
        redmap = new Map();
        *redmap = mat_coll->reduction_map(); 
        
        std::cout << "Subcomplex #" << i << " Collapsed" << std::endl;
        totAssembleTime += twr_assembler.build_tower_for_two_cmplxs(*mat_prev_coll, *mat_coll, *redmap, threshold, "./PersistenceOutput/CollapsedTowerRips_manual.txt");
        std::cout << "Tower updated for subcomplex #" << i << std::endl; 
        
        delete mat_prev_coll;
        mat_prev_coll = new FlagComplexSpMatrix();
        mat_prev_coll = mat_coll;
        mat_coll  = new FlagComplexSpMatrix();
        threshold = threshold+steps;
        i++;
        delete redmap;
    }

    sparse_distances = twr_assembler.distance_matrix();
    
    Rips_complex rips_complex_after_collapse(sparse_distances, end_threshold);
    Rips_complex rips_complex_before_collapse(*point_vector, end_threshold, Gudhi::Euclidean_distance());

    // Construct the Rips complex in a Simplex Tree
    
    Simplex_tree simplex_tree_aft, simplex_tree_bfr;
    rips_complex_before_collapse.create_complex(simplex_tree_bfr, dim_max);
    rips_complex_after_collapse.create_complex(simplex_tree_aft, dim_max);

    std::cout << "The complex contains " << simplex_tree_bfr.num_simplices() << " simplices before collapse. \n";
    std::cout << "   and has dimension " << simplex_tree_bfr.dimension() << " \n";

    std::cout << "The complex contains " << simplex_tree_aft.num_simplices() << " simplices  after collapse. \n";
    std::cout << "   and has dimension " << simplex_tree_aft.dimension() << " \n";

    // Sort the simplices in the order of the filtration
    simplex_tree_bfr.initialize_filtration();
    simplex_tree_aft.initialize_filtration();
    // Compute the persistence diagram of the complex
    Persistent_cohomology pcoh_bfr(simplex_tree_bfr);
    Persistent_cohomology pcoh_aft(simplex_tree_aft);
    // initializes the coefficient field for homology
    pcoh_bfr.init_coefficients(2);
    pcoh_aft.init_coefficients(2);

    pcoh_bfr.compute_persistent_cohomology(steps);
    pcoh_aft.compute_persistent_cohomology(0);
    // Output the diagram in filediag
    if (filediag_bfr.empty()) {
        pcoh_bfr.output_diagram();
    } 
    else {
        std::ofstream out(filediag_bfr);
        pcoh_bfr.output_diagram(out);
        out.close();
      }

    if (filediag_aft.empty()) {
        pcoh_aft.output_diagram();
    } 
    else {
        std::ofstream out(filediag_aft);
        pcoh_aft.output_diagram(out);
        out.close();
      }
    // for(auto & x : sparse_distances){
    //     for( auto & y: x )
    //         std::cout << y << ", " ;
    //     std::cout << std::endl;
    // }
    return 0;

}
  