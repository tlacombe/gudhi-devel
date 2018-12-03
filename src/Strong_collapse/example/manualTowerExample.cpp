
#include <gudhi/TowerAssembler_FlagComplex.h>
// #include <gudhi/VertexMaximalMaps.h>
// #include <gudhi/Fake_simplex_tree.h>


// using Fake_simplex_tree  = Gudhi::Fake_simplex_tree ;
typedef std::size_t Vertex;
typedef std::vector< std::tuple< double, Vertex, Vertex > > Filtered_sorted_edge_list;
using Distance_matrix       = std::vector<std::vector<double>>;

int main()
{
    Map map_empty;
		
    std::size_t tot_num_points = 10;

	std::vector<std::vector<std::size_t>> contents;
    // std::vector<FlagComplexSpMatrix> filtered_complex;
    std::vector<FlagComplexSpMatrix *> filtered_complex_collapsed;

    // Fake_simplex_tree * stree            = new Fake_simplex_tree();
    Filtered_sorted_edge_list edge_list;

    FlagComplexSpMatrix * mat_coll       = new FlagComplexSpMatrix(); 
    FlagComplexSpMatrix * mat_prev_coll  = new FlagComplexSpMatrix(tot_num_points); 

    TowerAssembler_FlagComplex towerFormater(tot_num_points) ;
    TowerAssembler_FlagComplex origTowerFormater(tot_num_points);

    // filtered_complex.emplace_back(FlagComplexSpMatrix(tot_num_points));

    std::cout << "Please enter the #simplicial complexes and then the skeleton in the following format" <<std::endl;
    std::cout << "Please enter 2 and then the vertices of the edges in the following format" <<std::endl;
    std::cout << "2" << std::endl;
    std::cout << "2 3" << std::endl;
    std::cout << "Enter -1 after each complex" << std::endl;
    std::cout << "Enter s to end the input" << std::endl;

    int number, complex_nums;
    int i = 0;
    while(std::cin >> complex_nums)  {
        for(; i < complex_nums; i++) {
            while(true) {
                std::cin >> number ;
                if (number > 0 ) {    
                    std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n'); // skip eol
                    std::string line;
                    std::getline(std::cin, line);
                    if (std::cin) {
                        contents.emplace_back(number);
                        std::istringstream iss(line);
                        std::copy_n(std::istream_iterator<double>(iss), number, contents.back().begin());
                    }
                    else
                        break;
                }
                else
                    break ;   
            }
            for (auto & row : contents){
                // std::cout << "Edge list updated with the edge: [" << row[0] << ", " << row[1] << "]" << std::endl; 
                edge_list.emplace_back(0,row[0], row[1]);
            }


            std::cout << "Updated the skeleton tree" << std::endl;
            // filtered_complex.emplace_back(FlagComplexSpMatrix(tot_num_points, edge_list));
            filtered_complex_collapsed.emplace_back(new FlagComplexSpMatrix(tot_num_points, edge_list));
            // contents.clear();
            edge_list.clear();
        }
    }
    contents.clear();
    // FlagComplexSpMatrix * mat;      = new FlagComplexSpMatrix(); 
    // for(int i = 0; i < 3; i++){
    //     mat;      = new FlagComplexSpMatrix()

    // }
    // delete edge_list;
    complex_nums = i;
    double filt_val = 0.1;
    std::cout << "Total number of input complexes are: " << complex_nums << std::endl;
    for(int i = 0 ; i < complex_nums ; ++i) {       
        
        std::cout << "Going for collapse and tower assembly" << std::endl;
        // if( i == complex_nums-1)
            // origTowerFormater.build_tower_for_two_cmplxs(filtered_complex.at(i),  filtered_complex.at(i+1), map_empty, i, "./PersistenceOutput/OriginalTowerRips_manual.txt");
        // std::cout << "Going for collapse and tower assembly" << std::endl;
        mat_coll = filtered_complex_collapsed.at(i);
        // mat_coll->print_sparse_skeleton();
        mat_coll->strong_collapse();
        std::cout << "Subcomplex #" << (i+1) << " Collapsed" << std::endl;
        
        Map redmap = mat_coll->reduction_map();  
            
        towerFormater.build_tower_for_two_cmplxs(*mat_prev_coll, *mat_coll, redmap, filt_val, "./PersistenceOutput/CollapsedTowerRips_manual.txt");
        std::cout << "Tower updated for subcomplex #" << i+1 << std::endl; 
        // delete mat_prev_coll;
        // mat_prev_coll = new FlagComplexSpMatrix();
        mat_prev_coll = mat_coll;
        mat_coll  = new FlagComplexSpMatrix();
        filt_val = filt_val+0.1;
    }
    // delete mat_coll;
    Distance_matrix dist = towerFormater.distance_matrix();
    for(auto & x : dist){
        for( auto & y: x )
            std::cout << y << ", " ;
        std::cout << std::endl;
    }
    towerFormater.print_sparse_matrix();
    return 0;
}


