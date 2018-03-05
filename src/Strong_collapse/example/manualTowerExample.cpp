
#include <gudhi/FormatTower.h>
// #include <gudhi/VertexMaximalMaps.h>
#include <gudhi/Fake_simplex_tree.h>


using Fake_simplex_tree  = Gudhi::Fake_simplex_tree ;

int main()
{
    Map map_empty;
		
    std::size_t tot_num_points = 7;

	std::vector<std::vector<std::size_t>> contents;
    std::vector<SparseMsMatrix> filtered_complex;
    std::vector<SparseMsMatrix> filtered_complex_collapsed;

    Fake_simplex_tree * stree       = new Fake_simplex_tree();
    SparseMsMatrix * mat_coll       = new SparseMsMatrix(); 
    SparseMsMatrix * mat_prev_coll  = new SparseMsMatrix(tot_num_points, 100*tot_num_points); 

    FormatTower towerFormater(tot_num_points) ;
    FormatTower origTowerFormater(tot_num_points);

    filtered_complex.emplace_back(SparseMsMatrix(tot_num_points, 100*tot_num_points));

    std::cout << "Please enter the #simplicial complexes and then the complexes in the following format" <<std::endl;
    std::cout << "Please enter the #vertices and then the vertices of the maximal simplices in the following format" <<std::endl;
    std::cout << "3" << std::endl;
    std::cout << "1 2 3" << std::endl;
    std::cout << "Enter -1 after each complex" << std::endl;
    std::cout << "Enter s to end the input" << std::endl;

    int number, complex_nums;
    
    while(std::cin >> complex_nums)  {
        for(int i = 0; i < complex_nums; i++) {
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
            for (auto & row : contents)
                stree->insert_simplex_and_subfaces(row);

            std::cout << "Updated the simplex tree" << std::endl;
            filtered_complex.emplace_back(SparseMsMatrix(*stree));
            filtered_complex_collapsed.emplace_back(SparseMsMatrix(*stree));
            contents.clear();
        }
    }
    delete stree;

    for(int i = 0 ; i < complex_nums ; ++i) {       
        
        // if( i == complex_nums-1)
        origTowerFormater.print_tower_for_two_cmplxs(filtered_complex.at(i),  filtered_complex.at(i+1), map_empty, i, "./PersistenceOutput/OriginalTowerRips_manual.txt");
        
        mat_coll = & filtered_complex_collapsed.at(i);
        mat_coll->strong_collapse();
        std::cout << "Subcomplex #" << (i+1) << " Collapsed" << std::endl;
        
        Map redmap = mat_coll->reduction_map();  
            
        towerFormater.print_tower_for_two_cmplxs(*mat_prev_coll, *mat_coll, redmap, i, "./PersistenceOutput/CollapsedTowerRips_manual.txt");
        std::cout << "Tower updated for subcomplex #" << i+1 << std::endl; 
        // delete mat_prev_coll;
        // mat_prev_coll = new SparseMsMatrix();
        mat_prev_coll = mat_coll;
        mat_coll  = new SparseMsMatrix();
    }
    // delete mat_coll;
    return 0;
}


