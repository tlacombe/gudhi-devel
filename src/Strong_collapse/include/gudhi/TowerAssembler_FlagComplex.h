#pragma once
#include <gudhi/FlagComplexSpMatrix.h>
#include <gudhi/Rips_edge_list.h>

#include <set>
#include <fstream>
#include <string>
#include <algorithm>


typedef std::size_t Vertex;
using Edge                  = std:pair<Vertex,Vertex>;
using Edge_list             = std:vector<Edge>;
using Simplex               = std::vector<Vertex>;

using vectorVertex          = std::vector<Vertex>;
using vert_unSet            = std::unordered_set<Vertex>;


// assumptions : (1) K1 and K2 have the same vertex set
//               (2) The set of simplices of K1 is a subset of set of simplices of K2
// K1  ->  K2    [Original Simplicial Complexes]
// |       |
// |       |
// K1c ->  K2c   [Strongly Collapsed Flag Complexes]

class TowerAssembler_FlagComplex
{
  private:
	Map renamedVertices; 
	std::size_t current_rename_counter;
	
    struct {
	        bool operator()(std::size_t a, std::size_t b) const
	        {   
	            return a < b;
	        }   
	    } vertex_compare;

    typedef std::vector< std::tuple<Filtration_value, Vertex_handle, Vertex_handle > > Filtered_sorted_edge_list;

	// template <typename Input_vertex_range>
	// std::vector<Simplex> all_faces(const Input_vertex_range &vertex_range){
	//     int set_size = vertex_range.size();
	//     unsigned int pow_set_size = pow(2, set_size);
	//     unsigned int counter, j;
	//     std::vector<Simplex> facets;
	//     std::vector<Vertex> maxSimplex(vertex_range.begin(), vertex_range.end());
	//     // std::sort(maxSimplex.begin(), maxSimplex.end(), vertex_compare);

	//     /*Run from counter 000..0 to 111..1*/
	//     for(counter = 1; counter < pow_set_size; counter++)
	//     {
	//       Simplex f;
	//       for(j = 0; j < set_size; j++)
	//        {          
	//           if(counter & (1<<j))                    /* Check if jth bit in the counter is set/true If set then inserts jth element from vertex_range */
	//             f.insert(maxSimplex[j]);
	//        }
	//        facets.emplace_back(f);
	//        f.clear();
	//     }
	//     return facets;
	// }
 //    template< typename Input_vertex_range>
 //    std::vector<Vertex> sort(const Input_vertex_range & vertex_range) {
 //        std::vector<Vertex> soreted_simplex(vertex_range.begin(), vertex_range.end());
 //        std::sort(soreted_simplex.begin(), soreted_simplex.end(), vertex_compare);
 //        return soreted_simplex;
 //    }


  public:
    
    TowerAssembler_FlagComplex(std::size_t numVert)
    {
    	for (std::size_t i = 0; i < numVert; ++i){
    		renamedVertices[i] = i;
    	}
    	current_rename_counter = numVert;
        Filtered_sorted_edge_list * edge_t = new Filtered_sorted_edge_list();
        FlagComplexSpMatrix * flag_Filtration = new FlagComplexSpMatrix();


    }
    
    ~TowerAssembler_FlagComplex(){};
    double print_tower_for_two_cmplxs(FlagComplexSpMatrix mat_1, const FlagComplexSpMatrix & mat_2,  Map redmap_2,  double filtration_value, std::string outFile) // mat_1 and mat_2 are simplex_trees of K1c and K2c (the collapsed ones), redmap_2 is the map of K2 -> K2c
    {
        auto begin_print  = std::chrono::high_resolution_clock::now();
        std::ofstream myfile (outFile, std::ios::app);
        if (myfile.is_open())
        {   
            for (auto & v : mat_1.vertex_set()) {
                auto collapsed_to = redmap_2.find(v); 
                if(collapsed_to != redmap_2.end()) {  // Collapse happend, there is an existing vertex in the map
                    if(mat_1.membership(collapsed_to->second)) { // Collapsed to an existing vertex.

                    	myfile << filtration_value  << " c " << renamedVertices.at(v) << " " << renamedVertices.at(collapsed_to->second) << std::endl; 
                    	// std::cout << filtration_value << " c " << renamedVertices.at(v) << " " << renamedVertices.at(collapsed_to->second) << std::endl;
                        auto contracted = flag_Filtration->active_strong_expansion(renamedVertices.at(v), renamedVertices.at(collapsed_to->second));
                        renamedVertices.at(contracted) = current_rename_counter;
                    	current_rename_counter++;
                    }
                    else {
	                    myfile << filtration_value << " i " << renamedVertices.at(collapsed_to->second) << std::endl;
	                    myfile  << filtration_value << " c " << renamedVertices.at(v) << " " << renamedVertices.at(collapsed_to->second) << std::endl; 
	                    // std::cout << filtration_value << " i " << renamedVertices.at(collapsed_to->second) << std::endl;
	                    // std::cout  << filtration_value << " c " << renamedVertices.at(v) << " " << renamedVertices.at(collapsed_to->second) << std::endl; 
	                    flag_Filtration->insert(renamedVertices.at(collapsed_to->second));
                        auto contracted = flag_Filtration->active_strong_expansion(renamedVertices.at(v), renamedVertices.at(collapsed_to->second));
                        renamedVertices.at(contracted)= current_rename_counter;
	                    current_rename_counter++;
                    }
                    mat_1.contraction(v, collapsed_to->second);  // If the vertex "collapsed_to->second" is not a member of mat_1, the contraction function will simply add and then collapse
                }
            }

            //The core K1c (mat_1) has gone through the transformation(re-labeling)/collapse and it is now a subcomplex of K2c, the remaining simplices need to be included
            // Writing the inclusion of all remaining simplices...
            std::vector<std::size_t>  renamed_simplex;
            for( const Edge & e  : mat_2.all_edges()) {
                auto u = e.begin();
                auto v = e.end();
                if(!mat_1.membership(u)) {
                    flag_Filtration.insert_vertex(renamedVertices.at(u));
                    myfile << filtration_value << " i";
                    myfile  << " " << renamedVertices.at(u);
                    myfile  << std::endl;
                    mat_1.insert_vertex(u);
                }

                if(!mat_1.membership(v)) {
                    flag_Filtration.insert_vertex(renamedVertices.at(v));
                    myfile << filtration_value << " i";
                    myfile  << " " << renamedVertices.at(v);
                    myfile  << std::endl;
                    mat_1.insert_vertex(v);
                }

                if(!mat_1.membership(e)){
                    flag_Filtration.insert_new_edges({renamedVertices.at(u),renamedVertices.at(v)});
                    myfile << filtration_value << " i";
                    myfile  << " " <<  renamedVertices.at(u) << renamedVertices.at(v);
                    myfile  << std::endl;
                    mat_1.insert_new_edges(e);
                
                }         
            }                        
            myfile << "# Tower updated for the additional subcomplex.\n";
            myfile.close();
        }
        else {
            std::cerr << "Unable to open file";
            exit(-1) ;
        }
        
        auto end_print  = std::chrono::high_resolution_clock::now();
        auto printTime = std::chrono::duration<double, std::milli>(end_print- begin_print).count();
        // std::cout << " Time to print the tower : " << printTime << " ms\n" << std::endl;
        return printTime;
    }    
    
};
