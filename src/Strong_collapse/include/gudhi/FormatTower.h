#pragma once
#include <gudhi/SparseMsMatrix.h>

#include <set>
#include <fstream>
#include <string>
#include <algorithm>


using Fake_simplex_tree     = Gudhi::Fake_simplex_tree;
using Vertex                = Fake_simplex_tree::Vertex;
using Simplex               = Fake_simplex_tree::Simplex;

using vectorVertex          = std::vector<Vertex>;
using vert_unSet            = std::unordered_set<Vertex>;
using simplexVector         = std::vector<Simplex>;


// assumptions : (1) K1 and K2 have the same vertex set
//               (2) The set of simplices of K1 is a subset of set of simplices of K2
// K1  ->  K2    [Original Simplicial Complexes]
// |       |
// |       |
// K1c ->  K2c   [Strongly Collapsed Simplicial Complexes]

class TowerAssembler
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

	template <typename Input_vertex_range>
	std::vector<Simplex> all_faces(const Input_vertex_range &vertex_range){
	    int set_size = vertex_range.size();
	    unsigned int pow_set_size = pow(2, set_size);
	    unsigned int counter, j;
	    std::vector<Simplex> facets;
	    std::vector<Vertex> maxSimplex(vertex_range.begin(), vertex_range.end());
	    // std::sort(maxSimplex.begin(), maxSimplex.end(), vertex_compare);

	    /*Run from counter 000..0 to 111..1*/
	    for(counter = 1; counter < pow_set_size; counter++)
	    {
	      Simplex f;
	      for(j = 0; j < set_size; j++)
	       {          
	          if(counter & (1<<j))                    /* Check if jth bit in the counter is set/true If set then inserts jth element from vertex_range */
	            f.insert(maxSimplex[j]);
	       }
	       facets.emplace_back(f);
	       f.clear();
	    }
	    return facets;
	}
    template< typename Input_vertex_range>
    std::vector<Vertex> sort(const Input_vertex_range & vertex_range) {
        std::vector<Vertex> soreted_simplex(vertex_range.begin(), vertex_range.end());
        std::sort(soreted_simplex.begin(), soreted_simplex.end(), vertex_compare);
        return soreted_simplex;
    }


  public:
    
    TowerAssembler(std::size_t numVert)
    {
    	for (std::size_t i = 0; i < numVert; ++i){
    		renamedVertices[i] = i;
    	}
    	current_rename_counter = numVert;
    }
    
    ~TowerAssembler(){};

    double print_tower_for_two_cmplxs(SparseMsMatrix mat_1, const SparseMsMatrix & mat_2,  Map redmap_2,  double filtration_value, std::string outFile) // mat_1 and mat_2 are simplex_trees of K1c and K2c (the collapsed ones), redmap_2 is the map of K2 -> K2c
    {
        auto begin_print  = std::chrono::high_resolution_clock::now();
        std::ofstream myfile (outFile, std::ios::app);
        if (myfile.is_open())
        {   
            for (auto & v : mat_1.vertex_set()) {
                auto collapsed_to = redmap_2.find(v); 
                if(collapsed_to != redmap_2.end()) {
                    if(mat_1.membership(collapsed_to->second)) {
                    	myfile << filtration_value  << " c " << renamedVertices.at(v) << " " << renamedVertices.at(collapsed_to->second) << std::endl; 
                    	// std::cout << filtration_value << " c " << renamedVertices.at(v) << " " << renamedVertices.at(collapsed_to->second) << std::endl;
                    	renamedVertices.at(v) = current_rename_counter;
                    	current_rename_counter++;
                    }
                    else {
	                    myfile << filtration_value << " i " << renamedVertices.at(collapsed_to->second) << std::endl;
	                    myfile  << filtration_value << " c " << renamedVertices.at(v) << " " << renamedVertices.at(collapsed_to->second) << std::endl; 
	                    // std::cout << filtration_value << " i " << renamedVertices.at(collapsed_to->second) << std::endl;
	                    // std::cout  << filtration_value << " c " << renamedVertices.at(v) << " " << renamedVertices.at(collapsed_to->second) << std::endl; 
	                    renamedVertices.at(v)= current_rename_counter;
	                    current_rename_counter++;
                    }
                    mat_1.contraction(v, collapsed_to->second);  // If the vertex "collapsed_to->second" is not a member of mat_1, the contraction function will simply add and then collapse
                }
            }

            //The core K1c (mat_1) has gone through the transformation(re-labeling)/collapse and it is now a subcomplex of K2c, the remaining simplices need to be included
            // Writing the inclusion of all remaining simplices...
            std::vector<std::size_t>  renamed_simplex;
            for( const Simplex & m  : mat_2.max_simplices()) {
                if(!mat_1.membership(m)) {
                    for(const Simplex & s: all_faces(m)) {
                        if(!mat_1.membership(s)) {
                            for(auto & v : s) {
                               renamed_simplex.push_back(renamedVertices.at(v)); 
                            }
                            myfile << filtration_value << " i";
                            // std::cout << filtration_value << " i";
                            for (Vertex v : sort(renamed_simplex)) {
                                myfile  << " " << v;
                                // std::cout << " " << v;
                            }
                            // std::cout << std::endl;
                            myfile  << std::endl;
                            renamed_simplex.clear();
                        }
         			}
                    mat_1.insert_maximal_simplex_and_subfaces(m);
                }                        
            }   

            myfile << "# Tower updated for the additional new complex.\n";
            myfile.close();
        }
        else {
            std::cerr << "Unable to open file";
            exit(-1) ;
        }
        
        auto end_print  = std::chrono::high_resolution_clock::now();
        auto printTime = std::chrono::duration<double, std::milli>(end_print- begin_print).count();
        std::cout << " Time to print the tower : " << printTime << " ms\n" << std::endl;
        return printTime;
    }
};
