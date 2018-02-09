#include <gudhi/SparseMsMatrix.h>
// #include <gudhi/Fake_simplex_tree.h>
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

class FormatTower
{
private:
	Map renamedVertices; 
	std::size_t renameCounter;
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
	    std::sort(maxSimplex.begin(), maxSimplex.end(), vertex_compare);
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


public:
    
    FormatTower(std::size_t numVert)
    {
    	for (std::size_t i = 0; i < numVert; ++i)
    	{
    		renamedVertices[i] = i;
    	}

    	renameCounter = numVert;
    };
    
    ~FormatTower(){};

    void one_step_tower(SparseMsMatrix mat_1, const SparseMsMatrix &mat_2,  Map collmap, std::string outFile) // mat_1 and mat_2 are simplex_trees of K1c and K2c (the collapsed ones), collmap is the map of K2 -> K2c
    {
        std::ofstream myfile (outFile, std::ios::app);
        if (myfile.is_open())
        {   
            for (auto & v : mat_1.vertex_set())
            {
                auto collapsed_to = collmap.find(v); 
                if(collapsed_to != collmap.end())
                {
                    if(mat_1.membership(collapsed_to->second))
                    {
                    	myfile  << "c " << renamedVertices[v] << " " << renamedVertices[collapsed_to->second] << std::endl; 
                    	renamedVertices[v] = renameCounter;
                    	renameCounter++;
                    }               

                    else
                    {
	                    myfile << "i " << renamedVertices[collapsed_to->second] << std::endl;
	                    myfile  << "c " << renamedVertices[v] << " " << renamedVertices[collapsed_to->second] << std::endl; 
	                    renamedVertices[v] = renameCounter;
	                    renameCounter++;
                    }
                    mat_1.contraction(v, collapsed_to->second);  // If the vertex collapsed_to->second is not a vertex of mat_1, the contraction function will simply add 
                }             
            }                                            

            //The core K1c (mat_1) has gone through the transformation(re-labeling)/collapse and it is now a subcomplex of K2c, the remaining simplices need to be included
            // Writing the inclusion of all remaining simplices...
            vectorVertex invertedV, renamed_maxSimp;
            vectorVertex::iterator iter;
            for( const Simplex & m  : mat_2.max_simplices() )
            {
                for(auto & v_m : m )
                {
                	renamed_maxSimp.push_back(renamedVertices[v_m]);
                }
                for(const Simplex & s: all_faces(renamed_maxSimp))
                {
                    if(!mat_1.membership(s))
                    {
                        // for (Vertex v : s)
                        // {
                        //     invertedV.push_back(renamedVertices[v]);
                        // }
                        // sort(invertedV.begin(), invertedV.end(), vertex_compare);
                        
                        myfile  << "i";
                        // for(iter = s.begin(); iter != s.end(); iter++)
                        for (Vertex v : s)
                        {
                            myfile  << " " << v;
                        }
                        myfile  << std::endl;
                        // invertedV.clear();
                    }
     			}
     			renamed_maxSimp.clear();
                mat_1.insert_maximal_simplex_and_subfaces(m);    
            }   

            myfile << "# This is one step tower.\n";
            myfile.close();
            // std::cout << "The rename-counter is currently set to: " << renameCounter << std::endl;
        }
        else
        {
            std::cerr << "Unable to open file";
            exit(-1) ;
        }
        return;
    }

    
};
