#include <gudhi/SparseMsMatrix.h>
#include <gudhi/Fake_simplex_tree.h>
#include <set>
#include <fstream>
#include <string>

using Fake_simplex_tree     = Gudhi::Fake_simplex_tree;
using Vertex                = Fake_simplex_tree::Vertex;
using Simplex               = Fake_simplex_tree::Simplex;

using Vert_Set              = std::set<Vertex>;
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
// Map renamedVertices; 
// std::size_t renameCounter;

Vert_Set complex_vertex_range(const simplexVector& maxSimplices)
{
    Vert_Set set;
    set.clear();
    for(const Simplex& s : maxSimplices)
        for (Vertex v : s)
            set.emplace(v);
    return set;  
}
bool vertex_compare (std::size_t  i, std::size_t  j) { return (i<j); }

template <typename Input_vertex_range>
std::vector<Simplex> all_faces(const Input_vertex_range &vertex_range){
    int set_size = vertex_range.size();
    unsigned int pow_set_size = pow(2, set_size);
    unsigned int counter, j;
    std::vector<Simplex> facets;
    std::vector<Vertex> maxSimplex(vertex_range.begin(), vertex_range.end());

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
    
    FormatTower() //std::size_t numVert)
    {
    	// for (std::size_t i = 0; i < numVert; ++i)
    	// {
    	// 	renamedVertices[i] = i;
    	// }

    	// renameCounter = numVert;
    };
    
    ~FormatTower(){};

    void one_step_tower(SparseMsMatrix st1, const SparseMsMatrix &st2,  Map collmap, std::string outFile) // st1 and st2 are simplex_trees of K1c and K2c (the collapsed ones), collmap is the map of K2 -> K2c
    {

        simplexVector maxSimplices1, maxSimplices2, contractedEdges;
        Vert_Set set1; 
        Vert_Set set2;
        vert_unSet contracted_edge;
        int ifcount     = 0;
        int elsecount   = 0 ;
        std::ofstream myfile (outFile, std::ios::app);
        if (myfile.is_open())
        {

            for(auto &v : st1.vertex_set())
                set1.emplace(v); 
            for(auto &v : st2.vertex_set())
                set2.emplace(v);
            // std::size_t set1_size = set1.size();
            // std::size_t set1_size = set2.size();
            set1.sort(vertex_compare);

            // std::sort (set1.begin(), set1.end(),this->vertex_compare());
            // std::sort (set2.begin(), set2.end(),this->vertex_compare());

            Vert_Set diff;
            diff.clear();
            std::set_difference(set1.begin(), set1.end(), set2.begin(), set2.end(), std::inserter(diff, diff.end())); // It computes the set (V1c - V2c) // These vertices are the ones which go through collapse. 
            
            																										  // We write them one by one and transform the st1 (simplicial complex K1) accordingly, which will finally be a subcomplex of K2c.			
            
            for (Vert_Set::iterator iter=diff.begin(); iter!=diff.end(); iter++)
            {
                auto collapsed_to = collmap[*iter];
                Simplex s_iter;
                s_iter.insert(*iter);
                Simplex s_ct;
                s_ct.insert(collapsed_to);
                if(set1.count(collapsed_to)!=0)
                {
                    st1.contraction(*iter, collapsed_to);  // If the vertex collapsed_to is not a vertex of st1, the contraction function will simply add               
                    contracted_edge.insert(*iter);
                    contracted_edge.insert(collapsed_to);
                    ifcount++;
                    contractedEdges.emplace_back(contracted_edge);
                    contracted_edge.clear();
                    std::cout<< "Entered in to the IF control block! Motive: contraction from "<<*iter << " to " << collapsed_to << std::endl;
                    if(st1.membership(s_iter))
                        std::cout << "Error: the vertex "<< *iter << " has been contracted and it's still a member!! :(" << std::endl; 
                    if(!st1.membership(s_iter))
                        std::cout << "The vertex "<< *iter << " has been contracted and it's not a member anymore :) " << std::endl; 

                    if(!st1.membership(s_ct))
                        std::cout << "Error: the vertex : " <<collapsed_to<< " is not a member! :(" <<std::endl;
                    if(st1.membership(s_ct))
                        std::cout << "The vertex : " <<collapsed_to<< " is still a member :)" <<std::endl;

                    myfile  << "c " << *iter << " " << collapsed_to << std::endl; 

                    // myfile  << "c " << renamedVertices[*iter] << " " << renamedVertices[collapsed_to] << std::endl; 
                    // renamedVertices[*iter] = renameCounter;
                    // renameCounter++;
                }
                else 
                {
                    std::cout<< "Else block Motive: contraction from "<<*iter << " to " << collapsed_to << std::endl;

                    st1.insert_simplex_and_subfaces(s_ct);
                    
                    if(!st1.membership(s_ct))
                        std::cout << "Error in Else: In insertion, the vertex " << collapsed_to << " should be a member!! :( " <<std::endl;

                    myfile  << "i " << collapsed_to << std::endl;  // This is a non-dominated vertex so it will never be collapsed in this iteration.

                    // myfile  << "i " << renamedVertices[collapsed_to] << std::endl;  // This is a non-dominated vertex so it will never collapse in this iteration.
                    
                    st1.contraction(*iter, collapsed_to);
                    contracted_edge.insert(*iter);
                    contracted_edge.insert(collapsed_to);
                    elsecount++;
                    contractedEdges.emplace_back(contracted_edge);
                    contracted_edge.clear();
                    
                    if(st1.membership(s_iter))
                        std::cout << "Error in Else: the vertex "<< *iter << " has been contracted and it's still a member!! :( " << std::endl; 
                    if(!st1.membership(s_iter))
                        std::cout << "In Else: The vertex "<< *iter << " has been contracted and it's not a member anymore :) " << std::endl; 

                    if(!st1.membership(s_ct))
                        std::cout << "Error in Else: the vertex : " <<collapsed_to<< " is not a member :(" <<std::endl;
                    if(st1.membership(s_ct))
                        std::cout << "In Else, the vertex : " <<collapsed_to<< " is still a member :)" <<std::endl;


                    myfile  << "c " << *iter << " " << collapsed_to << std::endl; 
      
                    // myfile  << "c " << renamedVertices[*iter] << " " << renamedVertices[collapsed_to] << std::endl; 
                    // renamedVertices[*iter] = renameCounter;
                    // renameCounter++;

                }
                
            }                                            

            //The core K1c (st1) has gone through the transformation(re-labeling)/collapse and it is now a subcomplex of K2c, the remaining simplices need to be included
            // Writing the inclusion of all remaining simplices...
            Vert_Set invertedV;
            Vert_Set::iterator iter;

            for (const Simplex& m : contractedEdges)
            {   std::cout << "Current dimension of the maximal simplex is :" <<(m.size()-1)<<std::endl;
                if(st1.membership(m))
                    std::cout << "Bug big bug!!! :( " <<std::endl;
            }
            std::cout << "If count is: " << ifcount << " Else count is: " << elsecount << std::endl;
            for( const Simplex& m  : maxSimplices2 )
            {
                for(const Simplex& s: all_faces(m))
                {
                    if(!st1.membership(s))
                    {
                        myfile  << "i";
                        for (Vertex v : s)
                        {
                            invertedV.insert(v);
                        }
                        for(iter = invertedV.begin(); iter != invertedV.end(); iter++)
                        {
                            myfile  << " " << *iter;
                        }
                        myfile  << std::endl;
                        invertedV.clear();
                    }
                }
                st1.insert_simplex_and_subfaces(m);        
            }   

            myfile << "# This is one step tower.\n";
            myfile.close();
            // std::cout << "The rename-counter is currently set to: " << renameCounter << std::endl;
        }
        else
        {
            std::cout << "Unable to open file";
            return ;
        }
        return;
    }

    
};
