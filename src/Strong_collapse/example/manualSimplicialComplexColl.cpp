#include <gudhi/FlagComplexSpMatrix.h>
// #include <gudhi/VertexMaximalMaps.h>
// #include <gudhi/Fake_simplex_tree.h>

#include <gudhi/Rips_complex.h>
#include <gudhi/distance_functions.h>

#include <CGAL/Epick_d.h>
#include <CGAL/point_generators_d.h>
#include <CGAL/Random.h>
#include <cmath>

#include <iostream>
#include <string>
#include <vector>
#include <limits>  // for std::numeric_limits
#include <sstream>      // for istringstream
#include <algorithm>    // for copy, copy_n
#include <iterator>     // for istream_iterator<>, ostream_iterator<>
#include <set>

using Point = CGAL::Epick_d< CGAL::Dimension_tag<20> >::Point_d;
// using Filtration_value = Fake_simplex_tree::Filtration_value;
// using Rips_complex = Gudhi::rips_complex::Rips_complex<Filtration_value>;

// using Fake_simplex_tree  = Gudhi::Fake_simplex_tree ;
// using Vertex             = Fake_simplex_tree::Vertex;
typedef std::size_t Vertex;
using Vector_of_points = std::vector<Point>;
typedef std::vector< std::tuple< double, Vertex, Vertex > > Filtered_sorted_edge_list;

int main()
{
	
	Vector_of_points points;
	//Vector_of_points noisePoints;
	Filtered_sorted_edge_list edge_list;
	std::vector<std::vector<std::size_t>> contents;
    // std::vector<std::vector<std::size_t>> edges;


    std::cout << "Please enter the vertices of the edge in the following format" <<std::endl;
    std::cout << "2" << std::endl;
    std::cout << "2 3" << std::endl;
    std::cout << "Enter s to stop" << std::endl;
    int number;
    // std::vector<std::vector<std::size_t>> random = { {1,4},{3,4},{4,8}, {1,8},{3,8},{89}, {6}, {2},{4},{5} }; //{{89,77,29},{2,4,6}, {2,4,5},{3,4},{2,3},{6},{2,7},{7,5}, {2,3,4,5}};
    while (std::cin >> number)
    {
        std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n'); // skip eol
        std::string line;
        std::getline(std::cin, line);
        if (std::cin)
        {
            contents.emplace_back(number);
            std::istringstream iss(line);
            std::copy_n(std::istream_iterator<double>(iss), number, contents.back().begin());
        }
        else
        {
            return 255;
        }
    }

    // if (!std::cin.eof())
    //     std::cout << "Warning: end of file not reached\n";
    
    for (auto & row : contents)
    {
        edge_list.emplace_back(0,row[0], row[1]);
        // edges.emplace_back(row);
        //std::copy(row.begin(), row.end(), std::ostream_iterator<double>(std::cout," "));
        //std::cout << "\n";
    }
    // for (auto& row : random)
        // edges.emplace_back(row);
	
	// std::cout << "Input complex is of dimension " << edge_list.dimension() << " with " << edge_list.num_edges() << " maximal edges, " << edge_list.filtration_simplex_range().size() << " edges and " << edge_list.num_vertices() << " vertices." << std::endl;

	clock_t edge_list_formed = clock();
	std::cout << "Now going for matrix formation" << std::endl;

	FlagComplexSpMatrix mat(5,edge_list);
    // mat.contraction(2,6);
    // mat.contraction(6,4);
    // mat.contraction(5,8);
    
    // SparseMsMatrix mat2(10,20);
	clock_t matrix_formed = clock();
	std::cout << "Matrix formed ... Now going for collapse" << std::endl;
 
    // for (auto& row : edges)
    // {
        
    //     std::cout << " checking membership for the simplex " << ": " ;
    //     for(auto & v: row)
    //             std::cout << v << ", " ;
    //     std::cout << std::endl;   
        
    //     if(mat.membership(row))
  
    //         std::cout << " It exists..." <<std::endl;
    //     else
    //         std::cout << " It doesn't  exist..." <<std::endl;
    //     //std::copy(row.begin(), row.end(), std::ostream_iterator<double>(std::cout," "));
    //     //std::cout << "\n";
    // }  
    // std::cout << "****************************************************" << std::endl; 
    // auto v = mat.active_strong_expansion(2,5,.2);
    // std::cout << "The contracted vertex is :" << v <<"; " << std::endl;
    // auto w = mat.active_strong_expansion(2,6,.2);
    // std::cout << "The contracted vertex is :" << w <<"; " << std::endl;
    // mat.contraction(4,3);
    // std::cout << "Manually contracted the edge [4,3]" << std::endl;
    mat.strong_collapse();
    // Fake_simplex_tree coll_tree = mat.collapsed_tree();
    clock_t collapse_done = clock();
    // for(auto & v: mat.vertex_set())
    // {   
    //     std::cout << "Active relative neighbors of the vertex 2 and : " << v << " are, " ;
    //     for(auto & w: mat.active_strong_expansion(2,5))
    //         std::cout << w << " ";
    //     std::cout << std::endl;
    // }
     for(auto & v: mat.all_edges())
    {   
        std::cout << "The current edge in the complex after the collapse are: " ;
        // for(auto & w: v)
            std::cout << std::get<0>(v) << " " << std::get<1>(v);
        std::cout << std::endl;
    } 
  

    // for (auto& row : edges)
    // {
        
    //     std::cout << " checking membership for the simplex " << ": " ;
    //     for(auto & v: row)
    //             std::cout << v << ", " ;
    //     std::cout << std::endl;   
        
    //     if(mat.membership(row))
  
    //         std::cout << " It exists..." <<std::endl;
    //     else
    //         std::cout << " It doesn't  exist..." <<std::endl;
    //     //std::copy(row.begin(), row.end(), std::ostream_iterator<double>(std::cout," "));
    //     //std::cout << "\n";
    // } 
    
	std::cout << "Collapse done !" << std::endl;

	std::cout << "Time for formation of Matrix : " << (matrix_formed - edge_list_formed)/CLOCKS_PER_SEC << " seconds" << std::endl;
	std::cout << "Time for Collapse : " << (collapse_done - matrix_formed)/CLOCKS_PER_SEC  << " seconds" << std::endl;

	// std::cout << "Collapsed Rips complex is of dimension " << coll_tree.dimension() << " with " << coll_tree.num_edges() << " maximal edges " << coll_tree.filtration_simplex_range().size() << " edges and "  << coll_tree.num_vertices() << " vertices." << std::endl;

    // mat.contraction(8,19);
	return 0;
}


