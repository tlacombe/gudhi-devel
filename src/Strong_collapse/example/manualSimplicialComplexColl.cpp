#include <gudhi/SparseMsMatrix.h>
#include <gudhi/VertexMaximalMaps.h>
#include <gudhi/Fake_simplex_tree.h>

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
using Filtration_value = Fake_simplex_tree::Filtration_value;
using Rips_complex = Gudhi::rips_complex::Rips_complex<Filtration_value>;

using Fake_simplex_tree  = Gudhi::Fake_simplex_tree ;
using Vertex             = Fake_simplex_tree::Vertex;

using Vector_of_points = std::vector<Point>;

int main()
{
	
	Vector_of_points points;
	//Vector_of_points noisePoints;
	Fake_simplex_tree stree;
	std::vector<std::vector<std::size_t>> contents;
    std::vector<std::vector<std::size_t>> simplices;
    std::cout << "Please enter the #vertices and then the vertices of the maximal simplices in the following format" <<std::endl;
    std::cout << "3" << std::endl;
    std::cout << "1 2 3" << std::endl;
    std::cout << "Enter s to stop" << std::endl;
    int number;
    std::vector<std::vector<std::size_t>> random = { {1,4},{3,4},{4,8}, {1,8},{3,8},{89}, {6}, {2},{4},{5} }; //{{89,77,29},{2,4,6}, {2,4,5},{3,4},{2,3},{6},{2,7},{7,5}, {2,3,4,5}};
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
        stree.insert_simplex_and_subfaces(row);
        simplices.emplace_back(row);
        //std::copy(row.begin(), row.end(), std::ostream_iterator<double>(std::cout," "));
        //std::cout << "\n";
    }
    for (auto& row : random)
        simplices.emplace_back(row);
	
	std::cout << "Input complex is of dimension " << stree.dimension() << " with " << stree.num_simplices() << " maximal simplices, " << stree.filtration_simplex_range().size() << " simplices and " << stree.num_vertices() << " vertices." << std::endl;

	clock_t stree_formed = clock();
	std::cout << "Now going for matrix formation" << std::endl;

	SparseMsMatrix mat(stree);
    VertexMaximalMaps vertMxMaps(stree);
    mat.contraction(2,6);
    mat.contraction(6,4);
    mat.contraction(5,8);
    
    // SparseMsMatrix mat2(10,20);
	clock_t matrix_formed = clock();
	std::cout << "Matrix formed ... Now going for collapse" << std::endl;
 
    // for (auto& row : simplices)
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
	
    mat.strong_collapse();
    // Fake_simplex_tree coll_tree = mat.collapsed_tree();
    clock_t collapse_done = clock();
  

    for (auto& row : simplices)
    {
        
        std::cout << " checking membership for the simplex " << ": " ;
        for(auto & v: row)
                std::cout << v << ", " ;
        std::cout << std::endl;   
        
        if(mat.membership(row))
  
            std::cout << " It exists..." <<std::endl;
        else
            std::cout << " It doesn't  exist..." <<std::endl;
        //std::copy(row.begin(), row.end(), std::ostream_iterator<double>(std::cout," "));
        //std::cout << "\n";
    } 
    
	std::cout << "Collapse done !" << std::endl;

	std::cout << "Time for formation of Matrix : " << (matrix_formed - stree_formed)/CLOCKS_PER_SEC << " seconds" << std::endl;
	std::cout << "Time for Collapse : " << (collapse_done - matrix_formed)/CLOCKS_PER_SEC  << " seconds" << std::endl;

	// std::cout << "Collapsed Rips complex is of dimension " << coll_tree.dimension() << " with " << coll_tree.num_simplices() << " maximal simplices " << coll_tree.filtration_simplex_range().size() << " simplices and "  << coll_tree.num_vertices() << " vertices." << std::endl;

    // mat.contraction(8,19);
	return 0;
}


