#include <iostream>
#include <ctime>

#include "iofile.h"

//#include "Simplex_tree.h"

#include <boost/pending/disjoint_sets.hpp>
#include <boost/iterator/counting_iterator.hpp>
#include "boost/range/counting_range.hpp"

#include "boost/container/flat_map.hpp"

#include "Euclidean_geometry.h"
#include "Rips_graph_naive.h"
#include "Simplex_tree.h"


using namespace std;



int main (int argc, char * const argv[]) 
{

	// Extract data points from file file_name.
	// Turn them into a Point_range object:	   <--- has to be templated; best, use a istream_iterator
	typedef std::vector<double> Point;           // in order not to copy many times the Points.
  typedef std::vector<Point>  Point_range;

  Point_range points;
	string file_name = "/Users/cmaria/Desktop/Points.txt";
	read_points(file_name,points);

  // Create a metric space from the points, with euclidean metric:
	Euclidean_geometry< Point > ms;
  ms.init(points);                 

  // Create a NeighborGraph with the space:
  double threshold = 100;
  Rips_graph_naive< Euclidean_geometry< Point > > ng(ms,threshold);  // <--- merge Euclidean_geom and Rips_graph ?

  // Create a simplex_tree
  Simplex_tree< Euclidean_geometry< Point > > st(ms); //constructor

  st.insert_graph(ng); //insert the graph

  int max_dim = 10;
  std::cout << "Expand the flag complex \n";
  st.expansion(max_dim);

  std::cout << st << std::endl;



  /****************************************/
/*	cout << "Points from file: ---------- \n";
	for(vector< vector< double > >::iterator it = points.begin();
			it!=points.end(); it++)
	{
		for(vector<double>::iterator p_it = it->begin();
				p_it!=it->end(); p_it++)
		{
			cout << *p_it << " ";	
		}
		cout << endl;
	}
	cout << "---------------------------- \n";
*/	/****************************************/
	
	/****** Construct a Rips complex: ******/
/*	int	dim_max			= 100;
	double	rho_max			= 100;

	// Creates an empty simplex tree
	Flag_simplex_tree< Euclidean_rips_naive_geometry_traits > st;
	
	// Initializes the geometry traits with the Point_range
	st.gt()->init(points);
*/	
	/****************************************/
/*	cout << "Points from metric space: -- \n";
	for(vector< vector< double > >::iterator it = st.gt()->point_range().begin();
			it != st.gt()->point_range().end(); ++it)
	{
		for(vector<double>::iterator p_it = it->begin();
				p_it!=it->end(); p_it++)
		{
			cout << *p_it << " ";	
		}
		cout << endl;
	}
	cout << "Nb elements MP: " << st.gt()->nb_elements() << endl;
	cout << "---------------------------- \n";
*/	/****************************************/
	
	// Constructs the Rips complex
	
/*	st.init(//points,			//Point_range, subset of the points in st.gt()
					dim_max,			// dimension maximal of the expansion
					rho_max);			// threshold for the Rips
	
	// Prints the simplex tree in std::cout
	cout << st.size_complex() << endl;
		
	
//	return 0;
	
//	std::cout << st;
	
	vector< int > v;
	for(int i = 1; i < 4 ; ++i) v.push_back(i);
	Flag_simplex_tree< Euclidean_rips_naive_geometry_traits >::Simplex_handle sh; 
	sh = st.find(v);
	st.print(sh);
	std::cout << "\n \n";


Flag_simplex_tree< Euclidean_rips_naive_geometry_traits >::
Boundary_simplex_range rg = st.boundary_simplex_range(sh);
	
Flag_simplex_tree< Euclidean_rips_naive_geometry_traits >::
		Boundary_simplex_iterator it = rg.begin();



	for(; it != rg.end(); ++it)
	{
		st.print(*it);
	}

*/	
	
	/* TO DO
	 PcoH<Flag_simplex_tree> co;
	 co.init(rho_max);	// CAM->does nothing
											// PH initialize the matrix
	 intervals_begin();	// <- does the work
	//end
	*/
	return EXIT_SUCCESS;
}
