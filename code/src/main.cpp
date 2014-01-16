#include <iostream>
#include <ctime>
#include "iofile.h"
#include "Flag_simplex_tree.h"

#include "Euclidean_rips_naive_geometry_traits.h"

using namespace std;


int main (int argc, char * const argv[]) 
{
	// Extract data points from file file_name.
	// Turn them into a Point_range object:	
	
	Euclidean_rips_naive_geometry_traits::Point_range points;
	string file_name = "/Users/cmaria/Desktop/Bro_fp.txt";  
 	read_points(file_name,points);                         //<- make different versions of read_points

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
	int			dim_max			= 10;
	double	rho_max			= 0.77;

	// Creates an empty simplex tree
	Flag_simplex_tree< Euclidean_rips_naive_geometry_traits > st;
	
	// Initializes the geometry traits with the Point_range
	st.gt()->init(points);
	
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
	
	st.init(//points,			//Point_range, subset of the points in st.gt()
					dim_max,			// dimension maximal of the expansion
					rho_max);			// threshold for the Rips
	
	// Prints the simplex tree in std::cout
	cout << st.size_complex() << endl;
		
	
	
	
	
	/* TO DO
	 PcoH<Flag_simplex_tree> co;
	 co.init(rho_max);	// CAM->does nothing
											// PH initialize the matrix
	 intervals_begin();	// <- does the work
	//end
	*/
	return EXIT_SUCCESS;
}
