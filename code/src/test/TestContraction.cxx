/*
 * TestContraction.cxx
 *
 *  Created on: Feb 11, 2014
 *      Author: dsalinas
 */

#include "Skeleton_blocker/Simple_skeleton_blockers_traits.h"
#include "Skeleton_blocker/Simplex.h"
#include "contraction/Skeleton_blocker_contractor.h"
#include "Utils.h"
#include "iofile.h"
#include "Test.h"
#include "Skeleton_blocker_geometric_complex.h"
//#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>


#include "contraction/Edge_profile.h"

#include "contraction/policies/Cost_policy.h"
#include "contraction/policies/Edge_length_cost.h"
#include "contraction/policies/Placement_policy.h"
#include "contraction/policies/Middle_placement.h"

#include "contraction/policies/Valid_contraction_policy.h"
#include "contraction/policies/Dummy_valid_contraction.h"
#include "contraction/policies/Link_condition_valid_contraction.h"


using namespace std;

struct Geometry_trait{
	typedef std::vector<double> Point;
};

typedef Geometry_trait::Point Point;

typedef Simple_complex_geometry_traits<Geometry_trait> Complex_geometric_traits;

typedef Skeleton_blocker_geometric_complex< Complex_geometric_traits > Complex;

typedef typename Complex::Vertex_handle Vertex_handle;
typedef typename Complex::Root_vertex_handle Root_vertex_handle;

using namespace contraction;

typedef Skeleton_blocker_contractor<Complex> Complex_contractor;

typedef Edge_profile<Complex> Profile;




int main (int argc, char *argv[])
{
	if (argc!=2){
		std::cerr << "Usage "<<argv[0]<<" ../src/data/min_cube3.off if the file src/data/min_cube3.off exists.\n";
		return -1;
	}
	Complex complex;
	std::string file_name(argv[1]);
	system("pwd");
	bool loaded = read_off_file<Complex>(file_name,complex);
	PRINT(complex.to_string());
	Complex_contractor contractor(complex);

	contractor.contract_edges(complex.get_num_vertices()-1);

	PRINT(complex.to_string());

}


