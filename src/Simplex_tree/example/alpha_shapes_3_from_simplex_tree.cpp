#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Alpha_shape_3.h>
#include <CGAL/iterator.h>

#include <fstream>


// Alpha_shape_3 templates type definitions
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Alpha_shape_vertex_base_3<Kernel>             Vb;
typedef CGAL::Alpha_shape_cell_base_3<Kernel>               Fb;
typedef CGAL::Triangulation_data_structure_3<Vb,Fb>         Tds;
typedef CGAL::Delaunay_triangulation_3<Kernel,Tds>          Triangulation_3;
typedef CGAL::Alpha_shape_3<Triangulation_3>                Alpha_shape_3;

// From file type definition
typedef Kernel::Point_3                                     Point;

// filtration with alpha values needed type definition
typedef Alpha_shape_3::FT Alpha_value_type;
typedef CGAL::Object      Object;
typedef CGAL::Dispatch_output_iterator<
  CGAL::cpp11::tuple<Object, Alpha_value_type>,
  CGAL::cpp11::tuple<std::back_insert_iterator< std::vector<Object> >, std::back_insert_iterator< std::vector<Alpha_value_type> >
                     > > Dispatch;
typedef Alpha_shape_3::Cell_handle   Cell_handle;
typedef Alpha_shape_3::Facet         Facet;
typedef Alpha_shape_3::Edge          Edge;
typedef Alpha_shape_3::Vertex_handle Vertex_handle;

int main (int argc, char * const argv[])
{
	// program args management
  if (argc != 2) {
    std::cerr << "Usage: " << argv[0]
      << " path_to_file_graph \n";
    return 0; // ----- >>
  }

  // Read points from file
  std::string filegraph   = argv[1];
  std::list<Point> lp;
  std::ifstream is(filegraph.c_str());
  int n;
  is >> n;
  std::cout << "Reading " << n << " points " << std::endl;
  Point p;
  for( ; n>0 ; n--)    {
    is >> p;
    lp.push_back(p);
  }

  // alpha shape construction from points
  Alpha_shape_3 as(lp.begin(),lp.end());
  std::cout << "Alpha shape computed in REGULARIZED mode by default" << std::endl;

  // filtration with alpha values from alpha shape
  std::vector<Object> the_objects;
  std::vector<Alpha_value_type> the_alpha_values;

  Dispatch disp = CGAL::dispatch_output<Object, Alpha_value_type>( std::back_inserter(the_objects), std::back_inserter(the_alpha_values));

  as.filtration_with_alpha_values(disp);
  std::cout << "filtration_with_alpha_values returns : " << the_objects.size() << " objects" << std::endl;

  Alpha_shape_3::size_type count_vertices = 0;
  Alpha_shape_3::size_type count_edges    = 0;
  Alpha_shape_3::size_type count_facets   = 0;
  Alpha_shape_3::size_type count_cells    = 0;

  // Loop on objects vector
  for(auto object_iterator: the_objects)
  {
	  if (const Cell_handle* cell = CGAL::object_cast<Cell_handle>(&object_iterator))
	  {
		  std::cout << "Cell[" << (*cell)->vertex(0)->point() << ";" << (*cell)->vertex(1)->point() << ";" << (*cell)->vertex(2)->point() << ";" << (*cell)->vertex(3)->point() << "]" << std::endl;
		  count_cells++;
	  }
	  if (const Facet* facet = CGAL::object_cast<Facet>(&object_iterator))
	  {
		  count_facets++;
	  }
	  if (const Edge* edge = CGAL::object_cast<Edge>(&object_iterator))
	  {
		  count_edges++;
	  }
	  if (const Vertex_handle* vertex = CGAL::object_cast<Vertex_handle>(&object_iterator))
	  {
		  count_vertices++;
		  std::cout << "intersection object is a vertex=" << (*vertex)->point() << std::endl;
	  }

  }
  std::cout << "vertices \t\t" << count_vertices << std::endl;
  std::cout << "edges \t\t"    << count_edges << std::endl;
  std::cout << "facets \t\t"   << count_facets << std::endl;
  std::cout << "cells \t\t"    << count_cells << std::endl;

  /*for(auto alpha_value_iterator: the_alpha_values)
  {
    if (alpha_value_iterator < *opt)
	  std::cout << "Approx=" << alpha_value_iterator << std::endl;
  }*/


  return 0;
}
