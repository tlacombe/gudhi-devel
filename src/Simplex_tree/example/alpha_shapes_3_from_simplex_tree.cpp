#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Alpha_shape_3.h>
#include <CGAL/iterator.h>

#include <fstream>
#include <list>
#include <cassert>

#include <vector>
#include <boost/variant.hpp>
#include <boost/optional.hpp>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;

typedef CGAL::Alpha_shape_vertex_base_3<Kernel>          Vb;
typedef CGAL::Alpha_shape_cell_base_3<Kernel>            Fb;
typedef CGAL::Triangulation_data_structure_3<Vb,Fb>  Tds;
typedef CGAL::Delaunay_triangulation_3<Kernel,Tds>       Triangulation_3;
typedef CGAL::Alpha_shape_3<Triangulation_3>         Alpha_shape_3;

typedef Kernel::Point_3                                  Point;
typedef Alpha_shape_3::Alpha_iterator               Alpha_iterator;


int main (int argc, char * const argv[])
{
  if (argc != 2) {
    std::cerr << "Usage: " << argv[0]
      << " path_to_file_graph \n";
    return 0;
  }
  std::string filegraph   = argv[1];

  std::list<Point> lp;

  //read input
  std::ifstream is(filegraph.c_str());
  int n;
  is >> n;
  std::cout << "Reading " << n << " points " << std::endl;
  Point p;
  for( ; n>0 ; n--)    {
    is >> p;
    lp.push_back(p);
  }

  // compute alpha shape
  Alpha_shape_3 as(lp.begin(),lp.end());
  std::cout << "Alpha shape computed in REGULARIZED mode by default" << std::endl;

  std::vector<CGAL::Object> the_objects;
  std::vector<Alpha_shape_3::FT> the_ft;

  typedef CGAL::Dispatch_output_iterator<
    CGAL::cpp11::tuple<CGAL::Object, Alpha_shape_3::FT>,
    CGAL::cpp11::tuple<std::back_insert_iterator< std::vector<CGAL::Object> >,
                       std::back_insert_iterator< std::vector<Alpha_shape_3::FT> >
                       > > Dispatch;

  Dispatch disp = CGAL::dispatch_output<CGAL::Object, Alpha_shape_3::FT>(
    std::back_inserter(the_objects),
    std::back_inserter(the_ft));

  as.filtration_with_alpha_values(disp);
  std::cout << "filtration_with_alpha_values returns : " << the_objects.size() << " objects" << std::endl;
  std::cout << "filtration_with_alpha_values returns : " << the_ft.size() << " FT" << std::endl;

  typename Alpha_shape_3::Cell_handle cell;
  typename Alpha_shape_3::Facet      facet;
  typename Alpha_shape_3::Edge       edge;
  typename Alpha_shape_3::Vertex_handle  vertex;

  typename Alpha_shape_3::size_type count_vertices = 0;
  typename Alpha_shape_3::size_type count_edges = 0;
  typename Alpha_shape_3::size_type count_facets = 0;
  typename Alpha_shape_3::size_type count_cells = 0;

  for(auto object_iterator: the_objects)
  {
	  if(CGAL::assign(cell,object_iterator))
		  count_cells++;
		  //std::cout << "intersection object is a cell=" << std::endl;
	  if(CGAL::assign(facet,object_iterator))
		  count_facets++;
		  //std::cout << "intersection object is a facet=" << std::endl;
	  if(CGAL::assign(edge,object_iterator))
		  count_edges++;
		  //std::cout << "intersection object is a edge=" << std::endl;
	  if(CGAL::assign(vertex,object_iterator))
		  count_vertices++;
		  //std::cout << "intersection object is a edge=" << std::endl;

  }
  std::cout << "vertices \t\t" << count_vertices << std::endl;
  std::cout << "edges \t\t"    << count_edges << std::endl;
  std::cout << "facets \t\t"   << count_facets << std::endl;
  std::cout << "cells \t\t"    << count_cells << std::endl;

  /*for(auto ft_iterator: the_ft)
  {
    if (ft_iterator < *opt)
	  std::cout << "Approx=" << ft_iterator << std::endl;
  }*/


  return 0;
}
