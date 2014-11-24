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

typedef CGAL::Exact_predicates_inexact_constructions_kernel Gt;

typedef CGAL::Alpha_shape_vertex_base_3<Gt>          Vb;
typedef CGAL::Alpha_shape_cell_base_3<Gt>            Fb;
typedef CGAL::Triangulation_data_structure_3<Vb,Fb>  Tds;
typedef CGAL::Delaunay_triangulation_3<Gt,Tds>       Triangulation_3;
typedef CGAL::Alpha_shape_3<Triangulation_3>         Alpha_shape_3;

typedef Gt::Point_3                                  Point;
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
  std::cout << "Alpha shape computed in REGULARIZED mode by default"
	    << std::endl;

  // find optimal alpha value
  Alpha_iterator opt = as.find_optimal_alpha(1);
  std::cout << "Optimal alpha value to get one connected component is "
	    <<  *opt    << std::endl;
  as.set_alpha(*opt);

  assert(as.number_of_solid_components() == 1);

  std::cout << "nb of alphas:" << as.number_of_alphas() << std::endl;

  std::list<Alpha_shape_3::Cell_handle> cells;
  std::list<Alpha_shape_3::Facet>       facets;
  std::list<Alpha_shape_3::Edge>        edges;
  std::list<Alpha_shape_3::Vertex_handle>     vertices;
  as.get_alpha_shape_cells(std::back_inserter(cells),Alpha_shape_3::INTERIOR);

  as.get_alpha_shape_facets(std::back_inserter(facets),Alpha_shape_3::REGULAR);

  as.get_alpha_shape_facets(std::back_inserter(facets),Alpha_shape_3::SINGULAR);

  as.get_alpha_shape_edges(std::back_inserter(edges),Alpha_shape_3::SINGULAR);
  as.get_alpha_shape_vertices(std::back_inserter(vertices),Alpha_shape_3::EXTERIOR);

  std::cout << " The 0-shape has : " << std::endl;
  std::cout << cells.size() << " interior tetrahedra" << std::endl;
  std::cout << facets.size() << " boundary facets" << std::endl;
  std::cout << edges.size()  << " singular edges" << std::endl;
  std::cout << vertices.size() << " vertices" << std::endl;

  //i'm not sure if this is the correct way to access the filtration
  std::list<CGAL::Object>       objects;
  as.filtration(std::back_inserter(objects));
  std::cout << "filtration returns : " << objects.size() << " objects" << std::endl;
  std::list<CGAL::Object>::iterator pos;
  pos = objects.begin();
  while (pos != objects.end()){
	  pos++;
  }

  as.get_alpha_shape_cells(std::back_inserter(cells),Alpha_shape_3::INTERIOR);

  as.get_alpha_shape_facets(std::back_inserter(facets),Alpha_shape_3::REGULAR);

  as.get_alpha_shape_facets(std::back_inserter(facets),Alpha_shape_3::SINGULAR);

  as.get_alpha_shape_edges(std::back_inserter(edges),Alpha_shape_3::SINGULAR);
  as.get_alpha_shape_vertices(std::back_inserter(vertices),Alpha_shape_3::EXTERIOR);

  std::cout << " The 0-shape has : " << std::endl;
  std::cout << cells.size() << " interior tetrahedra" << std::endl;
  std::cout << facets.size() << " boundary facets" << std::endl;
  std::cout << edges.size()  << " singular edges" << std::endl;
  std::cout << vertices.size() << " vertices" << std::endl;


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

  as.get_alpha_shape_cells(std::back_inserter(cells),Alpha_shape_3::INTERIOR);

  as.get_alpha_shape_facets(std::back_inserter(facets),Alpha_shape_3::REGULAR);

  as.get_alpha_shape_facets(std::back_inserter(facets),Alpha_shape_3::SINGULAR);

  as.get_alpha_shape_edges(std::back_inserter(edges),Alpha_shape_3::SINGULAR);
  as.get_alpha_shape_vertices(std::back_inserter(vertices),Alpha_shape_3::EXTERIOR);

  std::cout << " The 0-shape has : " << std::endl;
  std::cout << cells.size() << " interior tetrahedra" << std::endl;
  std::cout << facets.size() << " boundary facets" << std::endl;
  std::cout << edges.size()  << " singular edges" << std::endl;
  std::cout << vertices.size() << " vertices" << std::endl;


  return 0;
}
