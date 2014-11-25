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

  /*for(auto object_iterator: the_objects)
  {*/
  std::vector<CGAL::Object>::iterator object_iterator = the_objects.begin();
  std::cout << "object_iterator type=" << (*object_iterator).type().name() << std::endl;
  Vb type1;
  if(CGAL::assign(type1,*object_iterator))
	  std::cout << "intersection object is a type1=" << std::endl;
  Fb type2;
  if(CGAL::assign(type2,*object_iterator))
	  std::cout << "intersection object is a type2=" << std::endl;
  Tds type3;
  if(CGAL::assign(type3,*object_iterator))
	  std::cout << "intersection object is a type3=" << std::endl;
   /* }
   */
  /*for(auto ft_iterator: the_ft)
  {
    if (ft_iterator < *opt)
	  std::cout << "Approx=" << ft_iterator << std::endl;
  }*/


  return 0;
}
