#include <sys/types.h>
#include <sys/stat.h>

#include <gudhi/Simplex_tree.h>
#include <gudhi/Witness_complex.h>
#include <gudhi/pick_n_random_points.h>
#include <gudhi/Points_off_io.h>

#include <gudhi/Persistent_cohomology.h>
#include <gudhi/Graph_matching.h>

#include <CGAL/Epick_d.h>

#include <iostream>
#include <fstream>
#include <ctime>
#include <string>
#include <vector>
#include <fstream>
#include <stdio.h>
#include <time.h>

typedef CGAL::Epick_d<CGAL::Dynamic_dimension_tag> K;
typedef typename K::Point_d Point_d;
typedef typename Gudhi::witness_complex::Witness_complex<K> Witness_complex;
typedef std::vector< Vertex_handle > typeVectorVertex;
typedef std::vector< Point_d > Point_vector;

typedef Gudhi::Simplex_tree<> SimplexTree;
typedef typename Gudhi::persistent_cohomology::Persistent_cohomology<SimplexTree, Gudhi::persistent_cohomology::Field_Zp> PersistentCohomology;
typedef std::vector<std::pair<double, double> > GraphTable;

GraphTable read_diagram_from_file( const char* filename )
{
    std::ifstream in;
    in.open( filename );
    std::vector< std::pair<double, double> > result;
    if ( !in.is_open() )
    {
        std::cerr << "File : " << filename << " do not exist. The program will now terminate \n";
        throw "File do not exist \n";
    }

    std::string line;
    while (!in.eof())
    {
        getline(in,line);
        if ( line.length() != 0 )
        {
            std::stringstream lineSS;
            lineSS << line;
            double beginn, endd;
            lineSS >> beginn;
            lineSS >> endd;
            result.push_back( std::make_pair( beginn , endd ) );
        }
    }
    in.close();
    return result;
}  //read_diagram_from_file


void write_graph(SimplexTree& simplextree, PersistentCohomology& pers, GraphTable& graph)
{
  std::time_t raw_time;
  raw_time = time(NULL);
  char* filename;
  sprintf(filename, "temp_%s", asctime(localtime(&raw_time)));
  //std::ofstream f_out(filename);
  //graph = read_diagram_from_file(filename);
  for (auto pair: pers.get_persistent_pairs())
    if (simplextree.dimension(std::get<0>(pair)) == 1) {
      double alpha1 = simplextree.filtration(std::get<0>(pair));
      double alpha2 = simplextree.filtration(std::get<1>(pair));
      graph.push_back(std::make_pair(alpha1,alpha2));
    }
}

int main(int argc, char * const argv[]) {
  if (argc != 5) {
    std::cerr << "Usage: " << argv[0]
        << " path_to_point_file number_of_landmarks max_squared_alpha limit_dimension\n";
    return 0;
  }

  std::string file_name = argv[1];
  int nbL = atoi(argv[2]), lim_dim = atoi(argv[4]);
  double alpha2 = atof(argv[3]);
  clock_t start, end;
  SimplexTree simplex_tree1, simplex_tree2;

  // Read the point file
  Point_vector point_vector, landmarks1, landmarks2;
  Gudhi::Points_off_reader<Point_d> off_reader(file_name);
  if (!off_reader.is_valid()) {
      std::cerr << "Witness complex - Unable to read file " << file_name << "\n";
      exit(-1);  // ----- >>
    }
  point_vector = Point_vector(off_reader.get_point_cloud());
  
  std::cout << "Successfully read " << point_vector.size() << " points.\n";
  std::cout << "Ambient dimension is " << point_vector[0].dimension() << ".\n";

  // Choose landmarks
  Gudhi::subsampling::pick_n_random_points(point_vector, nbL, std::back_inserter(landmarks1));
  Gudhi::subsampling::pick_n_random_points(point_vector, nbL, std::back_inserter(landmarks2));
  
  // Compute witness complex
  start = clock();
  Witness_complex witness_complex1(landmarks1,
                                   point_vector);

  witness_complex1.create_complex(simplex_tree1, alpha2, lim_dim);
  end = clock();
  std::cout << "Witness complex 1 took "
      << static_cast<double>(end - start) / CLOCKS_PER_SEC << " s. \n";
  std::cout << "Number of simplices is: " << simplex_tree1.num_simplices() << "\n";

  start = clock();
  Witness_complex witness_complex2(landmarks2,
                                   point_vector);

  witness_complex2.create_complex(simplex_tree2, alpha2, lim_dim);
  end = clock();
  std::cout << "Witness complex 2 took "
      << static_cast<double>(end - start) / CLOCKS_PER_SEC << " s. \n";
  std::cout << "Number of simplices is: " << simplex_tree2.num_simplices() << "\n";

  PersistentCohomology pers1(simplex_tree1), pers2(simplex_tree2);
  GraphTable graph1, graph2;
  write_graph(simplex_tree1, pers1, graph1);
  write_graph(simplex_tree2, pers2, graph2);

  double bd = Gudhi::bottleneck::bottleneck_distance(graph1, graph2);
  std::cout << "The bottleneck distance is: " << bd << "\n";
  
}
