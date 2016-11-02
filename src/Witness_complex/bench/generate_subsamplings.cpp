#include <iostream>
#include <string>
#include <cstdlib>
//#include <ctime>
#include <stdio.h>
#include "output.h"

#include <gudhi/choose_n_farthest_points.h>
#include <gudhi/Points_off_io.h>

int main(int argc, char* const argv[])
{
  if (argc != 3) {
    std::cout << "No filename provided.\n";
    return 1;
  }
  char* file_name = argv[1];
  unsigned nb_files = atoi(argv[2]);

  // Read the point file
  Point_Vector point_vector;
  Gudhi::Points_off_reader<Point_d> off_reader(file_name);
  if (!off_reader.is_valid()) {
      std::cerr << "Witness complex - Unable to read file " << file_name << "\n";
      exit(-1);  // ----- >>
    }
  point_vector = Point_Vector(off_reader.get_point_cloud());

  std::size_t nbW = point_vector.size();
  std::cout << "|W|=" << point_vector.size() << "\n";
  std::cout << "d=" << point_vector[0].dimension() << "\n";
  std::cout << "file name=" << file_name << "\n";
  
  for (unsigned i = 0; i <= nb_files; ++i) {
    // Choose landmarks
    Point_Vector landmarks;
    std::size_t nbL = (i*nbW)/nb_files;
    Gudhi::subsampling::choose_n_farthest_points(K(), point_vector, nbL, std::back_inserter(landmarks));
    char* new_file_name = (char*)std::malloc(255*sizeof(char));   // !!!! Potential error
    sprintf(new_file_name, "%s_%lu.off", file_name, landmarks.size());
    write_points_to_off(new_file_name, point_vector);
    std::free(new_file_name);
  }
}
