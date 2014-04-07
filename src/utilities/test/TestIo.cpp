
#include "io.h"

bool test_points_extraction(std::string filename)
{
  std::vector< std::vector< double > > points; //read the points from the file
  read_points( filename, points );             //and turn them into a Point_range
  return true;
}

int main()
{
  return 0;
}