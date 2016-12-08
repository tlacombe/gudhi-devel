/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Pawel Dlotko
 *
 *    Copyright (C) 2015  INRIA Saclay (France)
 *
 *    This program is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *
 *    This program is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */


#include <gudhi/reader_utils.h>
#include <gudhi/Bitmap_cubical_complex_base.h>
#include <gudhi/Bitmap_cubical_complex_periodic_boundary_conditions_base.h>
#include <gudhi/Bitmap_cubical_complex.h>
#include <gudhi/Compute_persistence_with_phat.h>

using namespace Gudhi;
using namespace Gudhi::cubical_complex;
using namespace Gudhi::phat_interface;

//standard stuff
#include <iostream>
#include <sstream>
#include <vector>

using namespace std;

int main(int argc, char** argv) 
{
  std::cout << "The program compute persistence of a bitmap cubical complex taking into account periodic boundary conditions. See documentation of a bitmap cubical complex for furhter details.\n";
  if (argc != 2) {
    cout << "Wrong number of parameters. Please provide the name of a file with a Perseus style bitmap at the "
        << "input. The program will now terminate.\n";
    return 1;
  }


  typedef Gudhi::Cubical_complex::Bitmap_cubical_complex_periodic_boundary_conditions_base<double> Bitmap_base;
  typedef Gudhi::Cubical_complex::Bitmap_cubical_complex< Bitmap_base > Bitmap_cubical_complex;
  Bitmap_cubical_complex b(argv[1]);  
  Compute_persistence_with_phat< Bitmap_cubical_complex > phat(&b); 
  
  phat.compute_persistence_pairs();
 
  std::pair< std::vector< std::vector<double> >,
  std::vector< std::vector< std::pair<double, double> > > > persistence = phat.get_the_intervals();

  //this is for Gudhi format of outputting intervals:
  std::stringstream ss;
  ss << argv[1] << "_persistence";

  std::cerr << "Filename : " << ss.str().c_str() << std::endl;

  write_intervals_to_file_Gudhi_format<double>(persistence, ss.str().c_str(), -1);

  return 0;
}

