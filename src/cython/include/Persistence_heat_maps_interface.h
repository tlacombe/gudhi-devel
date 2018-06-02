/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Pawel Dlotko
 *
 *    Copyright (C) 2018 Swansea University, UK
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

#ifndef PERSISTENCE_HEAT_MAPS_INTERFACE_H_
#define PERSISTENCE_HEAT_MAPS_INTERFACE_H_

#include <gudhi/Persistence_heat_maps.h>

namespace Gudhi {
namespace Persistence_representations {

class Persistence_heat_maps_interface : public Persistence_heat_maps<constant_scaling_function> {
 public:

  Persistence_heat_maps_interface():Persistence_heat_maps(){}

 
  Persistence_heat_maps_interface(const std::vector<std::pair<double, double> >& interval,
						size_t how_many_pixels_raidus_of_Gausian_kernel,
                        size_t number_of_pixels,
                        double min_,
                        double max_
                        ):
   Persistence_heat_maps(interval,create_Gaussian_filter(how_many_pixels_raidus_of_Gausian_kernel, 1),
						this->erase_below_diagonal,number_of_pixels,((min_ == max_) ?  -1 : min_) , ((min_ == max_) ? -1 : max_)){}   			
									                
  
  Persistence_heat_maps_interface(const char* filename, 
						size_t how_many_pixels_raidus_of_Gausian_kernel,
                        size_t number_of_pixels,
                        double min_,
                        double max_,
                        int dimension
                        ):
  Persistence_heat_maps(filename,create_Gaussian_filter(how_many_pixels_raidus_of_Gausian_kernel, 1),
						this->erase_below_diagonal,number_of_pixels, ((min_ == max_) ?  -1 : min_) , ((min_ == max_) ? -1 : max_),dimension){}
						
  //****************
  static Persistence_heat_maps_interface* construct_from_file(  const char* filename, size_t how_many_pixels_raidus_of_Gausian_kernel,
																size_t number_of_pixels, double min_ = -1,
																double max_ = -1, int dimensions = -1 )
  {	  	  
      Persistence_heat_maps_interface* result = new Persistence_heat_maps_interface(filename,how_many_pixels_raidus_of_Gausian_kernel,number_of_pixels,min_,max_,dimensions);
	  return result;  
  }
  static Persistence_heat_maps_interface* construct_from_vector_of_pairs( const std::vector<std::pair<double, double> >& interval, size_t how_many_pixels_raidus_of_Gausian_kernel, size_t number_of_pixels, double min_ = 0, double max_ = 0 )
  {
      Persistence_heat_maps_interface* result = new Persistence_heat_maps_interface(interval, how_many_pixels_raidus_of_Gausian_kernel, number_of_pixels, min_, max_);
	  return result;  
  }
  
  void compute_mean_interface(const std::vector<Persistence_heat_maps_interface*>& maps_)
  {
	  std::vector<Persistence_heat_maps*> maps;
	  maps.reserve( maps_.size() );
	  for ( size_t i = 0 ; i != maps_.size() ; ++i )
	  {
		  maps.push_back( (Persistence_heat_maps*)maps_[i] );
	  }
	  this->compute_mean(maps);
  }

  void compute_median_interface(const std::vector<Persistence_heat_maps_interface*>& maps_)
  {
	  std::vector<Persistence_heat_maps*> maps;
	  maps.reserve( maps_.size() );
	  for ( size_t i = 0 ; i != maps_.size() ; ++i )
	  {
		  maps.push_back( (Persistence_heat_maps*)maps_[i] );
	  }
	  this->compute_median(maps);
  }
  
  void compute_percentage_of_active_interface(const std::vector<Persistence_heat_maps_interface*>& maps_, size_t cutoff = 1)
  {
	  std::vector<Persistence_heat_maps*> maps;
	  maps.reserve( maps_.size() );
	  for ( size_t i = 0 ; i != maps_.size() ; ++i )
	  {
		  maps.push_back( (Persistence_heat_maps*)maps_[i] );
	  }
	  this->compute_percentage_of_active(maps,cutoff);
  }
  
  void compute_average_interface(const std::vector<Persistence_heat_maps_interface*>& to_average)
  {
	  std::vector<Persistence_heat_maps*> maps;
	  maps.reserve( to_average.size() );
	  for ( size_t i = 0 ; i != to_average.size() ; ++i )
	  {
		  maps.push_back( (Persistence_heat_maps*)to_average[i] );
	  }
	  this->compute_average( maps );
  }  
  
  bool compare( const Persistence_heat_maps_interface& second )const 
  {
	  return *this == second;
  }
  
  static bool erase_below_diagonal;
};

bool Persistence_heat_maps_interface::erase_below_diagonal = false;


}  // namespace Persistence_representations
}  // namespace Gudhi

#endif  // PERSISTENCE_HEAT_MAPS_INTERFACE_H_
