/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Pawel Dlotko
 *
 *    Copyright (C) 2018  Swansea University
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

#ifndef PSSK_INTERFACE_H_
#define PSSK_INTERFACE__H_


#include <gudhi/PSSK.h>


namespace Gudhi {
namespace Persistence_representations {

/**
* This is a version of a representation presented in https://arxiv.org/abs/1412.6821
* In that paper the authors are using the representation just to compute kernel. Over here, we extend the usability by
*far.
* Note that the version presented here is not exact, since we are discretizing the kernel.
* The only difference with respect to the original class is the method of creation. We have full (square) image, and for
*every point (p,q), we add a kernel at (p,q) and the negative kernel
* at (q,p)
**/

class PSSK_interface : public PSSK {
 public:
  PSSK_interface(){}
  				

  PSSK_interface(const std::vector<std::pair<double, double> >& interval,
       size_t how_many_pixels_raidus_of_Gausian_kernel, size_t number_of_pixels = 1000,
       double min_ = -1, double max_ = -1)
      :
      PSSK(interval,create_Gaussian_filter(how_many_pixels_raidus_of_Gausian_kernel, 1),
           number_of_pixels,((min_ == max_) ?  std::numeric_limits<double>::max() : min_),
           ((min_ == max_) ?  -1 : max_)){}

  PSSK_interface(const char* filename, size_t how_many_pixels_raidus_of_Gausian_kernel,
       size_t number_of_pixels = 1000, double min_ = -1, double max_ = -1,
       int dimension = -1)
      :PSSK(filename,create_Gaussian_filter(how_many_pixels_raidus_of_Gausian_kernel, 1),
            number_of_pixels,((min_ == max_) ?  std::numeric_limits<double>::max() : min_),
            ((min_ == max_) ?  -1 : max_),dimension){}
            
  
  
  static PSSK_interface* construct_from_file(  const char* filename, size_t how_many_pixels_raidus_of_Gausian_kernel,
																size_t number_of_pixels, double min_ = 0,
																double max_ = 0, int dimensions = -1 )
  {
      PSSK_interface* result = new PSSK_interface(filename,how_many_pixels_raidus_of_Gausian_kernel,number_of_pixels,min_,max_,dimensions);
	  return result;  
  }
  
  
  static PSSK_interface* construct_from_vector_of_pairs( const std::vector<std::pair<double, double> >& interval, size_t how_many_pixels_raidus_of_Gausian_kernel, size_t number_of_pixels, double min_ = 0, double max_ = 0 )
  {
      PSSK_interface* result = new PSSK_interface(interval, how_many_pixels_raidus_of_Gausian_kernel, number_of_pixels, min_, max_);
	  return result;  
  }    
  
  
  void compute_mean_interface(const std::vector<PSSK_interface*>& maps_)
  {
	  std::vector<Persistence_heat_maps*> maps;
	  maps.reserve( maps_.size() );
	  for ( size_t i = 0 ; i != maps_.size() ; ++i )
	  {
		  maps.push_back( (Persistence_heat_maps*)maps_[i] );
	  }
	  this->compute_mean(maps);
  }

  void compute_median_interface(const std::vector<PSSK_interface*>& maps_)
  {
	  std::vector<Persistence_heat_maps*> maps;
	  maps.reserve( maps_.size() );
	  for ( size_t i = 0 ; i != maps_.size() ; ++i )
	  {
		  maps.push_back( (Persistence_heat_maps*)maps_[i] );
	  }
	  this->compute_median(maps);
  }
  
  void compute_percentage_of_active_interface(const std::vector<PSSK_interface*>& maps_, size_t cutoff = 1)
  {
	  std::vector<Persistence_heat_maps*> maps;
	  maps.reserve( maps_.size() );
	  for ( size_t i = 0 ; i != maps_.size() ; ++i )
	  {
		  maps.push_back( (Persistence_heat_maps*)maps_[i] );
	  }
	  this->compute_percentage_of_active(maps,cutoff);
  }
  
  void compute_average_interface(const std::vector<PSSK_interface*>& to_average)
  {
	  std::vector<Persistence_heat_maps*> maps;
	  maps.reserve( to_average.size() );
	  for ( size_t i = 0 ; i != to_average.size() ; ++i )
	  {
		  maps.push_back( (Persistence_heat_maps*)to_average[i] );
	  }
	  this->compute_average( maps );
  }        
  
  bool compare( const PSSK_interface& second )const
  {
	  return ( *this == second );
  }
 
};

}  // namespace Persistence_representations
}  // namespace Gudhi

#endif  // PSSK_INTERFACE__H_
