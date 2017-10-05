/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Pawel Dlotko
 *
 *    Copyright (C) 2017 Swansea University UK
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
#ifndef FUNCTIONS_FOR_TOPOLOGICAL_INFERENCE_H
#define FUNCTIONS_FOR_TOPOLOGICAL_INFERENCE_H


//#define GUDHI_USE_CGAL


#include <cmath>
#include <vector>
#include <algorithm>
#include <limits>
#include <iterator>
#include <gudhi/Bitmap_cubical_complex/counter.h>
#ifdef GUDHI_USE_CGAL
	#include <gudhi/Kd_tree_search.h>
	#include <CGAL/Epick_d.h>
#endif

namespace Gudhi {

namespace Topological_inference_with_cubical_complexes {
	
	
	

//****************************************************************
//Sample kernels to be used in kernels_centerd_in_point_cloud
//****************************************************************
class Euclidan_distance_squared
{
public:
	Euclidan_distance_squared(){}
	double operator()( const std::vector<double>& point1 , const std::vector<double>& point2 )
	{
		if ( point1.size() != point2.size() )
		{
			std::cerr << "Incompatible dimensions of points.\n";
			throw "Incompatible dimensions of points.\n";
		}
		double result = 0;
		for ( size_t i = 0 ; i != point1.size() ; ++i )
		{
			result += (point1[i]-point2[i]) * (point1[i]-point2[i]);
		}
		return result;
	}
};


class Manhattan_distance
{
public:
	Manhattan_distance(){}
	double operator()( const std::vector<double>& point1 , const std::vector<double>& point2 )
	{
		if ( point1.size() != point2.size() )
		{
			std::cerr << "Incompatible dimensions of points.\n";
			throw "Incompatible dimensions of points.\n";
		}
		double result = 0;
		for ( size_t i = 0 ; i != point1.size() ; ++i )
		{
			result += fabs(point1[i]-point2[i]);
		}
		return result;
	}
};


class Max_norm_distance
{
public:
	Max_norm_distance(){}
	double operator()( const std::vector<double>& point1 , const std::vector<double>& point2 )
	{
		if ( point1.size() != point2.size() )
		{
			std::cerr << "Incompatible dimensions of points.\n";
			throw "Incompatible dimensions of points.\n";
		}
		double result = -std::numeric_limits<double>::max();
		for ( size_t i = 0 ; i != point1.size() ; ++i )
		{
			if ( result < fabs(point1[i]-point2[i]) )result = fabs(point1[i]-point2[i]);
		}
		return result;
	}
};

template <typename standard_distance>
class periodic_domain_distance
{
public:
	periodic_domain_distance( const std::vector< std::pair< double , double > >& coordinates_of_grid , standard_distance& dist_ ):dist( dist_ )
	{	
		bool dbg = false;
		this->translation_vector.reserve( coordinates_of_grid.size() );
		for ( size_t i = 0 ; i != coordinates_of_grid.size() ; ++i )
		{			
			this->translation_vector.push_back( coordinates_of_grid[i].second - coordinates_of_grid[i].first );
		}
		if ( dbg )
		{
			std::cout << "Here is the translation vector: \n";
			for ( size_t i = 0 ; i != this->translation_vector.size() ; ++i )
			{
				std::cout << this->translation_vector[i] << " ";
			}
		}
	}
	double operator()( const std::vector<double>& point1 , const std::vector<double>& point2 )
	{		
		bool dbg = false;
		if ( point1.size() != point2.size() )
		{
			std::cerr << "Incompatible dimensions of points.\n";
			throw "Incompatible dimensions of points.\n";
		}
		//transate in all possible directions and compute the minimal between translated
		//points by using this->dist. 
		std::vector< unsigned > beginn_counter(point1.size(),0);
		std::vector< unsigned > end_counter(point1.size(),2);
		Gudhi::cubical_complex::counter count(beginn_counter,end_counter);
		
		//structure to store shifted versions of point2:
		std::vector<double> shifted_point2(point2.size(),0);
		double new_distance;
		
		double minimal_distance = std::numeric_limits<double>::max();
		while ( true )
		{				
			if ( dbg )std::cerr << "Here is the current counter : " << count << std::endl << "And here is the shifted point " << std::endl;			
			//shift point2
			for ( size_t i = 0 ; i != point2.size() ; ++i )
			{
				shifted_point2[i] = point2[i]+((int)count[i]-1)*this->translation_vector[i];
				if ( dbg ) std::cerr << shifted_point2[i] << " ";
			}			
			if ( dbg ) std::cerr << std::endl;
			
			//compute distance between point1 and shited version of point2;
			new_distance = this->dist(point1,shifted_point2);
			if ( new_distance < minimal_distance )
			{
				minimal_distance = new_distance;
			}
			if ( dbg )std::cout << "new_distance : " << new_distance << std::endl;
			
			bool can_be_incremented = count.increment();
			if ( !can_be_incremented )break;			
		}
		if ( dbg )std::cout << "We are done, the minimal distance is : " << minimal_distance << std::endl;
		return minimal_distance;
	}
protected:
	//in this very simple case, when all the boxes are squares, translation vectors are zero
	//in all the coordinates except from the i-th one for i-th translation vector. 
	std::vector< double > translation_vector;
	standard_distance dist;
};
//to be continued. 
//****************************************************************
//End sample kernels to be used in kernels_centerd_in_point_cloud
//****************************************************************




/**
 * This class take a point cloud and compute a function which is a sum of kernels centered
 * at the points of the point cloud. A few example of kernel functions are available in this 
 * file. 
**/ 
template <typename Kernel>
class kernels_centerd_in_point_cloud
{
public:
   /**
    * A constructor of a kernels_centerd_in_point_cloud class taking as a parameter the point cloud and the kernel. Then this kernel is placed on the top of
    * every point from the point cloud. The clas compute the value of a function being sum of all the distributions 
    **/ 
	kernels_centerd_in_point_cloud( Kernel& kernel_, const std::vector< std::vector<double> >& point_cloud_ ):kernel(kernel_),point_cloud(point_cloud_){}
	
	/**
	 * A function to compute the value of the kernels_centerd_in_point_cloud at a given point.
	**/ 
    double operator() ( const std::vector<double>& point )const
    {
        double result = 0;                   
        for ( size_t i = 0 ; i != this->point_cloud.size() ; ++i )
        {			
			result += this->kernel( this->point_cloud[i] , point );
		}
		return result;
    }
private:
	Kernel& kernel;
	std::vector< std::vector<double> > point_cloud;	
};//kernels_centerd_in_point_cloud



/**
 * This class take a point cloud and implement a function which is a sum of distances
 * from a given point to the points in point cloud. 
**/  
class Sum_of_distances_from_points
{
public:
	Sum_of_distances_from_points( const std::vector< std::vector<double> >& point_cloud_ )
	{		
		this->point_cloud = point_cloud_;	
	}
    double operator() ( const std::vector<double>& point )const
    {
        double result = 0;                   
        for ( size_t i = 0 ; i != this->point_cloud.size() ; ++i )
        {		
			double distance_from_this_point = 0;	
			for ( size_t aa = 0  ; aa != this->point_cloud[i].size() ; ++aa )
			{
				distance_from_this_point += (this->point_cloud[i][aa] - point[aa])*(this->point_cloud[i][aa] - point[aa]);
			}			
			result += distance_from_this_point;		
		}
		return result;
    }
private:
	std::vector< std::vector<double> > point_cloud;
};


template < typename distance >
class Distance_to_k_th_closest_point
{
public:	
	Distance_to_k_th_closest_point( const std::vector< std::vector<double> >& point_cloud_ , distance& dist_ , unsigned k_ ):
	dist(dist_),k(k_)
	#ifdef GUDHI_USE_CGAL
		,points_ds( this->convert_points(point_cloud_) )
	#endif	
	{
		//check if all the points have the same dimension:
		if ( point_cloud_.empty() )return;
		for ( size_t i = 0 ; i != point_cloud_.size() ; ++i )
		{
			if ( point_cloud_[i].size() != point_cloud_[0].size() )
			{
				std::cerr << "Point cloud containing points of different dimension in Distance_to_closest_point constructor.\n";
				throw "Point cloud containing points of different dimension in Distance_to_closest_point constructor.\n";
			}
		}						
		this->point_cloud = point_cloud_;				
	}
    double operator()( const std::vector<double>& point )const
    {
		#ifdef GUDHI_USE_CGAL		
		std::cerr << "this->point_cloud.size() : " << this->point_cloud.size() << std::endl;
		std::cerr << "point.size()  : " << point.size() << std::endl;
		std::cerr << "this->k :  " << this->k << std::endl;
		
		std::cout << "depth : " << this->points_ds.tree_depth() << "\n";
		
			//I have tried to use kere std::advance, but it got too messy with the constant interators of unknown type. So I did it by a bruteforce. 
			auto knn_range = this->points_ds.query_k_nearest_neighbors(
			CGAL::Epick_d< CGAL::Dynamic_dimension_tag >::Point_d( point.begin() , point.end() ), this->k+1, true);										
			size_t counter = 0;
			double dist_ = 0;
			for (auto const& nghb : knn_range)
			{
				
				std::cout << nghb.first << " (sq. dist. = " << nghb.second << ")\n"; 
				++counter;
				if ( counter == k )
				{
					dist_ = nghb.second;
					break;
				}				
			}
			return dist_;
		 #else		 
			//this is brutal, version, in case we do not have a k-d tree from CGAL. 
			bool dbg = false;
			
			std::vector< double > heap( k , std::numeric_limits<double>::max() );
			std::make_heap (heap.begin(),heap.end());
			//now compute the distances between point and any other point, and put it into heap:			
			for ( size_t i = 0 ; i != this->point_cloud.size() ; ++i )
			{
				//compute distance
				double distance_ = this->dist( this->point_cloud[i] , point );
				//check if it is smaller than the largest element in the heap:
				
				if ( dbg )
				{
					std::cout << "distance_ : " << distance_ << std::endl;
					std::cout << "heap.front() : " << heap.front() << std::endl;
					std::cout << "Here is the whole heap: \n";
					for ( size_t i = 0 ; i != heap.size() ; ++i )
					{
						std::cout << heap[i] << " ";
					}
					getchar();
				}
				
				if ( distance_ < heap.front() )
				{				
					if ( dbg ){std::cerr << "Adding new element to the heap \n";}
						
					//remove the largest element from the heap.
					std::pop_heap (heap.begin(),heap.end()); 
					heap.pop_back();
					//put distance to the heap.
					heap.push_back( distance_ );				   
					//update the heap.
					std::push_heap (heap.begin(),heap.end());
				}			
			}
			
			if ( dbg )
			{
				for ( size_t i = 0 ; i != heap.size() ; ++i )
				{
					std::cout << heap[i] << " ";
				}
				getchar();
			}
			
			
			//and now we have a heap of smallest distances from point to the point in the point cloud. 
			//we will now find the minimal one:
			double result = -std::numeric_limits<double>::max();
			for ( size_t i = 0 ; i != heap.size() ; ++i )
			{
				if ( heap[i] > result )result = heap[i];
			}  
			return result;		
		#endif
    }
private:

	#ifdef GUDHI_USE_CGAL
		std::vector< CGAL::Epick_d< CGAL::Dynamic_dimension_tag >::Point_d > convert_points( const std::vector< std::vector<double> >& point_cloud_ )
		{
			std::vector< CGAL::Epick_d< CGAL::Dynamic_dimension_tag >::Point_d > points;
			points.reserve( point_cloud_.size() );
			for (size_t i = 0; i != point_cloud_.size() ; ++i)
			{				
				points.push_back( CGAL::Epick_d< CGAL::Dynamic_dimension_tag >::Point_d( point_cloud_[i].begin() , point_cloud_[i].end() ) );				
			}	
			
			for ( auto c : points[2] )
				std::cout << c << std::endl;
			
			std::cerr << "points.size() : "  << points.size() << std::endl;
						
			return points;
		}
	#endif

	std::vector< std::vector<double> > point_cloud;
	distance& dist;
	unsigned k;
	#ifdef GUDHI_USE_CGAL
		Gudhi::spatial_searching::Kd_tree_search< CGAL::Epick_d< CGAL::Dynamic_dimension_tag > , 
		std::vector< CGAL::Epick_d< CGAL::Dynamic_dimension_tag >::Point_d > > points_ds;		
	#endif
};



}  // namespace Topological_inference_with_cubical_complexes
}  // namespace Gudhi

#endif
