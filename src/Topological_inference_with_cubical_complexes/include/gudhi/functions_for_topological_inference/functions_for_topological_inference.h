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

#include <cmath>
#include <vector>
#include <algorithm>
#include <limits>
#include <iterator>
#include <gudhi/Bitmap_cubical_complex/counter.h>
#ifdef GUDHI_USE_CGAL
	#include <gudhi/Kd_tree_search.h>	
	#include <gudhi/Periodic_kd_tree_search.h>
	#include <CGAL/Epick_d.h>
#endif

namespace Gudhi {

namespace Topological_inference_with_cubical_complexes {
	
	
	

//****************************************************************
//Sample kernels to be used in kernels_centerd_in_point_cloud
//****************************************************************
/**
 * A class that compute square of Euclidean distance between points.
**/ 
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

/**
 * A class that compute square of Manhattan distance between points.
**/ 
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


/**
 * A class that compute square of Max norm distance between points.
**/ 
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
		double result = -std::numeric_limits<double>::infinity();
		for ( size_t i = 0 ; i != point1.size() ; ++i )
		{
			if ( result < fabs(point1[i]-point2[i]) )result = fabs(point1[i]-point2[i]);
		}
		return result;
	}
};

/**
 * A class that any type of distance in a periodic domain. The template parameter 
 * of that class is a distace in non periodic domain that is used to compute distance
 * in the periodic domain. See the description of constructor of this class
 * for further details.
**/ 
template <typename standard_distance>
class periodic_domain_distance
{
public:
	/**
	 * This is a constructor of a periodic_domain_distance class. It take as a parameter the 
	 * std::vector< std::pair< double , double > >& coordinates_of_grid describing coodinates of
	 * the grid. See coordinates_of_grid in the Topological_inference_with_cubical_complexes class
	 * for details. The second parameter of the constructor is a distance that is used to compute
	 * distance in periodic domain. 
	**/ 
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
	/**
	* Computations of distance between point1 and point2 on the predefined periodic domain.
	* The code keep point1 fixed, and copy point 2 in a covering space in all possible directions.
	* For each of them, standard distance in the covering space is computed, and the minimal one
	* is returned. 
	**/ 
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
		
		double minimal_distance = std::numeric_limits<double>::infinity();
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
//to be continued when new inspiration comes. 
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


/**
 * A class that take a point cloud P, and compute distance from any
 * given point to the k-th nearest neighbour in P. 
 * This is a general class that can be used for any distance, and both periodic 
 * and non-periodic domains. For more efficeint version (usng gudhi spatial searchig
 * requiering CGAL) designed only for Euclidean distance in non-periodic domain, 
 * please use Distance_to_k_th_closest_point_k_d_tree.
**/ 
template < typename distance >
class Distance_to_k_th_closest_point
{
public:	
	Distance_to_k_th_closest_point( const std::vector< std::vector<double> >& point_cloud_ , distance& dist_ , unsigned k_ ):
	dist(dist_),k(k_)	
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
		//this is brutal, version, in case we do not have a k-d tree from CGAL. 
		bool dbg = false;
		
		std::vector< double > heap( k , std::numeric_limits<double>::infinity() );
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
		double result = -std::numeric_limits<double>::infinity();
		for ( size_t i = 0 ; i != heap.size() ; ++i )
		{
			if ( heap[i] > result )result = heap[i];
		}  
		return result;	
	
    }
private:
	std::vector< std::vector<double> > point_cloud;
	distance& dist;
	unsigned k;	
};



#ifdef GUDHI_USE_CGAL
/**
 * A class that take a point cloud P, and compute Euclidean distance from any
 * given point to the k-th nearest neighbour in P. 
 * CGAL is required to compile it. 
**/ 
class Distance_to_k_th_closest_point_k_d_tree
{
public:	
	Distance_to_k_th_closest_point_k_d_tree( const std::vector< std::vector<double> >& point_cloud_ , unsigned k_ ):
	points(convert_points(point_cloud_)),points_ds(this->points),k(k_)
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
		auto pt = CGAL::Epick_d< CGAL::Dynamic_dimension_tag >::Point_d( point.begin() , point.end() );				
		//I have tried to use kere std::advance, but it got too messy with the constant interators of unknown type. So I did it by a bruteforce. 
		double dist_ = std::numeric_limits< double >::infinity();		
		auto knn_range = this->points_ds.query_k_nearest_neighbors(	pt, this->k+1, true);					
		size_t counter = 0;			
		for (auto const& nghb : knn_range)
		{						
			++counter;
			if ( counter == k )
			{
				dist_ = nghb.second;
				break;
			}				
		}					
		return dist_;	
    }
private:
	std::vector< CGAL::Epick_d< CGAL::Dynamic_dimension_tag >::Point_d >	
	convert_points( const std::vector< std::vector<double> >& point_cloud_ )
	{		
		std::vector< CGAL::Epick_d< CGAL::Dynamic_dimension_tag >::Point_d > points_;
		points_.reserve( point_cloud_.size() );
		for (size_t i = 0; i != point_cloud_.size() ; ++i)
		{							
			points_.push_back( CGAL::Epick_d< CGAL::Dynamic_dimension_tag >::Point_d( point_cloud_[i].begin() , point_cloud_[i].end() ) );											
		}
		return points_;		    
	}

	std::vector< std::vector<double> > point_cloud;
	std::vector< CGAL::Epick_d< CGAL::Dynamic_dimension_tag >::Point_d > points;
	
	Gudhi::spatial_searching::Kd_tree_search< CGAL::Epick_d< CGAL::Dynamic_dimension_tag > , 
	std::vector< CGAL::Epick_d< CGAL::Dynamic_dimension_tag >::Point_d > > points_ds;			
	
	unsigned k;		
};
#endif




#ifdef GUDHI_USE_CGAL
/**
 * A class that implements distace in periodic domain using periodic
 * kd-trees from spatial seraching procedure.
**/
class Distance_to_k_th_closest_point_periodic_k_d_tree
{
public:
	/**
	 * This is a constructor of a Distance_to_k_th_closest_point_periodic_k_d_tree class. 
	 * It take as a parameter the  std::vector< std::pair< double , double > >& coordinates_of_grid 
	 * describing coodinates of the periodic cuboid. See coordinates_of_grid in the 
	 * Topological_inference_with_cubical_complexes class for details. 
	 * The second parameter is the point cloud in which k-th nearest neighbor will be searched for.
	 * The last parameter is the number k, such that the distance to the k-th nearest neighbor is
	 * found by the procedure.
	**/ 
	Distance_to_k_th_closest_point_periodic_k_d_tree( 
							const std::vector< std::vector<double> >& point_cloud_,
							const std::vector< std::pair< double , double > >& coordinates_of_grid,							
							unsigned k_
	                        ):
	lower_left(compute_coordinates_of_ISO_box(coordinates_of_grid,true)),
	top_right(compute_coordinates_of_ISO_box(coordinates_of_grid,false)),
	points(this->convert_points(point_cloud_)),
	points_ds( this->points,K::Iso_box_d( this->lower_left,
	                                      this->top_right ) ),
	                                      k(k_){}
	    
	
	/**
	* Computations of distance between point1 and point2 on the predefined periodic domain.
	* The code keep point1 fixed, and copy point 2 in a covering space in all possible directions.
	* For each of them, standard distance in the covering space is computed, and the minimal one
	* is returned. 
	**/ 
	double operator()( const std::vector<double>& point )
	{		
		auto pt = CGAL::Epick_d< CGAL::Dynamic_dimension_tag >::Point_d( point.begin() , point.end() );							
		//I have tried to use kere std::advance, but it got too messy with the constant interators of unknown type. So I did it by a bruteforce. 
		double dist_ = std::numeric_limits< double >::infinity();		
		auto knn_range = this->points_ds.query_k_nearest_neighbors(	pt, this->k+1, true);		
		size_t counter = 0;		
		for (auto const& nghb : knn_range)
		{						
			++counter;
			if ( counter == k )
			{
				dist_ = nghb.second;
				break;
			}				
		}
		return dist_;	
	}
	
    typedef CGAL::Epick_d< CGAL::Dynamic_dimension_tag > K;
    typedef typename K::Point_d Point;
    typedef Gudhi::spatial_searching::Periodic_kd_tree_search<K, std::vector<Point> > Points_ds;    
  
protected:
	
	std::vector< Point >	
	convert_points( const std::vector< std::vector<double> >& point_cloud_ )
	{		
		std::vector< CGAL::Epick_d< CGAL::Dynamic_dimension_tag >::Point_d > points_;
		points_.reserve( point_cloud_.size() );
		for (size_t i = 0; i != point_cloud_.size() ; ++i)
		{							
			points_.push_back( Point( point_cloud_[i].begin() , point_cloud_[i].end() ) );											
		}
		return points_;		    
	}//convert_points
	
	Point compute_coordinates_of_ISO_box( const std::vector< std::pair< double , double > >& coordinates_of_grid , double shall_I_take_first_component )
	{		
		std::vector<double> vect( coordinates_of_grid.size() , 0 );	
		for ( size_t i = 0 ; i != coordinates_of_grid.size() ; ++i )
		{
			if ( shall_I_take_first_component )
			{
				vect[i] = coordinates_of_grid[i].first;
			}
			else
			{
				vect[i] = coordinates_of_grid[i].second;
			}
		}
		Point p( vect.begin() , vect.end() );
		return p;
	}//compute_coordinates_of_ISO_box
		
	Point lower_left;		
	Point top_right;
	std::vector< Point > points;		
	Points_ds points_ds;
	unsigned k;			
};
#endif




}  // namespace Topological_inference_with_cubical_complexes
}  // namespace Gudhi

#endif
