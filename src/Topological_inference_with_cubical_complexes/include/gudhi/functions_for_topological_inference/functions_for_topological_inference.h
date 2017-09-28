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
//#ifdef GUDHI_USE_CGAL
	#include <gudhi/Kd_tree_search.h>
	#include <CGAL/Epick_d.h>
//#endif

namespace Gudhi {

namespace Topological_inference_with_cubical_complexes {



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
	kernels_centerd_in_point_cloud( const Kernel& kernel_, const std::vector< std::vector<double> >& point_cloud_ ):point_cloud(point_cloud_),kernel(kernel_){}
	
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
	Kernel kernel;
	std::vector< std::vector<double> > point_cloud;	
};//kernels_centerd_in_point_cloud



/**
 * This class take a point cloud and implement a function which is a sum of distances
 * from a given point to the points in point cloud. 
**/  
class Sum_of_distances_from_points
{
public:
	Sum_of_distances_from_points( const std::vector< std::vector<double> >& point_cloud_ , double shift_ = 1 )
	{		
		this->shift = shift_;
		this->point_cloud = point_cloud_;	
	}
    double operator() ( const std::vector<double>& point )const
    {
		bool dbg = false;
        double result = 0;                   
        for ( size_t i = 0 ; i != this->point_cloud.size() ; ++i )
        {		
			double distance_from_this_point = 0;	
			for ( size_t aa = 0  ; aa != this->point_cloud[i].size() ; ++aa )
			{
				distance_from_this_point += (this->point_cloud[i][aa] - point[aa])*(this->point_cloud[i][aa] - point[aa]);
			}			
			result += distance_from_this_point;
			if ( dbg )
			{
				std::cout << "The distance is : " << distance_from_this_point << " and 1/(distance + shift) = " << 1.0/(sqrt( distance_from_this_point )+this->shift) << std::endl;
				getchar();
			}
		}
		return result;
    }
private:
	std::vector< std::vector<double> > point_cloud;
	double shift;
};



//template < int dimension >
class Distance_to_k_th_closest_point
{
public:
	//#ifdef GUDHI_USE_CGAL
		//typedef CGAL::Epick_d<CGAL::Dimension_tag< dimension > > K;
		typedef CGAL::Epick_d< CGAL::Dynamic_dimension_tag > K;
		typedef typename K::Point_d Point;
		typedef std::vector<Point> Points;
		typedef Gudhi::spatial_searching::Kd_tree_search<K, Points> k_d_tree;							
	//#endif







	Distance_to_k_th_closest_point( const std::vector< std::vector<double> >& point_cloud_ )
	{
		//check if all the points have the same dimension:
		if ( point_cloud_.empty() )return;
		for ( size_t i = 0 ; i != point_cloud_.size() ; ++i )
		{
			if ( point_cloud_[i].size() != dimension )
			{
				std::cerr << "Point cloud containing points of different dimension in Distance_to_closest_point constructor.\n";
				throw "Point cloud containing points of different dimension in Distance_to_closest_point constructor.\n";
			}
		}						
		this->point_cloud = point_cloud_;	
		//#ifdef GUDHI_USE_CGAL
			Points points;
			for (size_t i = 0; i != point_cloud_.size() ; ++i)
				points.push_back( Point( point_cloud_[i] ) );
			this->points_ds(points);
		//#endif	
	}
    double operator()( const std::vector<double>& point , unsigned k )const
    {
		//#ifdef GUDHI_USE_CGAL
			auto knn_range = this->points_ds.query_k_nearest_neighbors(point, k, true);
			return knn_range[k-1]; //test it if it should not be k!!!
		/*#else
			//this is brutal, n^2 version, in case we do not have a k-d tree from CGAL. 
			double result = std::numeric_limits<double>::max();                   
			for ( size_t i = 0 ; i != this->point_cloud.size() ; ++i )
			{					
				double distance_from_this_point = 0;	
				for ( size_t aa = 0  ; aa != this->point_cloud[i].size() ; ++aa )
				{
					distance_from_this_point += (this->point_cloud[i][aa] - point[aa])*(this->point_cloud[i][aa] - point[aa]);
				}
				if ( result > distance_from_this_point )result = distance_from_this_point;
			}
			return result;				  
		#endif*/	
    }
private:
	std::vector< std::vector<double> > point_cloud;
	//#ifdef GUDHI_USE_CGAL
		k_d_tree points_ds;
	//#endif
};

//****************************************************************
//Sample kernels to be used in kernels_centerd_in_point_cloud
//****************************************************************
double Euclidan_distance_squared( const std::vector<double>& point1 , const std::vector<double>& point2 )
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

double Manhattan_distance_squared( const std::vector<double>& point1 , const std::vector<double>& point2 )
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
//to be continued. 

}  // namespace Topological_inference_with_cubical_complexes
}  // namespace Gudhi

#endif
