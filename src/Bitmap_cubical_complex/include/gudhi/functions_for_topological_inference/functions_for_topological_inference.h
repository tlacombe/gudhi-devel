//for normal_distribution
#include <cmath>


//This extra functions are placed in Gudhi stat, and when this is all merged to the same branch, they should be taken from there!
/**
 * This is a simple procedure to create n by n (or 2*pixel_radius times 2*pixel_radius.
**/ 
std::vector< std::vector<double> > create_Gaussian_filter( size_t pixel_radius , double sigma )
{
	bool dbg = false;
	//we are computing the kernel mask to 2 standard deviations away from the center. We discretize it in a grid of a size 2*pixel_radius times 2*pixel_radius.
	
    double r = 0;
    double sigma_sqr = sigma * sigma;
 
    // sum is for normalization
    double sum = 0;
    
    //initialization of a kernel:
    std::vector< std::vector<double> > kernel( 2*pixel_radius +1 );
    for ( size_t i = 0 ; i != kernel.size() ; ++i )
    {
		std::vector<double> v( 2*pixel_radius +1 , 0 );
		kernel[i] = v;
	}
	
	if ( dbg )
	{
		std::cerr << "Kernel initalize \n";	
		cerr << "pixel_radius : " << pixel_radius << endl; 
		std::cerr << "kernel.size() : " << kernel.size() << endl;
		getchar();
	}
 
    for (int x = -pixel_radius; x <= (int)pixel_radius; x++)
    {
        for(int y = -pixel_radius; y <= (int)pixel_radius; y++)
        {			
			double real_x = 2*sigma*x/pixel_radius;
			double real_y = 2*sigma*y/pixel_radius;						
            r = sqrt(real_x*real_x + real_y*real_y);						
            kernel[x + pixel_radius][y + pixel_radius] = (exp(-(r*r)/sigma_sqr))/(3.141592 * sigma_sqr);
            sum += kernel[x + pixel_radius][y + pixel_radius];            
        }
    }      
 
    // normalize the kernel
    for( size_t i = 0; i != kernel.size() ; ++i)
    {
        for( size_t j = 0; j != kernel[i].size() ; ++j)
        {
            kernel[i][j] /= sum;
		}
            
    }
    
    if ( dbg )
    {
		std::cerr << "Here is the kernel : \n";
		for( size_t i = 0; i != kernel.size() ; ++i)
		{
			for( size_t j = 0; j != kernel[i].size() ; ++j)
			{
				std::cerr << kernel[i][j] << " ";
			}
			std::cerr << std::endl;
		}
	}
    return kernel;
}






//we assume that there are no spaces at the end of lines.
template<typename T>
std::vector< std::vector<T> > read_points_from_file( const char* filename , bool compute_ranges = true )
{
	bool dbg = false;
	std::vector< std::vector<T> > result;
	
	std::ifstream in;
	in.open( filename );
	
	
	
	std::string line;
	while (std::getline(in, line))
	{
		std::istringstream iss(line);
		std::vector<T> point;
		size_t i = 0;
		while ( true )
		{
			T a;
			iss >> a;			
			point.push_back( a );		
			if (dbg)std::cout << a << " ";
			++i;
			if ( !iss.good() )break; 			
		}
		if (dbg)
		{
			std::cout << std::endl;
			//getchar();
		}
		result.push_back( point );
	}
	
	if ( compute_ranges )
	{
		std::vector< T > mins( result[0].size() , INT_MAX );
		std::vector< T > max( result[0].size() , INT_MIN );
		for ( size_t i = 0 ; i != result.size() ; ++i )
		{
			for ( size_t pt = 0 ; pt != result[i].size() ; ++pt )
			{
				if ( mins[pt] > result[i][pt] )mins[pt] = result[i][pt];
				if ( max[pt] < result[i][pt] )max[pt] = result[i][pt];
			}
		}
		
		std::cout << "Ranges: " << std::endl;
		for ( size_t i = 0 ; i != mins.size() ; ++i )
		{
			std::cout << mins[i] << " "  << max[i] << std::endl;
		}
	}
	
	in.close();
	return result;
}

//***************END OF FUNCTIONS THAT SHOULD BE COMMON FOR THE GUDHI STAT.




//*****************************************************************************************//
//*****************************************************************************************//
//*****************************************************************************************//
//*****************************************************************************************//
//*********************Functions for topological inference*********************************//
//*****************************************************************************************//
//*****************************************************************************************//
//*****************************************************************************************//
//*****************************************************************************************//




class Gaussian_kernels_centerd_in_point_cloud
{
public:
	Gaussian_kernels_centerd_in_point_cloud( const std::vector< std::vector<double> >& point_cloud_ , double standard_deviation_ )
	{
		this->standard_deviation = standard_deviation_;
		this->point_cloud = point_cloud_;
	}
    double operator() ( const std::vector<double>& point )const
    {
        double result = 0;                   
        for ( size_t i = 0 ; i != this->point_cloud.size() ; ++i )
        {
			//compute value of a Gaussian kernel centered in this->point_cloud[i] on a point 
			for ( size_t aa = 0  ; aa != this->point_cloud[i].size() ; ++aa )
			{
				 //TODO
				//result += ( 1 / ( s * sqrt(2*M_PI) ) ) * exp( -0.5 * pow( (x-)/s, 2.0 ) );
			}
		}
		return result;
    }
private:
	std::vector< std::vector<double> > point_cloud;
	double standard_deviation;
};


class Sum_of_inverse_of_distances_from_points
{
public:
	Sum_of_inverse_of_distances_from_points( const std::vector< std::vector<double> >& point_cloud_ , double shift_ = 1 )
	{
		bool dbg = false;
		this->shift = shift_;
		this->point_cloud = point_cloud_;
		if ( dbg )
		{
			std::cout << "Point cloud \n";
			for ( size_t i = 0 ; i != this->point_cloud.size() ; ++i )
			{
				for ( size_t j = 0 ; j != this->point_cloud[i].size() ; ++j )
				{
					std::cerr << this->point_cloud[i][j] << " ";
				}	
				std::cerr << std::endl;		
			}
		}
	}
    double operator() ( const std::vector<double>& point )const
    {
		bool dbg = false;
        double result = 0;                   
        for ( size_t i = 0 ; i != this->point_cloud.size() ; ++i )
        {		
			if ( dbg )
			{
				std::cout << "Computing a distance of a point : " << std::endl;
				for ( size_t yy = 0 ; yy != point.size() ; ++yy )
				{
					std::cout << point[yy] << " ";
				}
				std::cout << std::endl;
				std::cout << "From a point : " << std::endl;
				for ( size_t yy = 0 ; yy != this->point_cloud[i].size() ; ++yy )
				{
					std::cout << this->point_cloud[i][yy] << " ";
				}
				std::cout << std::endl;
			}
			double distance_from_this_point = 0;	
			for ( size_t aa = 0  ; aa != this->point_cloud[i].size() ; ++aa )
			{
				distance_from_this_point += (this->point_cloud[i][aa] - point[aa])*(this->point_cloud[i][aa] - point[aa]);
			}
			//result += 1.0/(sqrt( distance_from_this_point )+this->shift);
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





class Distance_to_closest_point
{
public:
	Distance_to_closest_point( const std::vector< std::vector<double> >& point_cloud_ )
	{
		bool dbg = false;
		this->point_cloud = point_cloud_;
		if ( dbg )
		{
			std::cout << "Point cloud \n";
			for ( size_t i = 0 ; i != this->point_cloud.size() ; ++i )
			{
				for ( size_t j = 0 ; j != this->point_cloud[i].size() ; ++j )
				{
					std::cerr << this->point_cloud[i][j] << " ";
				}	
				std::cerr << std::endl;		
			}
		}
	}
    double operator() ( const std::vector<double>& point )const
    {
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
    }
private:
	std::vector< std::vector<double> > point_cloud;
};
