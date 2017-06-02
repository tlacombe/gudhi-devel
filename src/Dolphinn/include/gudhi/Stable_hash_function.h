#ifndef DOLPHINN_HASH_H
#define DOLPHINN_HASH_H


#include <vector>
#include <chrono>
#include <random>
#include <cmath>
#include <iostream>
#include <numeric>
#include <unordered_map>
#include <string>
#include <thread>
#include <utility>
#include <Eigen/Dense>
#include <gudhi/euclidean_dist.h>

/**
 * Currently we have 2 possible families of hashing functions:
 * - https://en.wikipedia.org/wiki/Locality-sensitive_hashing#Stable_distributions
 * - hyperplane hashing
 */

namespace Gudhi {
namespace dolphinn {
/**  
 * \brief Hash functions used to build Dolphinn's hypercube.
 * 
 * \details
 * The class contains the methods and the informations to hash a point/dataset.
 *
 *\remark This class is meant to be abstracted.
 */

template <class T, typename Point>
class Stable_hash_function
{
		typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> Matrix;
    // of original pointset
    const int dimension; 
    double r;
		//for line LSH
    double b;
    // direction of the line for the stable distribution or offset for the hyperplanes
    Point a;
    //for Hyperplane LSH
    Matrix m;
    std::uniform_real_distribution<double> uni_distribution;
    std::uniform_int_distribution<int> uni_bit_distribution;
    std::default_random_engine generator;
    // key and a vector of the indices of the associated points
    std::unordered_map<int, std::vector<int> > hashtable;
    // for every key remember its random bit
    std::unordered_map<int, char> hashtable_for_random_bit;
    // Hamming cube vertex and vertices of assigned points.
    // This is used *only* by one of the hash.
    std::unordered_map< std::string, std::vector<int> > hashtable_cube;
  public:
  
  	//accessors
  	Point get_a() {
  		return a;
  	}
  	double get_r() {
  		return r;
  	}
  	double get_b() {
  		return b;
  	}
  	std::unordered_map< std::string, std::vector<int> > get_hashtable() {
  		return hashtable_cube;
  	}
  	Matrix get_m() {
  		return m;
  	}
  	
  	/** \brief Constructor that creates a 
  	 * vector from a stable distribution.
 	 *
	 * Assign 'D' random values from a normal distribution
	 * N(0,1/sqrt(D)).
	 *
	 * @param D  		 dimension of points
	 * @param r  		 parameter of Stable Distribution
	 * @param mean  	 optional parameter of Normal Distribution. Default is 0.0.
	 * @param deviation	optional parameter of Normal Distribution. Default is 1.0.
	 */
  	Stable_hash_function(const int D, const double r, const double mean = 0.0, const double deviation = 1.0)
  		: dimension(D), r(r), uni_distribution(0, r), uni_bit_distribution(0, 1),
  		generator(std::chrono::system_clock::now().time_since_epoch().count())
  	{
  		std::normal_distribution<typename std::conditional<std::is_same<T, double>::value, double, T>::type> distribution(mean, deviation);
      for(int i = 0; i < D; ++i)
      {
  			a.push_back(distribution(generator));
      }
      b = uni_distribution(generator);
  	}
  	
  	
  	/** \brief Constructor that creates a 
  	 * hashing function based on hyperplanes
 	 *
	 * @param K  		 dimension of the hypercube
	 * @param D  		 dimension of the points
	 * @param r  		 parameter hashing function (here r=0) 
	 */
  	Stable_hash_function(const int K, const int D, const double r)
  		: dimension(D), r(r), m(Matrix(K,D)),
  		generator(std::chrono::system_clock::now().time_since_epoch().count())
  	{
  		std::normal_distribution<double> distribution(0.0,1.0/std::sqrt((double)D));
  		for(int i =0; i<K;++i) {
  			for(int j =0; j<D;++j) {
  				m(i,j)=distribution(generator);
  			}
  		}
  	}

    /** \brief Constructor that creates a 
     * vector from a stable distribution,
     * in a distributed environment.
     *
     * Assign 'D' random values from a normal distribution
     * N(0,1/sqrt(D)).
     *
     * @param D            dimension of points
     * @param r            parameter of Stable Distribution
     * @param thread_info  Something that identifies the thread, so that every thread
     *                     creates its own random numbers, and not the same - as is the
     *                     case with just seeding the random generator with the time.
     * @param mean         optional parameter of Normal Distribution. Default is 0.0.
     * @param deviation    optional parameter of Normal Distribution. Default is 1.0.
    */
    Stable_hash_function(const int D, const double r, const int thread_info, const double mean = 0.0, const double deviation = 1.0)
      : dimension(D), r(r), uni_distribution(0, r), uni_bit_distribution(0, 1),
      generator(thread_info + std::chrono::system_clock::now().time_since_epoch().count())
    {
      std::normal_distribution<typename std::conditional<std::is_same<T, double>::value, double, T>::type> distribution(mean, deviation);
      for(int i = 0; i < D; ++i)
      {
        a.push_back(distribution(generator));
      }
      if(r!=0){b = uni_distribution(generator);}
    }
    
    
    /** \brief Hash a pointset with the hyperplanes method and fills the hypercube.
	 *
	 * @param v  	vector of points
	 */
    void hyperplane_hash(const std::vector<Point>& v) {
    	for(int j=0;j<dimension;++j){
    		a.push_back(v[0][j]);
    	}
    	for(size_t i=1;i<v.size();++i){
    		for(int j=0;j<dimension;++j){
		  		a[j]+=v[i][j];
		  	}
    	}
    	for(int j=0;j<dimension;++j){
    		a[j] = a[j]/v.size();
    	}
    	int j=0;
    	const double* aptr = &a[0];
    	for(const auto& x:v) {
    		const double* ptr = &x[0];
    		Eigen::Matrix<double, Eigen::Dynamic, 1> vector =  m * 
    										(Eigen::Map<const Eigen::Matrix<double ,Eigen::Dynamic, 1>>(ptr,dimension,1)-Eigen::Map<const Eigen::Matrix<double ,Eigen::Dynamic, 1>>(aptr,dimension,1));
    		std::vector<char> hash_key(vector.size());
    		for(int i =0; i<vector.size();++i) {
    			if(vector(i,0)>0){
    				hash_key[i]=1;
    			} else {
    				hash_key[i]=0;
    			}
    		}
    		hashtable_cube[std::string(hash_key.begin(),hash_key.end())].push_back(j);
    		++j;
    	}
    }
    
    /** \brief Hash a point with the hyperplanes method.
	 *
	 * @param x 				point to hash
	 * @param hash_key  stores the result of the hashing
	 */
    template<typename bitT>
    void hyperplane_hash(const Point& x, std::vector<bitT>& hash_key) {
    	const double* aptr = &a[0];
   		const double* ptr = &x[0];
  		Eigen::Matrix<double, Eigen::Dynamic, 1> vector =  m * 
  										(Eigen::Map<const Eigen::Matrix<double ,Eigen::Dynamic, 1>>(ptr,dimension,1)-Eigen::Map<const Eigen::Matrix<double ,Eigen::Dynamic, 1>>(aptr,dimension,1));
  		for(int i =0; i<vector.size();++i) {
  			if(vector(i,0)>0){
  				hash_key[i]=1;
  			} else {
  				hash_key[i]=0;
  			}
  		}
  	}

  	/** \brief Hash a pointset with the stable distribution method.
	 *
	 * @param v  	vector of points
	 * @param N   number of points
	 * @param D   dimension of points
	 */
  	void hash(const std::vector<Point>& v, const int N, const int D)
  	{
  		for(int i = 0; i < N; ++i)
  		{
  			hashtable[hash(v[i])].push_back(i);
  		}
  	}

  	/** \brief Hash a point.
	 *
	 * @param v  			 the point to be hashed
	 * @return 		     result of hash function
	 */
  	
  	int hash(const Point v)
  	{
  		double scalar_product = 0;
  		for(size_t i=0;i<v.size();++i){
				scalar_product += v[i] * a[i];
			}
  		//std::cout << scalar_product << " " << b << " " << r << std::endl;
  		//std::cout << (scalar_product + b) << " " << (scalar_product + b) / r << " " << (scalar_product + b) / (double)r << std::endl;
  		if(r>0){
  			return floor((scalar_product + b) / r);
  		} else {
  			if(scalar_product>0){return 1;} else {return 0;}
  		}
  	}

   	/** \brief Assing random bit for every key.
 	 *
 	 * @param v 	vector of (to be) mapped points
 	 * @param k		iteration (assign the k-th bit of the points)
 	 * @param K   dimension of the cube that awaits for
 	 *				  the points to be mapped on its vertices
	 */
    template<typename bitT>
  	void assign_random_bit(std::vector<bitT>& v, const int k, const int K)
  	{
  		bitT random_bit;
  		if(r>0) {
				for(auto& key_value: hashtable)
				{
					random_bit = uni_bit_distribution(generator);
		      hashtable_for_random_bit[key_value.first] = random_bit;
					for(auto const& point_idx: key_value.second)
					{
						v[k + point_idx * K] = random_bit;
					}
				}
  		} else {
  			for(auto& key_value: hashtable)
				{
					if(key_value.first>0) {
						random_bit = 1;
					}	else {
						random_bit = 0;
					}
		      hashtable_for_random_bit[key_value.first] = random_bit;
					for(auto const& point_idx: key_value.second)
					{
						v[k + point_idx * K] = random_bit;
					}
				}
  		}
  	} 

    /** \brief Assing random bit for every key and fill cube's hashtable.
     *
     * @param v   vector of (to be) mapped points
     * @param K   dimension of the cube that awaits for
     *          the points to be mapped on its vertices
    */
    template<typename bitT>
    void assign_random_bit_and_fill_hashtable_cube(std::vector<bitT>& v, const int K)
    {
      bitT random_bit;
      if(r>0){
		    for(auto& key_value: hashtable)
		    {
		      random_bit = uni_bit_distribution(generator);
		      hashtable_for_random_bit[key_value.first] = random_bit;
		      for(auto const& point_idx: key_value.second)
		      {
		        v[(K - 1) + point_idx * K] = random_bit;
		        //std::cout << (int)(std::string(v.begin() + point_idx * K, v.begin() + (point_idx + 1) * K))[2] << std::endl;
		        hashtable_cube[std::string(v.begin() + point_idx * K, v.begin() + (point_idx + 1) * K)].push_back(point_idx);
		      }
		    }
      } else {
      	for(auto& key_value: hashtable)
		    {
		      if(key_value.first>0) {
		      	random_bit = 1;
		      } else {
		      	random_bit = 0;
		      }
		      hashtable_for_random_bit[key_value.first] = random_bit;
		      for(auto const& point_idx: key_value.second)
		      {
		        v[(K - 1) + point_idx * K] = random_bit;
		        //std::cout << (int)(std::string(v.begin() + point_idx * K, v.begin() + (point_idx + 1) * K))[2] << std::endl;
		        hashtable_cube[std::string(v.begin() + point_idx * K, v.begin() + (point_idx + 1) * K)].push_back(point_idx);
		      }
		    }
      }
    }	

    /** \brief Assing random bit for queries.
   *
   * @param q   							query
   * @param mapped_q_begin   	(to be) mapped query
   * @param k   							iteration (assign the k-th bit of the query)
   */
    template <typename bit_iterator>
    void assign_random_bit_query(const Point q, bit_iterator mapped_q_begin, const int k)
    {
      int q_key = hash(q);
      const auto& q_key_it = hashtable_for_random_bit.find(q_key);
    	if(q_key_it != hashtable_for_random_bit.end())
      {
      	*(mapped_q_begin + k) = q_key_it->second;
    	}
    	else
    	{
    		*(mapped_q_begin + k) = uni_bit_distribution(generator);
    	}
    }   
    
   // hash(const std::vector<Point>& v, const int N, const int D)
    //assign_random_bit(std::vector<bitT>& v, const int k, const int K)

    /** \brief Radius query the Hamming cube.
      *
      * @param mapped_query        mapped query
      * @param radius              find a point within r with query
      * @param K                   dimension of the mapped query
      * @param MAX_PNTS_TO_SEARCH  threshold
      * @param pointset            original points
      * @param query_point         original query
      * @return                    index of a point, where Eucl(point[i], query_point) <= r
    */
    template <typename iterator>
    int radius_query(std::string mapped_query, const double radius, const int K, const int MAX_PNTS_TO_SEARCH, iterator pointset, iterator query_point)
    {
      int points_checked = 0;
      int answer_point_idx = -1;
      double squared_radius = radius * radius;
      const auto& q_key_it = hashtable_cube.find(mapped_query);
      // search query's cube vertex, if pointsets' points exist there
      if(q_key_it != hashtable_cube.end())
      {
        //print_string_cast_int(q_key_it->first); std::cout << " " << q_key_it->second.size() << std::endl;
        answer_point_idx = Euclidean_distance_within_radius<iterator>(pointset, q_key_it->second, dimension, query_point, squared_radius, MAX_PNTS_TO_SEARCH);
        //if(answer_point_idx!=-1) std::cout << squared_Eucl_distance(query_point, query_point + dimension, pointset + answer_point_idx * dimension) << std::endl;
        points_checked += q_key_it->second.size();
      }
      // check neighboring vertices from query's cube vertex
      int Hamming_dist = 1;
      while (points_checked < MAX_PNTS_TO_SEARCH && answer_point_idx == -1)
      {
        find_strings_with_fixed_Hamming_dist_for_radius_query<iterator>(mapped_query, K - 1, Hamming_dist++, points_checked, MAX_PNTS_TO_SEARCH, squared_radius, pointset, query_point, answer_point_idx);
      }
      //std::cout << "ANSWER = " << answer_point_idx << ", checked points = " << points_checked << std::endl;
      return answer_point_idx;
    }

    /** \brief Find strings within a given Hamming distance. Used by 'radius_query()'.
      *
      * @param str                 given string
      * @param i                   index
      * @param changesLeft         changes left to make
      * @param points_checked      current points checked
      * @param MAX_PNTS_TO_SEARCH  threshold
      * @param squared_radius      check if any original point lies in r Euclidean distance from the original query
      * @param pointset						 original dataset
      * @param query_point				 pointer to the queried point
      * @param answer_point_idx    index of point that has distance less or equal than r with the query
    */
    template <typename iterator>
    bool find_strings_with_fixed_Hamming_dist_for_radius_query(std::string& str, const int i, const int changesLeft, 
      int& points_checked, const int MAX_PNTS_TO_SEARCH, const double squared_radius, iterator& pointset, 
      iterator& query_point, int& answer_point_idx)
    {
      bool stop = false;
      if (changesLeft == 0) {
        //print_string_cast_int(str); std::cout << std::endl;
        const auto& key_value_it = hashtable_cube.find(str);
        //std::cout << "start checking\n";
        if(key_value_it != hashtable_cube.end())
        {
          //std::cout << " " << key_value_it->second.size() << std::endl;
          answer_point_idx = Euclidean_distance_within_radius<iterator>(pointset, key_value_it->second, dimension, query_point, squared_radius, MAX_PNTS_TO_SEARCH);
          //if(answer_point_idx!=-1) std::cout << squared_Eucl_distance(query_point, query_point + dimension, pointset + answer_point_idx * dimension) << std::endl;
          points_checked += key_value_it->second.size();
          //std::cout << "check: " << points_checked << " " << MAX_PNTS_TO_SEARCH << std::endl;
          stop = (answer_point_idx != -1 || points_checked > MAX_PNTS_TO_SEARCH);
          //std::cout << "stop = " << stop << std::endl;
        }
        return stop;
      }
      if (i < 0)
        return 0;
      // flip current bit
      if(!stop)
      {
        str[i] ^= 1;
        stop = find_strings_with_fixed_Hamming_dist_for_radius_query(str, i-1, changesLeft-1, points_checked, MAX_PNTS_TO_SEARCH, squared_radius, pointset, query_point, answer_point_idx);
      }
      // or don't flip it (flip it again to undo)
      if(!stop)
      {
        str[i] ^= 1;
        stop = find_strings_with_fixed_Hamming_dist_for_radius_query(str, i-1, changesLeft, points_checked, MAX_PNTS_TO_SEARCH, squared_radius, pointset, query_point, answer_point_idx);
      }
      return stop;
    }

    /** \brief Nearest Neighbor query the Hamming cube.
      *
      * @param mapped_query        mapped query
      * @param K                   dimension of the mapped query
      * @param m                   number of nearest neighbors to search
      * @param MAX_PNTS_TO_SEARCH  threshold
      * @param pointset            original points
      * @param query_point         original query
      * @return                    index and distance from query of (approximate) Nearest Neighbor
    */
    template <typename iterator>
    std::vector<std::pair<int, double>>m_nearest_neighbors_query(std::string mapped_query, const int K, const int m, const int MAX_PNTS_TO_SEARCH, iterator pointset, iterator query_point)
    {
      int points_checked = 0;
      std::vector<std::pair<int, double>> answer_point_idx_dist(m, std::make_pair(-1, 1000000.0));
      const auto& q_key_it = hashtable_cube.find(mapped_query);
      // search query's cube vertex, if pointsets' points exist there
      if(q_key_it != hashtable_cube.end())
      {
    	  find_M_Nearest_Neighbor_indices<iterator>(pointset, q_key_it->second, dimension, m, query_point, answer_point_idx_dist, MAX_PNTS_TO_SEARCH);
    	  points_checked += q_key_it->second.size();
      }
      // check neighboring vertices from query's cube vertex
      int Hamming_dist = 1;
      while (points_checked < MAX_PNTS_TO_SEARCH)
      {
        find_strings_with_fixed_Hamming_dist_for_nearest_neighbor_query<iterator>(mapped_query, K - 1, m, Hamming_dist++, points_checked, MAX_PNTS_TO_SEARCH, pointset, query_point, answer_point_idx_dist);
      }
      return answer_point_idx_dist;
    }

    /** \brief Find strings within a given Hamming distance. Used by 'nearest_neighbor_query()'.
      *
      * @param str                     given string
      * @param i                       index
      * @param m                       number of nearest neighbors to search
      * @param changesLeft             changes left to make
      * @param points_checked          current points checked
      * @param MAX_PNTS_TO_SEARCH      threshold
      * @param pointset           		 original points
      * @param query_point        		 original query
      * @param answer_point_idx_dist   index and distance of current best Nearest Neighbor
    */
    template <typename iterator>
    bool find_strings_with_fixed_Hamming_dist_for_nearest_neighbor_query(std::string& str, const int i, const int m,  const int changesLeft,
      int& points_checked, const int MAX_PNTS_TO_SEARCH, iterator& pointset, 
      iterator& query_point,std::vector<std::pair<int, double>>& answer_point_idx_dist)
    {
      bool stop = false;
      if (changesLeft == 0) {
        const auto& key_value_it = hashtable_cube.find(str);
        if(key_value_it != hashtable_cube.end())
        {
          find_M_Nearest_Neighbor_indices<iterator>(pointset, key_value_it->second, dimension, m, query_point, answer_point_idx_dist, MAX_PNTS_TO_SEARCH);
          points_checked += key_value_it->second.size();
          stop = (points_checked > MAX_PNTS_TO_SEARCH);
        }
        return stop;
      }
      if (i < 0)
        return 0;
      // flip current bit
      if(!stop)
      {
        str[i] ^= 1;
        stop=find_strings_with_fixed_Hamming_dist_for_nearest_neighbor_query(str, i-1, m, changesLeft-1, points_checked, MAX_PNTS_TO_SEARCH, pointset, query_point, answer_point_idx_dist);
      }
      // or don't flip it (flip it again to undo)
      if(!stop)
      {
        str[i] ^= 1;
        stop=find_strings_with_fixed_Hamming_dist_for_nearest_neighbor_query(str, i-1, m, changesLeft, points_checked, MAX_PNTS_TO_SEARCH, pointset, query_point, answer_point_idx_dist);
      }
      return stop;
    }

    /** \brief Check if vector is full of 'value'.
     *
     * @param vec   vector to be checked
     * @param value check if all elements of 'vec' have this value
     * @return      True if all elements of 'vev' are 'value'. False, otherwise.
    */
    bool check_vec(std::vector<int>& vec, int value)
    {
      for (auto& v: vec)
        if(v != value)
          return false;
      return true;
    }

    /** \brief Return first element of 'vec' that is different from 'value'.
     *
     * @param vec   vector to be checked
     * @param value check if all elements of 'vec' have this value
     * @return      Element of 'vec' different from 'value', if exists. 'value', otherwise.
    */
    int find_non_value_in_vec(std::vector<int>& vec, int value)
    {
      for (auto& v: vec)
        if(v != value)
          return v;
      return value;
    }

  	/** \brief Print 'hashtable'.
 	 *
	 */
  	void print()
  	{
  		std::cout << "Hashtable for h, is:\n";
  		for (auto& x: hashtable)
  		{
    		std::cout << x.first << ": ";
    		for(auto& y: x.second)
    		{
      			std::cout << y << ", ";
    		}
    		std::cout << "\n";
  		}
  	}

  	/** \brief Print number of keys and number of values per key.
 	 *
	 */
  	void print_stats()
  	{
  		std::cout << "Hashtable is of size (|keys|) = " << hashtable.size() << "\n";
  		for(auto& key_value: hashtable)
  			std::cout << key_value.first << " has " << key_value.second.size() << " values/points\n";
  		std::cout << "\n";
  	}

    /** \brief Print hashtable of Hamming cube. 
    * @param print_indices Print all the values of the unordered_map. Default false.
    *
   */
    void print_hashtable_cube(const bool print_indices = false)
    {
      if(!hashtable_cube.size())
      {
        std::cout << "Are you sure is this the last hash? The Hamming cube's hashtable is created ";
        std::cout << "only after the last, (K-1)-th, bit of the mapped point is set. Hashtable won't print\n";
        return;
      }
      std::cout << "Hashtable is of size (|keys|) = " << hashtable_cube.size() << "\n";
      for(auto& key_value: hashtable_cube)
      {
        print_string_cast_int(key_value.first); std::cout << " has " << key_value.second.size() << " values/points\n";
        if(print_indices)
          for(auto& point_idx: key_value.second)
            std::cout << point_idx << " ";
      }
      std::cout << "\n";
    }

  	/** \brief Print vector 'a'.
 	 *
	 */
  	void print_a()
  	{
  		std::cout << "'a' vector of the stable distribution for h, is:\n";
  		for(unsigned int i = 0; i < dimension; ++i)
  			std::cout << a[i] << " ";
  		std::cout << "\n";
  	}
  	
  	void print_string_cast_int(const std::string& str)
		{
			for(auto& character: str)
    		std::cout << (int)character;
		}

};
}
}
#endif /*DOLPHINN_HASH_H*/
