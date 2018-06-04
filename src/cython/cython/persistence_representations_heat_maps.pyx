from cython cimport numeric
from libcpp.vector cimport vector
from libcpp.utility cimport pair
from libcpp cimport bool  
from cython.operator cimport dereference as deref
import os
import sys

"""This file is part of the Gudhi Library. The Gudhi library
   (Geometric Understanding in Higher Dimensions) is a generic C++
   library for computational topology.

   Author(s):       Pawel Dlotko

   Copyright (C) 2018 Swansea University

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

__author__ = "Pawel Dlotko"
__copyright__ = "Copyright (C) 2018 Swansea University"
__license__ = "GPL v3"



cdef extern from "Persistence_heat_maps_interface.h" namespace "Gudhi::Persistence_representations":
	cdef cppclass Persistence_heat_maps_interface "Gudhi::Persistence_representations::Persistence_heat_maps_interface":
		Persistence_heat_maps_interface()
		#For some reasons such a constructors are not accepted. That is why we use the static methods from below instead. 
		#Persistence_heat_maps_interface(const vector[pair[double, double]]& interval,size_t how_many_pixels_raidus_of_Gausian_kernel,size_t number_of_pixels,double min_ = 0,double max_ = 0)																										  
		#Persistence_heat_maps_interface(const char* filename, size_t how_many_pixels_raidus_of_Gausian_kernel,size_t number_of_pixels,double min_ = 0,double max_ = 0,unsigned dimensions = 0)

		void compute_mean_interface(const vector[Persistence_heat_maps_interface*]& maps_)
		void compute_median_interface(const vector[Persistence_heat_maps_interface*]& maps_)
		void compute_percentage_of_active_interface(const vector[Persistence_heat_maps_interface*]& maps_, size_t cutoff )
		void print_to_file(const char* filename) const
		void load_from_file(const char* filename)
		bool check_if_the_same(const Persistence_heat_maps_interface& second) const 
		double get_min() const
		double get_max() const
		vector[double] vectorize(int number_of_function) const
		size_t number_of_vectorize_functions() const
		double project_to_R(int number_of_function) const
		size_t number_of_projections_to_R() const
		double distance(const Persistence_heat_maps_interface& second_, double power) const
		void compute_average_interface(const vector[Persistence_heat_maps_interface*]& to_average)
		double compute_scalar_product(const Persistence_heat_maps_interface& second_) const
		pair[double, double] get_x_range() const
		pair[double, double] get_y_range() const
		bool compare( const Persistence_heat_maps_interface& second )const 
		
		#**************
        #static methods
		@staticmethod
		Persistence_heat_maps_interface* construct_from_file( const char* filename, size_t how_many_pixels_raidus_of_Gausian_kernel,
							size_t number_of_pixels,
							double min_,
							double max_,
							int dimensions)
		@staticmethod
		Persistence_heat_maps_interface* construct_from_vector_of_pairs( const vector[pair[double, double]]& interval,
							size_t how_many_pixels_raidus_of_Gausian_kernel,
							size_t number_of_pixels,
							double min_,
							double max_)
        #***************

        
        
cdef class PersistenceHeatMaps:

	cdef Persistence_heat_maps_interface* thisptr

#Can we have only one constructor, or can we have more
	def __init__(self, how_many_pixels_raidus_of_Gausian_kernel=5, number_of_pixels=1000, 
	vector_of_intervals=None, file_with_intervals='', min_=0, max_=0,dimension=-1):
		"""
		This is a class implementing persistence heat maps. At every point of a diagram we place a Gaussiian kernel (standard deviation is given in terms of the number of pixels).
		This way we obtain yet another representation of persistence that can be used in data analysis. 
		"""	

	def __cinit__(self, how_many_pixels_raidus_of_Gausian_kernel=5, number_of_pixels=1000, 
	vector_of_intervals=None, file_with_intervals='', min_=0, max_=0,dimension=-1):
		"""
		This is a constructor of a class PersistenceHeatMaps.
		The persistence intervals can be input to the class by either providing a filename
		of a file storing them, or by giving them explicitelly as a vector of pairs.
		Additional parameters required to construct Persistence Heat Maps are:
		:param how_many_pixels_raidus_of_Gausian_kernel - indicate the diameter of the Gausian kernel
		to be placed on the points of persistence diagram.
		:param number_of_pixels - size of the persistence image.
		:param max_  an optional parameter determining the upper right corner of the image (max_,max_). 
		If not set, will be computed based on the input data. 
		:param min_ an optional parameter determining the lower left corner of the image (min_,min_). 
		If not set, will be computed based on the input data. 
		:param dimension - na optional parameter, a dimension of the intervals to be read from a file.           
		"""
		if ( (vector_of_intervals is None) and ( file_with_intervals is '' ) ):	
			print "Creating empty persistence heat map." 
			self.thisptr = new Persistence_heat_maps_interface()
		else:
			if (vector_of_intervals is None) and (file_with_intervals is not ''):
				if os.path.isfile(file_with_intervals):
					if (min_ is not None) and (max_ is not None):							
						self.thisptr = Persistence_heat_maps_interface.construct_from_file(str.encode(file_with_intervals), how_many_pixels_raidus_of_Gausian_kernel, number_of_pixels, min_, max_, dimension)
					else:
						self.thisptr = Persistence_heat_maps_interface.construct_from_file(str.encode(file_with_intervals), how_many_pixels_raidus_of_Gausian_kernel, number_of_pixels,0,0,0)
				else:
					print("file " + file_with_intervals + " not found."	)			
			else:            
				if (file_with_intervals is '') and (vector_of_intervals is not None):
					if (min_ is not None) and (max_ is not None):
						self.thisptr = Persistence_heat_maps_interface.construct_from_vector_of_pairs(vector_of_intervals, how_many_pixels_raidus_of_Gausian_kernel, number_of_pixels, min_, max_)
					else:
						self.thisptr = Persistence_heat_maps_interface.construct_from_vector_of_pairs(vector_of_intervals, how_many_pixels_raidus_of_Gausian_kernel, number_of_pixels,0,0)
				else:
					self.thisptr = new Persistence_heat_maps_interface()


	def __dealloc__(self):
		"""
		destructor
		"""
		if self.thisptr != NULL:
			del self.thisptr   
            
	def compute_mean( self,maps_=[] ):
		"""
		Compute a mean value of a collection of heat maps and store it in 
		the current object. Note that all the persistence maps send in a 
		vector to this procedure need to have the same parameters.
		If this is not the case, the program will throw an exception.
		:param vector of PersistenceHeatMaps
		:type PersistenceHeatMaps
		"""            
		cdef vector[Persistence_heat_maps_interface*] cpp_list    
		if ( self.thisptr != NULL ) and ( maps_ is not None ):	
			for elt in maps_: 
				cpp_list.push_back((<PersistenceHeatMaps>elt).thisptr)
			self.thisptr.compute_mean_interface( cpp_list )  
            
	def compute_median( self,maps_=[] ):
		"""
		Compute a median of a collection of heat maps and store it in 
		the current object. Note that all the persistence maps send in a 
		vector to this procedure need to have the same parameters.
		If this is not the case, the program will throw an exception.
		:param vector of PersistenceHeatMaps
		:type PersistenceHeatMaps
		"""            
		cdef vector[Persistence_heat_maps_interface*] cpp_list    
		if ( self.thisptr != NULL ) and ( maps_ is not None ):	
			for elt in maps_: 
				cpp_list.push_back((<PersistenceHeatMaps>elt).thisptr)
			self.thisptr.compute_median_interface( cpp_list ) 
            
	def compute_percentage_of_active( self,maps_=[] , cutoff = 1 ):
		"""
		Compute a median of a collection of heat maps and store it in 
		the current object. Note that all the persistence maps send in a 
		vector to this procedure need to have the same parameters.
		If this is not the case, the program will throw an exception.
		:param vector of PersistenceHeatMaps
		:type PersistenceHeatMaps
		"""            
		cdef vector[Persistence_heat_maps_interface*] cpp_list    
		if ( self.thisptr != NULL ) and ( maps_ is not None ):	
			for elt in maps_: 
				cpp_list.push_back((<PersistenceHeatMaps>elt).thisptr)
			self.thisptr.compute_percentage_of_active_interface( cpp_list , cutoff)  
                       
	def load_from_file(self,filename):
		"""
		This procedure loads a persistence heat map from file. It erase all the data
		that was previously stored in this vector.
		:param Name of the file.
		:type String
		"""
		if ( self.thisptr != NULL ) and ( filename is not None ):
			self.thisptr.load_from_file(filename)   
		        

	def print_to_file(self,filename) :
		"""
		The procedure stores a persistence heat map to a file. The file can be later
		used by a procedure load_from_file.
		:param Name of the file.
		:type String
		"""
		if ( self.thisptr != NULL ) and ( filename is not None ):
			self.thisptr.print_to_file(filename)                          
       
	def check_if_the_same( self, PersistenceHeatMaps second ):
		"""
		Check if this Persistence Heat Map and the second Persistence Heat Map
		are the same.
		:type bool
		"""
		if ( (self.thisptr != NULL) and (second is not None) and (second.thisptr != NULL) ):
			return self.thisptr.check_if_the_same( deref(second.thisptr) )
            
	def get_min( self ):
		"""
		Returns the minimal value of Persistence Heat Map.
		:type double
		"""
		if ( self.thisptr != NULL ):
			return self.thisptr.get_min()      
            
	def get_max( self ):
		"""
		Returns the minimal value of Persistence Heat Map.
		:type double
		"""
		if ( self.thisptr != NULL ):
			return self.thisptr.get_max()              
            
	def vectorize(self, number_of_function):
		"""
		This function returns a vector. It
		is required in a concept
		Vectorized_topological_data
		:param number of function
		:type vector of doubles
		"""
		if ( self.thisptr != NULL ) and ( number_of_function is not None ):
			return self.thisptr.vectorize(number_of_function)

	def number_of_vectorize_functions(self):
		"""
		The number of projections to R is defined to be the length
		of the vector. 
		"""
		if ( self.thisptr != NULL ):
			return self.thisptr.number_of_vectorize_functions()   
                        
	def project_to_R(self, number_of_function):
		"""
		The number of projections to R is defined to the number of nonzero
		entries in vector. I-th projection is the I-th entry
		This function is required by the Real_valued_topological_data concept.
		At the moment this function is not tested, since it is quite likely
		to be changed in the future. Given this, when
		using it, keep in mind that it
		will be most likely changed in the next versions.
		:param number of function
		:type doubles
		"""
		if ( self.thisptr != NULL ) and ( number_of_function is not None ):
			return self.thisptr.project_to_R(number_of_function)

	def number_of_projections_to_R(self):
		"""
		The function gives the number of possible projections to R. This
		function is required by the
		Real_valued_topological_data concept
		:type integer
		"""
		if ( self.thisptr != NULL ):
			return self.thisptr.number_of_projections_to_R() 
     
    #I havev problem with getting this function compiled due to a strange error                     
	def distance(self, PersistenceHeatMaps second, power = 1):
		"""
		A function to compute distance between Persistence Heat Maps.
		The parameter of this function is a PersistenceHeatMaps.
		This function is required in Topological_data_with_distances concept.
		For max norm distance, set power to numeric_limits<double>::max()
		:param the heat map to compute distance to
		:type double
		"""
		if ( self.thisptr != NULL ) and ( second is not None ) and ( power is not None ):
			return self.thisptr.distance(deref(second.thisptr), power)

	def compute_average( self,to_average=[] ):
		"""
		A function to compute averaged Persistence Heat Maps, based on vector
		of Persistence Heat Maps.
		This function is required by Topological_data_with_averages concept.
		:param vector of persistence vectors to average
		:type PersistenceHeatMaps
		"""            
		cdef vector[Persistence_heat_maps_interface*] cpp_list    
		if ( self.thisptr != NULL ) and ( to_average is not None ):	
			for elt in to_average: 
				cpp_list.push_back((<PersistenceHeatMaps>elt).thisptr)
			self.thisptr.compute_average_interface( cpp_list )   
                                 
	def compute_scalar_product(self, PersistenceHeatMaps second):
		"""
		A function to compute scalar product of Persistence Heat Maps.
		The parameter of this function is a PersistenceHeatMap.
		This function is required in Topological_data_with_scalar_product concept.
		:param the heat map to compute scalar product with
		:type double
		"""
		if ( self.thisptr != NULL ) and ( second is not None ):
			return self.thisptr.compute_scalar_product( deref(second.thisptr) )
                        
	def get_x_range( self ):
		"""
		Returns the x range of Persistence Heat Maps.
		:type pair[double,double]
		"""
		if ( self.thisptr != NULL ):
			return self.thisptr.get_x_range() 
                        
	def get_y_range( self ):
		"""
		Returns the y range of Persistence Heat Maps.
		:type pair[double,double]
		"""
		if ( self.thisptr != NULL ):
			return self.thisptr.get_y_range()                               
			
	def compare( self , PersistenceHeatMaps second ):
		"""
		Returns true if the second Persistence Heat Map is identical to
		the first one, and false in the other case. 
		"""
		if ( self.thisptr != NULL ): 
			return self.thisptr.compare( deref(second.thisptr) ) 
