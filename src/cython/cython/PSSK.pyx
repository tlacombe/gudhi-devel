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



cdef extern from "PSSK_interface.h" namespace "Gudhi::Persistence_representations":
	cdef cppclass PSSK_interface "Gudhi::Persistence_representations::PSSK_interface":
		PSSK()
		
		void compute_mean_interface(const vector[PSSK_interface*]& maps_)
		void compute_median_interface(const vector[PSSK_interface*]& maps_)
		void compute_percentage_of_active_interface(const vector[PSSK_interface*]& maps_, size_t cutoff )
		
		void print_to_file(const char* filename) const
		void load_from_file(const char* filename)
		bool check_if_the_same(const PSSK_interface& second) const 
		double get_min() const
		double get_max() const
		vector[double] vectorize(int number_of_function) const
		size_t number_of_vectorize_functions() const
		double project_to_R(int number_of_function) const
		size_t number_of_projections_to_R() const
		double distance(const PSSK_interface& second_, double power) const
		void compute_average_interface(const vector[PSSK_interface*]& to_average)
		double compute_scalar_product(const PSSK_interface& second_) const
		pair[double, double] get_x_range() const
		pair[double, double] get_y_range() const
		
		#**************
        #static methods
		@staticmethod
		PSSK_interface* construct_from_file( const char* filename, size_t how_many_pixels_raidus_of_Gausian_kernel,
							size_t number_of_pixels,
							double min_,
							double max_,
							int dimensions)
		@staticmethod
		PSSK_interface* construct_from_vector_of_pairs( const vector[pair[double, double]]& interval,
							size_t how_many_pixels_raidus_of_Gausian_kernel,
							size_t number_of_pixels,
							double min_,
							double max_)
        #***************

        
        
cdef class PSSK:

	cdef PSSK_interface* thisptr


#Can we have only one constructor, or can we have more
	def __init__(self, how_many_pixels_raidus_of_Gausian_kernel, number_of_pixels, 
	min_=None, max_=None,dimension=-1):
		"""
		This is a version of a representation presented in https://arxiv.org/abs/1412.6821
		In that paper the authors are using the representation just to compute kernel. Over here, we extend the usability by
		far.
		Note that the version presented here is not exact, since we are discretizing the kernel.
		The only difference with respect to the original class is the method of creation. We have full (square) image, and for
		every point (p,q), we add a kernel at (p,q) and the negative kernel
		at (q,p).
		"""	

	def __cinit__(self, how_many_pixels_raidus_of_Gausian_kernel, number_of_pixels, 
	vector_of_intervals=None, file_with_intervals='', min_=None, max_=None,dimension=-1):
		"""
		This is a constructor of a class PSSK.
		The persistence intervals can be input to the class by either providing a filename
		of a file storing them, or by giving them explicitelly as a vector of pairs.
		Additional parameters required to construct PSSK are:
		:param how_many_pixels_raidus_of_Gausian_kernel - indicate the diameter of the Gausian kernel
		to be placed on the points of persistence diagram.
		:param number_of_pixels - size of the persistence image.
		:param max_  an optional parameter determining the upper right corner of the image (max_,max_). 
		If not set, will be computed based on the input data. 
		:param min_ an optional parameter determining the lower left corner of the image (min_,min_). 
		If not set, will be computed based on the input data. 
		:param dimension - na optional parameter, a dimension of the intervals to be read from a file.           
		"""
		if ( (vector_of_intervals is None) or ( file_with_intervals is '' ) ):
			print "Please provide parameters to construct the persistence vectors." 
		else:
			if (vector_of_intervals is None) and (file_with_intervals is not ''):				
				if os.path.isfile(file_with_intervals):
					if (min_ is not None) and (max_ is not None):
						self.thisptr = PSSK_interface.construct_from_file(str.encode(file_with_intervals), how_many_pixels_raidus_of_Gausian_kernel, number_of_pixels, min_, max_, dimension)
					else:
						self.thisptr = PSSK_interface.construct_from_file(str.encode(file_with_intervals), how_many_pixels_raidus_of_Gausian_kernel, number_of_pixels,0,0,0)
				else:
					print("file " + file_with_intervals + " not found.")				
			else:            
				if (file_with_intervals is '') and (vector_of_intervals is not None):
					if (min_ is not None) and (max_ is not None):
						self.thisptr = PSSK_interface.construct_from_vector_of_pairs(vector_of_intervals, how_many_pixels_raidus_of_Gausian_kernel, number_of_pixels, min_, max_)
					else:
						self.thisptr = PSSK_interface.construct_from_vector_of_pairs(vector_of_intervals, how_many_pixels_raidus_of_Gausian_kernel, number_of_pixels,0,0)
				else:
					self.thisptr = new PSSK_interface()

	def __dealloc__(self):
		"""
		destructor
		"""
		if self.thisptr != NULL:
			del self.thisptr   
            
	def compute_mean( self,maps_=[] ):
		"""
		Compute a mean value of a collection of PSSK and store it in 
		the current object. Note that all the PSSK added to a 
		vector to this procedure need to have the same parameters.
		If this is not the case, the program will throw an exception.
		:param vector of PSSK
		:type PersistenceHeatMaps
		"""            
		cdef vector[PSSK_interface*] cpp_list    
		if ( self.thisptr != NULL ) and ( maps_ is not None ):	
			for elt in maps_: 
				cpp_list.push_back((<PSSK>elt).thisptr)
			self.thisptr.compute_mean_interface( cpp_list )  
            
	def compute_median( self,maps_=[] ):
		"""
		Compute a median of a collection of PSSK and store it in 
		the current object. Note that all the PSSK added to a 
		vector to this procedure need to have the same parameters.
		If this is not the case, the program will throw an exception.
		:param vector of PSSL
		:type PSSK
		"""            
		cdef vector[PSSK_interface*] cpp_list    
		if ( self.thisptr != NULL ) and ( maps_ is not None ):	
			for elt in maps_: 
				cpp_list.push_back((<PSSK>elt).thisptr)
			self.thisptr.compute_median_interface( cpp_list ) 
            
	def compute_percentage_of_active( self,maps_=[] , cutoff = 1 ):
		"""
		Compute a median of a collection of PSSK and store it in 
		the current object. Note that all the PSSK send in a 
		vector to this procedure need to have the same parameters.
		If this is not the case, the program will throw an exception.
		:param vector of PersistenceHeatMaps
		:type PSSK
		"""            
		cdef vector[PSSK_interface*] cpp_list    
		if ( self.thisptr != NULL ) and ( maps_ is not None ):	
			for elt in maps_: 
				cpp_list.push_back((<PSSK>elt).thisptr)
			self.thisptr.compute_percentage_of_active_interface( cpp_list , cutoff)  
                       
	def load_from_file(self,filename):
		"""
		This procedure loads a PSSK from file. It erase all the data
		that was previously stored in this vector.
		:param Name of the file.
		:type String
		"""
		if ( self.thisptr != NULL ) and ( filename is not None ):
			self.thisptr.load_from_file(filename)            

	def print_to_file(self,filename) :
		"""
		The procedure stores a PSSK to a file. The file can be later
		used by a procedure load_vector_from_file.
		:param Name of the file.
		:type String
		"""
		if ( self.thisptr != NULL ) and ( filename is not None ):
			self.thisptr.print_to_file(filename)                          
       
	def check_if_the_same( self, PSSK second ):
		"""
		Check if this PSSK and the second PSSK
		are the same.
		:type bool
		"""
		if ( (self.thisptr != NULL) and (second is not None) and (second.thisptr != NULL) ):
			return self.thisptr.check_if_the_same( deref(second.thisptr) )
            
	def get_min( self ):
		"""
		Returns the minimal value of PSSK.
		:type double
		"""
		if ( self.thisptr != NULL ):
			return self.thisptr.get_min()      
            
	def get_max( self ):
		"""
		Returns the minimal value of PSSK.
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
	#def distance(self, PSSK second, power):
	#	"""
	#	A function to compute distance between PSSK.
	#	The parameter of this function is a PSSK.
	#	This function is required in Topological_data_with_distances concept.
	#	For max norm distance, set power to numeric_limits<double>::max()
	#	:param the heat map to compute distance to
	#	:type double
	#	"""
	#	if ( self.thisptr != NULL ) and ( second is not None ) and ( power is not None ):
	#		return self.thisptr.distance(deref(second.thisptr), power) 	   
                        
	def compute_average( self,to_average=[] ):
		"""
		A function to compute averaged PSSK, based on vector
		of PSSK.
		This function is required by Topological_data_with_averages concept.
		:param vector of persistence vectors to average
		:type PersistenceHeatMaps
		"""            
		cdef vector[PSSK_interface*] cpp_list    
		if ( self.thisptr != NULL ) and ( to_average is not None ):	
			for elt in to_average: 
				cpp_list.push_back((<PSSK>elt).thisptr)
			self.thisptr.compute_average_interface( cpp_list )   
                                 
	def compute_scalar_product(self, PSSK second):
		"""
		A function to compute scalar product of PSSK.
		The parameter of this function is a PSSK.
		This function is required in Topological_data_with_scalar_product concept.
		:param the heat map to compute scalar product with
		:type double
		"""
		if ( self.thisptr != NULL ) and ( second is not None ):
			return self.thisptr.compute_scalar_product( deref(second.thisptr) )
                        
	def get_x_range( self ):
		"""
		Returns the x range of PSSK.
		:type pair[double,double]
		"""
		if ( self.thisptr != NULL ):
			return self.thisptr.get_x_range() 
                        
	def get_y_range( self ):
		"""
		Returns the y range of PSSK.
		:type pair[double,double]
		"""
		if ( self.thisptr != NULL ):
			return self.thisptr.get_y_range()                               
