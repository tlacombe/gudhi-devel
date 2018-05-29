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

   Author(s):	   Pawel Dlotko

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



cdef extern from "Persistence_vectors_interface.h" namespace "Gudhi::Persistence_representations":
	cdef cppclass Vector_distances_in_diagram_interface "Gudhi::Persistence_representations::Vector_distances_in_diagram_interface":
		Vector_distances_in_diagram_interface()
		Vector_distances_in_diagram_interface(vector[pair[double, double]], size_t where_to_cut ) 
		Vector_distances_in_diagram_interface(const char* filename, size_t where_to_cut, int dimension)				
		size_t size()const				
		void print_to_file(const char* filename)const
		void load_from_file(const char* filename)const		
		double project_to_R(size_t number_of_function)const			   
		size_t number_of_projections_to_R()const		
		vector[double] vectorize(int number_of_function) const
		size_t number_of_vectorize_functions() const 
		void compute_average_interface(const vector[Vector_distances_in_diagram_interface*] to_average)
		double distance(const Vector_distances_in_diagram_interface& second, double power)
		double compute_scalar_product(const Vector_distances_in_diagram_interface& second)
	 
		
		
cdef class PersistenceVectors:

	cdef Vector_distances_in_diagram_interface* thisptr

#Can we have only one constructor, or can we have more
	def __init__(self, where_to_cut=100, file_with_intervals='', 
	vector_of_intervals=None, dimension=-1):
		"""
		This is an implementation of idea presented in the paper 
		<i>Stable Topological Signatures for Points on 3D  Shapes</i> 
		\cite Carriere_Oudot_Ovsjanikov_top_signatures_3d .<br>
		The parameter of the class is the class that computes distance 
		used to construct the vectors. The typical function is either 
		Euclidean of maximum (Manhattan) distance.
 
		This class implements the following concepts: 
		Vectorized_topological_data, Topological_data_with_distances,
		Real_valued_topological_data, Topological_data_with_averages, 
		Topological_data_with_scalar_product
		"""	

	def __cinit__(self, where_to_cut=100, file_with_intervals='', 
	vector_of_intervals=None, dimension=-1):
		"""
		This is a constructor of a class PersistenceVectors.
		It either take text file and a positive integer, or a vector
		of pairs and a positive integer. In case of file, each line of 
		the input file is supposed to contain two numbers of a type double
		(or convertible to double) representing the birth and the death
		of the persistence interval. If the pairs are not sorted so that
		birth <= death, then the constructor will sort then that way.
		In case of vector of pairs, it simply accept vector of pair of
		doubles.
		:param vector_of_intervals -- vector of pairs of doubles with
			   birth-death pairs. None if we construct it from file.
		:type vector of pairs of doubles or None
		:param dimension -- diension of intervals to be extracted from file
		:type nonnegative integer or None
		:param file_with_intervals - a path to Gudhi style file with
			   persistence interfals.	   
		:param where_to_cut - number of elements to be generated in the 
		vector. 
		"""
		if ( (vector_of_intervals is None) and ( file_with_intervals is '' ) ):
			print "Please provide parameters to construct the persistence vectors." 
		else:
			if (vector_of_intervals is None) and (file_with_intervals is not ''):
				if os.path.isfile(file_with_intervals):
					self.thisptr = new Vector_distances_in_diagram_interface(str.encode(file_with_intervals), where_to_cut, dimension)
				else:
					print("file " + file_with_intervals + " not found.")
			else:			
				if (file_with_intervals is '') and (vector_of_intervals is not None):
					self.thisptr = new Vector_distances_in_diagram_interface(vector_of_intervals, where_to_cut)														
				else:
					self.thisptr = new Vector_distances_in_diagram_interface()


	def __dealloc__(self):
		"""
		destructor
		"""
		if self.thisptr != NULL:
			del self.thisptr	  

	def load_vector_from_file(self,filename):
		"""
		This procedure loads a persistence vector from file. It erase all the data
		that was previously stored in this vector.
		:param Name of the file.
		:type String
		"""
		if ( self.thisptr != NULL ) and ( filename is not None ):
			self.thisptr.load_from_file(filename)
			

	def print_to_file(self,filename) :
		"""
		The procedure stores a persistence vector to a file. The file can be later
		used by a procedure load_vector_from_file.
		:param Name of the file.
		:type String
		"""
		if ( self.thisptr != NULL ) and ( filename is not None ):
			self.thisptr.print_to_file(filename)

	def size( self ):
		"""
		Returns the size of the vector.
		"""
		if ( self.thisptr != NULL ):
			return self.thisptr.size()
  
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
		:type nonnegative integer
		"""
		if ( self.thisptr != NULL ) and ( number_of_function is not None ):
			return self.thisptr.project_to_R(number_of_function)


	def number_of_projections_to_R(self):
		"""
		The function gives the number of possible projections to R. This
		function is required by the
		Real_valued_topological_data concept
		"""
		if ( self.thisptr != NULL ):
			return self.thisptr.number_of_projections_to_R()

	def vectorize(self, number_of_function):
		"""
		This function returns a vector. It
		is required in a concept
		Vectorized_topological_data
		:param number of function
		:type nonnegative intege
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

	def compute_average( self,to_average=[] ):
		"""
		A function to compute averaged persistence vector, based on vector
		of persistence vectors.
		This function is required by Topological_data_with_averages concept.
		:param vector of persistence vectors to average
		:type vectors of references to persisence vectors. 
		"""			
		cdef vector[Vector_distances_in_diagram_interface*] cpp_list	
		if ( self.thisptr != NULL ) and ( to_average is not None ):	
			for elt in to_average: 
				cpp_list.push_back((<PersistenceVectors>elt).thisptr)
			self.thisptr.compute_average_interface( cpp_list )

			


	def distance(self, PersistenceVectors second, power):
		"""
		A function to compute distance between persistence vectirs.
		The parameter of this function is a Persistence_vector.
		This function is required in Topological_data_with_distances concept.
		For max norm distance, set power to numeric_limits<double>::max()
		:param the landascape to compute distance to
		:type PersistenceLandscape
		"""
		if ( self.thisptr != NULL ) and ( second is not None ) and ( power is not None ):
			return self.thisptr.distance(deref(second.thisptr), power)
		else:
			return self.thisptr.distance( deref(second.thisptr) ,0 )
 
	def compute_scalar_product(self, PersistenceVectors second):
		"""
		A function to compute scalar product of persistence vectors.
		The parameter of this function is a Persistence_vector.
		This function is required in Topological_data_with_scalar_product concept.
		:param the landascape to compute scalar product with
		:type PersistenceLandscape
		"""
		if ( self.thisptr != NULL ) and ( second is not None ):
			return self.thisptr.compute_scalar_product( deref(second.thisptr) )

 

