Syntax for (modified) off file :

OFF	# Header keyword
NVertices  NFaces  NEdges   # NEdges not used or checked
     
pt[0][0] ... pt[0][dim] # Coordinate of point of dimension 'dim'

...

pt[NVertices-1][0] ... pt[NVertices-1][dim] # Coordinate of last point

# plus eventually some faces      
# Nv = # vertices on this face
# v[0] ... v[Nv-1]: vertex indices
#		in range 0..NVertices-1
Nv  v[0] v[1] ... v[Nv-1]  # a face Nv-dimensional face 
     ...
         			
