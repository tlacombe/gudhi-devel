
/*! \mainpage Skeleton blockers data-structure

\author David Salinas

\section Introduction
The Skeleton-blocker data-structure had been introduced in the two papers 
\cite skbl_socg2011 \cite skbl_ijcga2012. 
It proposes a light encoding for simplicial complexes by storing only an *implicit* representation of its
simplices.
Intuitively, it just stores the 1-skeleton of a simplicial complex with a graph and the set of its "missing faces" that
is very small in practice (see next section for a formal definition).
This data-structure handles every classical operations used for simplicial complexes such as
 as simplex enumeration or simplex removal but operations that are particularly efficient 
 are operations that do not require simplex enumeration such as edge iteration, link computation or simplex contraction.


\todo{image}

\section Definitions

\section API

\subsection Overview

Four classes define  a simplicial complex namely :

\li <Code>Skeleton_blocker_complex</Code> : a simplicial complex with basic operations such as vertex/edge/simplex enumeration and construction
\li <Code>Skeleton_blocker_sub_complex</Code> : a sub complex that is used for representing a link of a simplex in a parent complex
\li <Code>Skeleton_blocker_simplifiable_complex</Code> : a simplicial complex with simplification operations such as edge contraction or simplex collapse
\li <Code>Skeleton_blocker_geometric_complex</Code> : a simplicial complex who has access to geometric points in  \f$R^d\f$ 

The two last classes are derived classes from the <Code>Skeleton_blocker_complex</Code> class. The class <Code>Skeleton_blocker_sub_complex</Code> inheritates from a template passed parameter
that may be either <Code>Skeleton_blocker_complex</Code> or <Code>Skeleton_blocker_geometric_complex</Code> (a sub-complex may have points or not).

\todo{include links}

\subsection Visitor

The class <Code>Skeleton_blocker_complex</Code> has a visitor that is called when usual operations such as adding an edge or remove a vertex are called.
You may want to use this visitor to compute statistics or to update another data-structure (for instance this visitor is heavily used in the 
<Code>Contraction</Code> package).



\section Example

 
\subsection s Iterating through vertices, edges, blockers and simplices	

  \code{.cpp}
  typedef Skeleton_blocker_complex<T> Complex;
  Complex complex;
  for(int i = 0 ; i < 10; ++i);
  	complex.add_vertex();
  for(auto v : complex.vertex_range())
  	std::cout << "Vertex "<<v<<std::endl;
  for(auto e : complex.edge_range())
  	std::cout << "Edge "<<e<<std::endl;  	
  \endcode

 
Iteration through simplices is straightforward with c++11 for range loops.
Note that simplex iteration with this implicit data-structure just takes
a few more time than iteration via an explicit representation 
such as the Simplex Tree. The following example computes the Euler Characteristic
of a simplicial complex.

  \code{.cpp}
  typedef Skeleton_blocker_complex<T> Complex;
  Complex complex;

  // build a full complex with 10 vertices and 2^10-1 simplices
  for(int i = 0 ; i < 10; ++i);
  	complex.add_vertex();
  for(int i = 0 ; i < 10; ++i);
  	for(int j = i+1 ; j < 10; ++j);
  	complex.add_edge(Vertex_handle(i),Vertex_handle(j));
  
  unsigned euler = 0;
  // we use a reference to a simplex instead of a copy
  // value here because a simplex is a set of integers 
  // and copying it cost time
  for(const auto & s : complex.simplex_range()){
	if(s.dimension()%2 == 0) euler += 1;
	else euler -= 1;
  }
  \endcode


\verbatim
This is verbatim
\endverbatim



\subsection Acknowledgements
The author wishes to thank Dominique Attali for leaving him use a prototype 
of the data-structure. 


\copyright GNU General Public License v3.                         
\verbatim  Contact: David Salinas,     david.salinas@inria.fr \endverbatim

*/
