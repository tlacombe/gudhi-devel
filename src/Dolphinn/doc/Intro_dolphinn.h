#ifndef DOC_DOLPHINN_INTRO_DOLPHINN_H_
#define DOC_DOLPHINN_INTRO_DOLPHINN_H_

// needs namespace for Doxygen to link on classes
namespace Gudhi {
// needs namespace for Doxygen to link on classes
namespace dolphinn {

/**  \defgroup dolphinn Dolphinn
 
 \author    Owen Rouill&eacute;, Georgios Samaras
 @{
 
\section dolphinndefinition Definition

Dolphinn is a <a target="_blank" href=https://en.wikipedia.org/wiki/Locality-sensitive_hashing>LSH</a>-based method used for approximate neighbor search (k-nearest neigbors and radius search).
The data structure used is a unit hypercube of dimension  \f$ K \f$  which stores the data points on its vertices.
The vertex on which a point is stored is determined by the result of  \f$ K \f$  LSH functions (or LSH based functions),
each determining a coordinate of the vertex (0 or 1 for each coordinate). Building the data structure is linear 
in time and space in both dimension and number of points. As it is LSH based, the method gives approached results.

The LSH functions currently implemented are the stable distribution and the hyperplane hashing. More details can be found in <a target="_blank" href=https://arxiv.org/abs/1612.07405>this paper</a> (hyperplane hashing was not implmented at the time).


\section lshsec Locality-Sensitive Hashing

\subsection lshdef Definition

An LSH function is a hashing function that tends to send points that are close to each other in the same bucket and points that are far from each other to different ones with good probabilities. They can be viewed as random partitions of the space. 

Let  \f$ (M,d) \f$  be a metric space,  \f$ 0<r_1<r_2 \f$  and  \f$ 0 \leq p_1 < p_2 \leq 1 \f$  be reals, a family  \f$ F \f$  of functions is said  \f$ (p_1,p_2,r_1,r_2) \f$ -sensitive if for all  \f$ (x,y)\in M \f$  and a  \f$ h\in F \f$  chosen uniformly at random:  \f$ d(x,y)\leq r_1 \implies P[h(x)=h(y)]\geq p_1 \f$  and  \f$ d(x,y)\geq r_2 \implies P[h(x)=h(y)]\leq p_2 \f$ .

There exists several such families. Here we are concerned by two of them: the line LSH and the hyperplane LSH.

\subsection lshline Line LSH

The stable distribution based LSH, called here "line LSH", is parametrized by a real  \f$ r \f$  (different  \f$ r \f$  give different subfamilies), they were designed to work with the Euclidean distance. Let  \f$ (M,d) \f$  be an Euclidean space,  \f$ r \f$  a real,  \f$ a \in M \f$  a vector with entries chosen independently from a stable distribution and  \f$ b \in [ 0;r ]  \f$  chosen uniform at random,  \f$  h_{a,b}:~x\in M \mapsto \lfloor \frac{<a,x> + b}{r} \rfloor  \f$  is a member of this family. This function consists in orthogonally projecting  \f$ x \f$  on the line passing by the origin and directed by  \f$ a \f$ , shifting the result by a constant quantity, and returning the part of the line the result landed in. The line is sliced in windows of size  \f$ r \f$ , each one of these windows is a possible output for the hashing function.

\image latex "linelsh.eps" "Illustration of the line LSH" width=10cm

In the figure above, the small spheres are elements of the dataset, the double lined line separated in windows of the same size ( \f$ r \f$  in the text) is the support of the projection, two points in the same windows (delimited by the dashed lines) will have the same image.

\subsection lshhyper Hyperplane LSH

The hyperplane LSH is a family initially designed for cosine similarity, it offers good properties for a dataset uniformly distributed on a sphere. The idea is to create a random hyperplane passing by the origin and the points that are on the same side of the hyperplane have the same image through the function. To do so, a vector with independent and identically distributed Gaussian coordinates is generated, this is the normal vector of the hyperplane, and the sign of the dot product with the point to hash gives the side of the hyperplane. 

This LSH family is not based on a distance and does not require any parameter making it easier to use. However, this is not a perfect family as the datasets are rarely uniformly distributed on a sphere. To avoid pathologic cases where all the points are on the same side of an hyperplane, the dataset is centered before the hashing. Instead of rewriting all the points, the centering is done by computing the center of the dataset and subtracting this vector from the points before the dot product. The formula obtained is for any point  \f$ p \f$ :  \f$ sgn(<p-c,v>) \f$  where  \f$ sgn \f$  is the function that returns 1 for any non-negative real and 0 otherwise,  \f$ c \f$  is the center of the dataset, and  \f$ v \f$  the normal vector of the hyperplane.

\image latex "hyplsh.eps" "Illustration of the line LSH, the arrow its normal vector, and the two types of points the two groups obtained." width=10cm

\section doldef Dolphinn

Dolphinn is a LSH based method that consists in using a vector of LSH functions to project the points on the vertices of an hypercube such that two points in the original space will be sent to neighboring vertices of the hypercube. It has three parameters: the choice of the LSH family, the dimension of the hypercube on which the points are sent, and a parameter that corresponds to a timeout.

\subsection dolprinc Principle

Let  \f$ K \f$  be a positive integer, and  \f$ \mathbb{F} \f$  an LSH family.

let  \f$ [h_1,...,h_K] \f$  be a vector of  \f$ K \f$  LSH functions of  \f$ \mathbb{F} \f$ . If the codomain of this LSH family is of size greater than  \f$ 2 \f$ : for each  \f$ h_i \f$ , let  \f$ f_i \f$  be a function that assign every possible output of  \f$ h_i \f$  to  \f$ 0 \f$  or  \f$ 1 \f$  uniformly (for each possible output  \f$ y \f$  of  \f$ h_i \f$ , flip a coin, if head  \f$ f_i(y)=0 \f$ ,   \f$ f_i(y)=1 \f$  otherwise), and let's redefine  \f$ h_i=f_i \circ h_i \f$ . 

Now each  \f$ h_i \f$  maps points into  \f$ \{0,1\} \f$  and  \f$ [h_1,...,h_K] \f$  into  \f$ \{0,1\}^K \f$  (which can be seen as the vertices of an hypercube of dimension  \f$ K \f$ ). Dolphinn is based on the observation that if two points have close image through the vector, it means several LSH families mapped them together, and thus the probability of the being close to each other is high.

Dolphinn's data structure is obtained by hashing all the points of the dataset and storing them on the vertices of an hypercube. The following figure illustrates the construction of an hypercube for a regular dataset in dimension 2, with  \f$ K=2 \f$  and a line LSH. However, for simplicity, the lines used for the LSH are perpendicular, which does not happen in practice. This building is linear in time and space in the size of the dataset, its dimension, and  \f$ K \f$ .

\image latex dol2.eps "Construction of the hypercube" width=10cm

The left part of the figure represents a regular dataset (with large dots) and two lines for the hashing, the light grey dashed lines corresponds to the windows of the LSH functions, the  \f$ 0 \f$  and  \f$ 1 \f$  corresponds to the result of the coin flip for each window. The right part represents the resulting hypercube of dimension 2 (a square), with the coordinates of each on of its vertices. The circled data points got the same results through the two hashing functions: they are sent to the same vertex (here  \f$ (1,1) \f$ ).

To query the data structure, the queried point is hashed with the vector of functions which gives a vertex  \f$ v \f$ . The hypercube is then searched, starting from  \f$ v \f$ , and following by the vertices at Hamming distance 1, 2... until the timeout is reached (or a witness is found for the range search). The timeout mechanism in Dolphinn is not a time limit, but a maximal number of points compared to the query (points of the original dataset, not the hypercube's vertices).

NB: the Hamming distance between two vectors is the number of differing coordinates.

\subsection dolparam Choosing the parameters

The LSH family: this parameter depends on the dataset, the hyperplane LSH has the advantage of not having to choose the size of the windows. Concerning the value of  \f$ r \f$ , if too small, the partition proposed by Dolphinn is close to random, if too large, the partition can be too uneven, up to the case where all the points lie in the same window. The optimal value of  \f$ r \f$  depends on the dataset. The following figure gives the accuracy obtained for the two categories of LSH families against the value of  \f$ r \f$ . Here the accuracy is the average percentage of 10-nearest neighbors found with the other parameters fixed. This graph was made with a single dataset and it is not guarantied that the hyperplane LSH has always better results than the line LSH.

\image html "./acc.png" "Accuracy comparison for line and hyperplane LSH"

Plot for the accuracy against the size of the windows for the line LSH (log scale), compared to the accuracy obtained with the hyperplanes. The optimal value for  \f$ r \f$  is around  \f$ 700 \f$  for an accuracy of  \f$ 50\% \f$ , which is significantly smaller than the result of the hyperplanes:  \f$ 70\% \f$ . For too large and too small values of  \f$ r \f$ , the results tend to be no better than by trying random points. These results where obtained with the small version of ANN\_SIFT1M (10000 points), a hypercube dimension of 14 and a maximum points to explore of 1000.

The Dimension of the hypercube and the timeout: these two parameters increases both the query time and the accuracy. There is no universal best solution. For  \f$ K \f$  fixed, the number of vertices that are at distance  \f$ d \f$  from a given one is  \f$ \binom{K}{d} \f$ , this means for a large  \f$ K \f$  finding vertices has an exponential cost. There is a tradeoff to find between doing a lot of distance computation (big "timeout") and spending a lot a time searching for new vertices of the hypercube (large  \f$ K \f$ ). Unfortunately Dolphinn suffers from accuracy problems, which means we would like to use greater  \f$ K \f$ . The search of the vertices is still worked on to be able to handle these big  \f$ K \f$ .

Starting values: to begin with Dolphinn it is advised to start with  \f$ K \f$  of the same order as the  \f$ log_2 \f$  of the dataset's size and with the hyperplane LSH.

\section dolhowto How to use Dolphinn?

Dolphinn Works with array of points, points being array of doubles. The first step to use Dolphinn is to fill the hypercube with hte initial dataset, this is done by creating an object of the Dolphinn class. To find the neighbors of a set of points, these have to be put together in an array, and the correponding method called (either k_nearest_neighbors or radius query). These methods fill an array given to them by reference. The following code is an example of use of Dolphinn for the search of the k-nearest neighbors.

\include Dolphinn/Dolphinn_example_knn_from_points.cpp

 */
/** @} */  // end defgroup dolphinn

}  // namespace dolphinn

}  // namespace Gudhi

#endif  // DOC_DOLPHINN_INTRO_DOLPHINN_H_
