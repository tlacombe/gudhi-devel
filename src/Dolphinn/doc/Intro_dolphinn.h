#ifndef DOC_DOLPHINN_INTRO_DOLPHINN_H_
#define DOC_DOLPHINN_INTRO_DOLPHINN_H_

// needs namespace for Doxygen to link on classes
namespace Gudhi {
// needs namespace for Doxygen to link on classes
namespace dolphinn {

/**  \defgroup dolphinn Dolphinn
 * 
 * \author    Owen Rouill&eacute;, Georgios Samaras
 * @{
 * 
 * \section dolphinndefinition Definition
 * 
 * Dolphinn is a <a target="_blank" href=https://en.wikipedia.org/wiki/Locality-sensitive_hashing>LSH</a>-based method used for approximate neighbor search (k-nearest neigbors and radius search).
 * The data structure used is a unit hypercube of dimension D which stores the data points on its vertices.
 * The vertex on which a point is stored is determined by the result of D LSH functions (or LSH based functions),
 * each determining a coordinate of the vertex (0 or 1 for each coordinate). Building the data structure is linear 
 * in time and space in both dimension and number of points. As it is LSH based, the method gives approached results.
 *
 * The LSH functions currently implemented are the <a target="_blank" href=https://en.wikipedia.org/wiki/Locality-sensitive_hashing#Stable_distributions>stable distribution</a>
 * and the hyperplane hashing.
 *
 * More details can be found in <a target="_blank" href=https://arxiv.org/abs/1612.07405>this paper</a>.
 *
 */
/** @} */  // end defgroup dolphinn

}  // namespace dolphinn

}  // namespace Gudhi

#endif  // DOC_DOLPHINN_INTRO_DOLPHINN_H_
