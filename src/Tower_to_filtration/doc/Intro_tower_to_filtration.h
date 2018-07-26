/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2018  INRIA Sophia Antipolis-Méditerranée (France)
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

#ifndef DOC_TOWER_TO_FILTRATION_H_
#define DOC_TOWER_TO_FILTRATION_H_

namespace Gudhi {
namespace tower_to_filtration {

/** \defgroup tower_to_filtration Tower to Filtration
 *
 * \author    Hannah Schreiber
 *
 * @{
 *
 * The module is based on @cite KerberS17. It transforms a tower into a filtration with the same barcode.
 *
 * \section background Background
 *
 * A @a tower of length @f$m@f$ is a collection of simplicial complexes @f$\mathbb{K}_0,\ldots,\mathbb{K}_m@f$
 * and simplicial maps @f$\phi_i:\mathbb{K}_i\rightarrow\mathbb{K}_{i+1}@f$ for @f$i=0,\ldots,m-1@f$:
 * @f[
 *      \mathbb{K}_0 \xrightarrow{\phi_0} \mathbb{K}_1 \xrightarrow{\phi_1} \cdots \xrightarrow{\phi_{m-1}} \mathbb{K}_m.
 * @f]
 * 
 * A map @f$\mathbb{K}\stackrel{\phi}{\rightarrow}\mathbb{L}@f$ between simplicial complexes is called @a simplicial
 * if with @f$\sigma=\{v_0,\ldots,v_k\}\in\mathbb{K}@f$, @f$\phi(\sigma)@f$ is equal to
 * @f$\{\phi(v_0),\ldots,\phi(v_k)\}@f$ and @f$\phi(\sigma)@f$ is a simplex in @f$\mathbb{L}@f$.
 * By definition, a simplicial map maps vertices to vertices and is completely determined by its action on the vertices.
 * Moreover, the composition of simplicial maps is again simplicial.
 *
 * A simple example of a simplicial map is the inclusion map @f$\mathbb{L}\stackrel{\phi}{\hookrightarrow}\mathbb{K}@f$
 * where @f$\mathbb{L}@f$ is a subcomplex of @f$\mathbb{K}@f$. If @f$\mathbb{K}=\mathbb{L}\cup\{\sigma\}@f$
 * with @f$\sigma\notin\mathbb{L}@f$, we call @f$\phi@f$ an @a elementary @a inclusion.
 * The simplest example of a non-inclusion simplicial map is @f$\mathbb{K}\stackrel{\phi}{\rightarrow}\mathbb{L}@f$
 * such that there exist two vertices @f$u,v\in\mathbb{K}@f$ with
 * @f$\mathcal{V}(\mathbb{L})=\mathcal{V}(\mathbb{K})\setminus\{v\}@f$,
 * @f$\phi(u)=\phi(v)=u@f$, and @f$\phi@f$ is the identity on all remaining vertices of @f$\mathbb{K}@f$.
 * We call @f$\phi@f$ an @a elementary @a contraction.
 * These notions were introduced by Dey, Fan and Wang in @cite DBLP:journals/corr/abs-1208-5018 and they also showed that
 * any simplicial map @f$\mathbb{K}\stackrel{\phi}{\rightarrow}\mathbb{L}@f$ can be written as the composition of
 * elementary contractions and inclusions.
 *
 * A tower is called a @a filtration if all @f$\phi_i@f$ are inclusion maps.
 *
 * TODO: Definition of Complex and Barcode ?
 *
 * \section Usage
 *
 * TODO
 * add_inclusion / add_contraction functions
 * >> and << functions
 * file format
 *
 * To come:
 *
 * \section Examples
 *
 * TODO
 *
 * @}
 */

}
}

#endif  // DOC_TOWER_TO_FILTRATION_H_
