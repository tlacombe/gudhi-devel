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

/** @defgroup tower_to_filtration Tower to Filtration
 *
 * @author    Hannah Schreiber
 *
 * @{
 *
 * The module is based on @cite KerberS17. It transforms a tower into a filtration with the same barcode.
 *
 * @section background Background
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
 * The primary aim of this module is to compute a filtration @f$\mathcal{F}@f$ from a given tower @f$\mathcal{T}@f$,
 * such that @f$\mathcal{F}@f$ and @f$\mathcal{T}@f$ have the same barcode, using @ref Tower_converter.
 * It can also compute the corresponding barcode using @ref Persistence.
 *
 * <span style="color:red;">TODO</span>: Definition of Complex and Barcode<span style="color:red;">?</span>
 *
 * @section usage Usage
 *
 * @ref Tower_converter and @ref Persistence can both be used via two functions:
 *
 * <span style="color:red;">TODO</span>: make @ref Persistence independent of @ref Tower_converter
 * (i.e. @ref Tower_converter has to be used externally).
 *
 * - @ref Tower_converter<ComplexStructure>::add_insertion /
 *	@ref Persistence<ComplexStructure,ColumnType>::add_insertion
 *
 * Add an elementary tower insertion operation, which is directly processed.
 * As result, the insertion is added to the output filtration.
 *
 * - @ref Tower_converter<ComplexStructure>::add_contraction /
 *	@ref Persistence<ComplexStructure,ColumnType>::add_contraction
 *
 * Add an elementary tower contraction operation, which is directly processed.
 * As result, all the sequence of insertions equivalent to the given contraction  is added to the output filtration.
 *
 * - <span style="color:red;">TODO</span>
 *
 * Allowing sequence of edges as input to build flag complex
 *
 * @subsection temp Templates
 *
 * @ref Tower_converter uses one template @ref ComplexStructure as complex type.
 * For now only Gudhi::tower_to_filtration::Hash_complex (hash_complex.h) and Gudhi::tower_to_filtration::Simplex_tree (simplex_tree.h)
 * satisfy the concept. <span style="color:red;">(To be rectified)</span>
 *
 * @ref Persistence additionnaly uses a template @ref ColumnType
 * for the columns of its internally stored boundary matrix.
 * @ref Heap_column and @ref List_column satisfy the concept.
 * The first one is based on a heap representation of the column, the second one on a list representation.
 *
 * @subsection formats Input file and output formats
 *
 * - Output format
 *
 * @ref Tower_converter outputs its result either in a file,
 * or in an output stream (std::stringstream). In both cases, it uses the same output format:
 *
 * Each new line corresponds to an inclusion of a @f$d@f$-simplex @f$s@f$:
 *
 * @f$d@f$ @f$v_1@f$ ... @f$v_{d+1}@f$ @f$ts@f$,
 *
 * where @f$(v_i)_{1 \leq i \leq d+1}@f$ is by default the set of vertices of @f$s@f$
 * (option @ref Tower_converter::VERTICES)
 * or the set of facets of @f$s@f$ (option @ref Tower_converter::FACES).
 * For each @f$i < j \in \{1,...,d+1\}@f$, @f$v_i < v_j@f$.
 * And @f$ts@f$ corresponds to the filtration value of @f$s@f$.
 *
 * @ref Persistence outputs the persistence barcode in a file,
 * where each new line corresponds to a persistence pair @f$(b,d)@f$ of dimension @f$dim@f$:
 *
 * @f$dim@f$ @f$b@f$ @f$d@f$.
 *
 * Essential cycles (i.e. paired with infinity) are not printed.
 * (<span style="color:red;">TODO</span> : add "finilize function" to enable the printing of those paires<span style="color:red;">?</span>)
 *
 * - Input format
 *
 * The file tc_reading_utilities.h includes two reading functions
 * `read_operation` and `>>`
 * (see @subpage Tower_to_filtration/example_tower_from_file_write_filtration_into_file.cpp and
 * @subpage Tower_to_filtration/example_tower_from_file_streamed_output.cpp
 * for examples of use).
 * The input format for both is the following.
 *
 * Each line is either a comment beginning with '#' or a tower operation:
 *
 * In a case of an inclusion of a @f$d@f$-simplex @f$s@f$:
 *
 * [@f$ts@f$] i @f$v_1@f$ ... @f$v_{d+1}@f$,
 *
 * where @f$(v_i)_{1 \leq i \leq d+1}@f$ is the set of vertices of @f$s@f$
 * and for each @f$i < j \in \{1, ..., d+1\}@f$, @f$v_i < v_j@f$.
 * And @f$ts@f$ is an optional time indicator.
 * (If there is no time indication, time will start at 0 and increase by one at each insertion or contraction
 * when using `>>`).
 *
 * In a case of a contraction:
 *
 * [@f$ts@f$] c @f$v_d@f$ @f$v_k@f$,
 *
 * where @f$v_d@f$ and @f$v_k@f$ are the vertices to be contracted.
 * From here on, the remaining vertex needs to be refered by @f$v_k@f$ and **NOT** @f$v_d@f$.
 * And @f$ts@f$ is again the optional time indicator.
 *
 * @section examples Examples
 *
 * Following examples are avaible in the 'example/Tower_to_filtration/' folder.
 *
 * - @subpage Tower_to_filtration/example_tower_from_file_write_filtration_into_file.cpp
 *
 * @code{.sh}
 *	./Tower_to_filtration_example_tower_from_file_write_filtration_into_file input_file_name output_file_name
 * @endcode
 * Simple example of how to use the module when reading the tower from a file
 * and writing the resulting filtration in a file.
 * See @ref formats.
 *
 * - @subpage Tower_to_filtration/example_tower_from_file_streamed_output.cpp
 *
 * @code{.sh}
 *	./Tower_to_filtration_example_tower_from_file_streamed_output input_file_name
 * @endcode
 * Simple example of how to use the module when reading the tower from a file and
 * writing the resulting filtration in an output stream.
 * See @ref formats.
 * The output stream should be processed/emptied regularly to avoid memory consummption.
 *
 * - @subpage Tower_to_filtration/example_elementary_input_streamed_output.cpp
 *
 * @code{.sh}
 *	./Tower_to_filtration_example_elementary_input_streamed_output
 * @endcode
 * Simple example of how to add tower operations to the module and
 * writing the resulting filtration in an output stream.
 * See @ref formats for stream output format.
 * The output stream should be processed/emptied regularly to avoid memory consummption.
 *
 * - @subpage Tower_to_filtration/example_elementary_input_write_filtration_into_file.cpp
 *
 * @code{.sh}
 *	./Tower_to_filtration_example_elementary_input_write_filtration_into_file output_file_name
 * @endcode
 * Simple example of how to add tower operations to the module
 * and writing the resulting filtration in a file.
 * See @ref formats for output format.
 *
 * - <span style="color:red;">TODO</span>: examples for persistence
 *
 * @}
 */

}
}

#endif  // DOC_TOWER_TO_FILTRATION_H_
