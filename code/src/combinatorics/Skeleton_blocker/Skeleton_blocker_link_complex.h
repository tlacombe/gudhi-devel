/*
 * Skeleton_blockers_link_complex.h
 *
 *  Created on: Feb 7, 2014
 *      Author: dsalinas
 */

#ifndef GUDHI_SKELETON_BLOCKERS_LINK_COMPLEX_H_
#define GUDHI_SKELETON_BLOCKERS_LINK_COMPLEX_H_

#include "utils/Utils.h"
#include "combinatorics/Skeleton_blocker/Skeleton_blocker_complex.h"



template<class ComplexType> class Skeleton_blocker_sub_complex;


/**
 *  \brief Class representing the link of a simplicial complex encoded by a skeleton/blockers pair.
 *  It inherits from Skeleton_blocker_sub_complex because such complex is a sub complex of a
 *  root complex.
 */
template<typename ComplexType>
class Skeleton_blocker_link_complex : public Skeleton_blocker_sub_complex<ComplexType>
{
	template <typename T> friend class Skeleton_blocker_link_superior;
	typedef typename ComplexType::Edge_handle Edge_handle;

	typedef typename ComplexType::boost_vertex_handle boost_vertex_handle;

	bool only_superior_vertices;

public:
	typedef typename ComplexType::Vertex_handle Vertex_handle;
	typedef typename ComplexType::Root_vertex_handle Root_vertex_handle;

	typedef typename ComplexType::Simplex_handle Simplex_handle;
	typedef typename ComplexType::Root_simplex_handle Root_simplex_handle;

	typedef typename ComplexType::BlockerMap BlockerMap;
	typedef typename ComplexType::BlockerPair BlockerPair;
	typedef typename ComplexType::BlockerMapIterator BlockerMapIterator;
	typedef typename ComplexType::BlockerMapConstIterator BlockerMapConstIterator;

	typedef typename ComplexType::Simplex_handle::const_iterator AddressSimplexConstIterator;
	typedef typename ComplexType::Root_simplex_handle::const_iterator IdSimplexConstIterator;

	typedef typename ComplexType::Complex_vertex_iterator Complex_vertex_iterator;



	Skeleton_blocker_link_complex(bool only_superior_vertices_=false):only_superior_vertices(only_superior_vertices_){
	}

	Skeleton_blocker_link_complex(const ComplexType & parent_complex, Simplex_handle& alpha_parent_adress,bool only_superior_vertices_ = false)
	:only_superior_vertices(only_superior_vertices_) {
		build_link(parent_complex,alpha_parent_adress);
	}

	Skeleton_blocker_link_complex(const ComplexType & parent_complex, Vertex_handle a_parent_adress, bool only_superior_vertices_ = false)
	:only_superior_vertices(only_superior_vertices_){
		only_superior_vertices = only_superior_vertices_ ;
		Simplex_handle alpha_simplex(a_parent_adress);
		build_link(parent_complex,alpha_simplex);
	}

	~Skeleton_blocker_link_complex(){
	}

protected:


	/**
	 * @brief compute vertices of the link.
	 * If the boolean only_superior_vertices is true, then only the vertices
	 * are greater than  vertices of alpha_parent_adress are added.
	 */
	void compute_link_vertices(const ComplexType & parent_complex, Simplex_handle& alpha_parent_adress,bool only_superior_vertices ){
		// we compute the intersection of neighbors of alpha and store it in link_vertices
		Simplex_handle link_vertices_parent;
		parent_complex.add_neighbours(alpha_parent_adress,link_vertices_parent,only_superior_vertices);
		// For all vertex 'v' in this intersection, we go through all its adjacent blockers.
		// If one blocker minus 'v' is included in alpha then the vertex is not in the link complex.
		for (auto v_parent : link_vertices_parent){
			bool new_vertex = true;
			for (BlockerMapConstIterator beta=parent_complex.blocker_map.lower_bound(v_parent);beta!=parent_complex.blocker_map.upper_bound(v_parent);++beta)
			{
				auto sigma_parent = beta->second;
				sigma_parent->remove_vertex(v_parent);
				if 	(alpha_parent_adress.contains(*sigma_parent)){
					new_vertex = false;
				}
				sigma_parent->add_vertex(v_parent);
				if(!new_vertex) break;
			}
			if (new_vertex) {
				this->add_vertex(parent_complex.get_id(v_parent));
			}
		}

	}

	void compute_link_edges(const ComplexType & parent_complex, Simplex_handle& alpha_parent_adress){
		AddressSimplexConstIterator y_link, x_parent , y_parent;
		// ----------------------------
		// Compute edges in the link
		// -------------------------
		for (auto x_link = this->vertex_range().begin() ;
				x_link!= this->vertex_range().end();
				++x_link
		)
		{
			auto y_link = x_link;
			for (++y_link ; y_link != this->vertex_range().end(); ++y_link){
				Vertex_handle x_parent = *parent_complex.get_address(this->get_id(*x_link));
				Vertex_handle y_parent = *parent_complex.get_address(this->get_id(*y_link));
				if (parent_complex.contains_edge(x_parent,y_parent) ){
					// we check that there is no blocker subset of alpha passing trough x and y
					bool new_edge=true;
					for (auto blocker_parent=parent_complex.blocker_map.lower_bound(x_parent);
							blocker_parent!=parent_complex.blocker_map.upper_bound(x_parent);
							++blocker_parent)
					{
						Simplex_handle* sigma_parent = (*blocker_parent).second;
						if 	(sigma_parent->contains(y_parent))
						{
							sigma_parent->remove_vertex(x_parent);
							sigma_parent->remove_vertex(y_parent);
							if(alpha_parent_adress.contains(*sigma_parent )){
								new_edge = false;
							}
							sigma_parent->add_vertex(x_parent);
							sigma_parent->add_vertex(y_parent);
							if (!new_edge) break;
						}
					}
					if (new_edge)
						this->add_edge(*x_link,*y_link);
				}
			}
		}
	}

	/**
	 * @brief : Given an address in the current complex, it returns the
	 * corresponding address in 'other_complex'.
	 * It assumes that other_complex have a vertex 'this.get_id(address)'
	 */
	boost::optional<Vertex_handle> give_equivalent_vertex(const ComplexType & other_complex,
			Vertex_handle address) const{
		Root_vertex_handle id((*this)[address].get_id());
		return other_complex.get_address(id);
	}

	void compute_link_blockers(const ComplexType & parent_complex,const Simplex_handle& alpha_parent){
		for (auto x_link : this->vertex_range()){

			Vertex_handle x_parent = * this->give_equivalent_vertex(parent_complex,x_link);

			for (BlockerMapConstIterator blocker_parent = parent_complex.blocker_map.lower_bound(x_parent); blocker_parent != parent_complex.blocker_map.upper_bound(x_parent); ++blocker_parent){

				Simplex_handle sigma_parent(*(*blocker_parent).second);
				sigma_parent.difference(alpha_parent);

				if (sigma_parent.dimension()>=2 && sigma_parent.first_vertex() == x_parent){
					Root_simplex_handle sigma_id(parent_complex.get_id(sigma_parent));
					auto  sigma_link = this->get_simplex_address(sigma_id);
					// ie if the vertices of sigma are vertices of the link
					if(sigma_link){
						bool is_new_blocker = true;
						for(auto a : alpha_parent){
							for(auto eta_parent : parent_complex.const_blocker_range(a)){
								Simplex_handle eta_minus_alpha(*eta_parent);
								eta_minus_alpha.difference(alpha_parent);
								if(eta_minus_alpha != sigma_parent && sigma_parent.contains(eta_minus_alpha)){
									is_new_blocker = false;
									break;
								}
							}
							if (!is_new_blocker) break;
						}
						if (is_new_blocker)
							this->add_blocker(new Simplex_handle(*sigma_link));

					}

				}
			}
		}
	}


public:
	/**
	 * @brief compute vertices, edges and blockers of the link.
	 * If the boolean only_superior_vertices is true, then the link is computed only
	 * with vertices that are greater than  vertices of alpha_parent_adress.
	 */
	void build_link(const ComplexType & parent_complex, Simplex_handle& alpha_parent_adress)
	{
		compute_link_vertices(parent_complex,alpha_parent_adress,only_superior_vertices);
		compute_link_edges(parent_complex,alpha_parent_adress);
		compute_link_blockers(parent_complex,alpha_parent_adress);
	}




};

#endif /* SKELETON_BLOCKERS_LINK_COMPLEX_H_ */
