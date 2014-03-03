#ifndef GUDHI_SKELETON_BLOCKER_SUB_COMPLEX_H
#define GUDHI_SKELETON_BLOCKER_SUB_COMPLEX_H

#include "Skeleton_blocker_complex.h"
#include "Skeleton_blocker_simplex.h"
#include "utils/Utils.h"

/**
 * @brief Simplicial subcomplex of a complex represented by a skeleton/blockers pair.
 *
 * @details Stores a subcomplex of a simplicial complex.
 * To simplify explanations below, we will suppose that :
 * - K is the root simplicial complex
 * - L is a subcomplex of K.
 *
 * One vertex of K may exists in L but with a different address.
 * To be able to locate the vertices in K from vertices of L, the class
 * stores a map 'adresses' between vertices of K and vertices of L.
 *
 * Note that the type for handle of vertices of L is 'Vertex_handle' and
 * the type for handle of vertices of K is 'Root_vertex_handle'.
 *
 * The template ComplexType is type of the root complex. It allows to know
 * if the subcomplex is geometric or not.
 * It has to be either 'Skeleton_blockers_complex' or 'Skeleton_blockers_geometric_complex'.
 *
 */
template<typename ComplexType>
class Skeleton_blocker_sub_complex : public ComplexType
{

public:
	typedef typename ComplexType::Graph Graph;
	typedef typename ComplexType::Edge_handle Edge_handle;

	typedef typename ComplexType::boost_vertex_handle boost_vertex_handle;

	typedef typename ComplexType::Vertex_handle Vertex_handle;
	typedef typename ComplexType::Root_vertex_handle Root_vertex_handle;
	typedef typename ComplexType::Simplex_handle Simplex_handle;
	typedef typename ComplexType::Root_simplex_handle Root_simplex_handle;

	template<class T> friend class Skeleton_blocker_link_complex;


	/**
	 * @brief Returns true iff the simplex formed by all vertices contained in 'addresses_sigma_in_link'
	 * but 'vertex_to_be_ignored' is in 'link'
	 */
	template<typename T> friend bool
	proper_face_in_union(
			Skeleton_blocker_sub_complex<T> & link,
			std::vector<boost::optional<typename T::Vertex_handle> > & addresses_sigma_in_link,
			int vertex_to_be_ignored);

	/**
	 * @brief Determines whether all proper faces of simplex 'sigma' belong to 'link1' \cup 'link2'
	 * where 'link1' and 'link2' are subcomplexes of the same complex of type ComplexType
	 */
	template<typename T> friend	bool
	proper_faces_in_union(Skeleton_blocker_simplex<typename T::Root_vertex_handle> & sigma, Skeleton_blocker_sub_complex<T> & link1, Skeleton_blocker_sub_complex<T> & link2);


protected:
	typedef std::map<Root_vertex_handle,Vertex_handle> IdAddressMap;
	typedef typename IdAddressMap::value_type AddressPair;
	typedef typename IdAddressMap::iterator IdAddressMapIterator;
	typedef typename IdAddressMap::const_iterator IdAddressMapConstIterator;
	std::map<Root_vertex_handle,Vertex_handle> adresses;

	/**
	 * Add a vertex 'global' of K to L. When added to L, this vertex will receive
	 * another number, addresses(global), its local adress.
	 * return the adress where the vertex lay on L.
	 */
	Vertex_handle add_vertex(Root_vertex_handle global);

private:
	/**
	 * Constructs a subcomplex which is empty
	 */
	Skeleton_blocker_sub_complex(){};


public:
	void clear();

	/**
	 * Compute the local vertex in L corresponding to the vertex global in K.
	 * runs in O(log n) if n = num_vertices()
	 */
	boost::optional<Vertex_handle> get_address(Root_vertex_handle global) const;

	//	/**
	//	 * Allocates a simplex in L corresponding to the simplex s in K
	//	 * with its local adresses and returns an AddressSimplex.
	//	 */
	//	boost::optional<Simplex_handle> get_address(const Root_simplex_handle & s) const;

private:
	/**
	 *  same as get_address except that it will return a simplex in any case.
	 *  The vertices that were not found are not added.
	 */
	std::vector<boost::optional<Vertex_handle> > get_addresses(const Root_simplex_handle & s) const{
		std::vector<boost::optional<Vertex_handle> > res;
		for (auto i : s)
		{
			res.push_back(get_address(i));
		}
		return res;
	}

};


template<typename ComplexType>
typename ComplexType::Vertex_handle
Skeleton_blocker_sub_complex<ComplexType>::add_vertex(Root_vertex_handle id)
{
	Vertex_handle address(boost::add_vertex(this->skeleton));
	(*this)[address].activate();
	(*this)[address].set_id(id);
	adresses.insert(AddressPair(id,address));
	this->num_vertices_++;
	this->degree.push_back(0);
	return address;
}


template<typename ComplexType>
boost::optional<typename ComplexType::Vertex_handle>
Skeleton_blocker_sub_complex<ComplexType>::get_address(typename ComplexType::Root_vertex_handle i) const
{
	boost::optional<Vertex_handle> res;
	IdAddressMapConstIterator it = adresses.find(i);
	if (it == adresses.end()) res.reset();
	else  res=(*it).second;
	return res;
}

template<typename ComplexType>
void
Skeleton_blocker_sub_complex<ComplexType>::clear(){
	adresses.clear();
	ComplexType::clear();
}


/**
 * @remark remarque perte de temps a creer un nouveau simplexe a chaque fois
 * alors qu'on pourrait utiliser a la place de 'addresses_sigma_in_link'
 * un simplex avec des valeurs spéciales ComplexDS::null_vertex par exemple
 * pour indiquer qu'un vertex n'appartient pas au complex
 */
template<typename ComplexType>
bool proper_face_in_union(
		Skeleton_blocker_sub_complex<ComplexType> & link,
		std::vector<boost::optional<typename ComplexType::Vertex_handle> > & addresses_sigma_in_link,
		int vertex_to_be_ignored)
{
	// we test that all vertices of 'addresses_sigma_in_link' but 'vertex_to_be_ignored'
	// are in link1 if it is the case we construct the corresponding simplex
	bool vertices_sigma_are_in_link = true;
	typename ComplexType::Simplex_handle sigma_in_link;
	for(int i = 0; i< addresses_sigma_in_link.size(); ++i){
		if (i!= vertex_to_be_ignored){
			if(!addresses_sigma_in_link[i]){
				vertices_sigma_are_in_link = false;
				break;
			}
			else sigma_in_link.add_vertex(*addresses_sigma_in_link[i]);
		}
	}
	// If one of vertices of the simplex is not in the complex then it returns false
	// Otherwise, it tests if the simplex is in the complex
	return vertices_sigma_are_in_link && link.contains(sigma_in_link);
}


template<typename ComplexType>
bool
proper_faces_in_union(Skeleton_blocker_simplex<typename ComplexType::Root_vertex_handle> & sigma, Skeleton_blocker_sub_complex<ComplexType> & link1, Skeleton_blocker_sub_complex<ComplexType> & link2)
{
	typedef typename ComplexType::Vertex_handle  Vertex_handle;
	std::vector<boost::optional<Vertex_handle> > addresses_sigma_in_link1 = link1.get_addresses(sigma);
	std::vector<boost::optional<Vertex_handle> > addresses_sigma_in_link2 = link2.get_addresses(sigma);

	for (int current_index = 0; current_index < addresses_sigma_in_link1.size() ; ++current_index)
	{

		if (!proper_face_in_union(link1,addresses_sigma_in_link1,current_index)
				&&	!proper_face_in_union(link2,addresses_sigma_in_link2,current_index)){
			return false;
		}
	}
	return true;
}





#endif

