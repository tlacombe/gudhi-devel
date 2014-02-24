#ifndef __GIGA_SIMPLEX_H
#define __GIGA_SIMPLEX_H

#include<cassert>
#include<iostream>
#include<set>
#include<vector>

using namespace std;

/**
 *@class Simplex
 *@brief Abstract simplex
 *
 * An abstract simplex is represented as an ordered set of T elements,
 * each element representing a vertex.
 * 
 * The element representing a vertex can be a 'ComplexDS::Vertex_handle' or a  'ComplexDS::Root_vertex_handle'
 */
template <typename T>

class Simplex {

private :
	set<T> simplex_set;

public:
	typedef typename T::boost_vertex_handle boost_vertex_handle;


	class Simplex_vertex_const_iterator{
	private:
		typename std::set<T>::iterator set_iterator;
		typedef typename std::set<T>::const_iterator const_set_iterator;
		friend class Simplex;

		//		Simplex_vertex_const_iterator(const_set_iterator& other):set_iterator(other){}
	public:

		explicit Simplex_vertex_const_iterator(typename std::set<T>::const_iterator other):set_iterator(other){}


		explicit Simplex_vertex_const_iterator() :	set_iterator(){}

		Simplex_vertex_const_iterator& operator=(const Simplex_vertex_const_iterator& other)
		{
			set_iterator = other.set_iterator;
			return(*this);
		}

		bool operator!= (const Simplex_vertex_const_iterator & other) const{
			if(set_iterator == other.set_iterator) return false;
			else return true;
		}

		T operator* () const{
			return *set_iterator;
		}

		const Simplex_vertex_const_iterator & operator++ (){
			set_iterator++;
			return *this;
		}

		Simplex_vertex_const_iterator  operator++ (int){
			Simplex_vertex_const_iterator tmp(*this);
			tmp.set_iterator++;
			return(tmp);
		}
	};


	typedef Simplex_vertex_const_iterator const_iterator;

	/** @name Constructors / Destructors / Initialization
	 */
	//@{

	/**
	 * Constructs the empty simplex {}
	 */
	Simplex():simplex_set() {}

	/**
	 * Clear the simplex
	 */
	inline void clear() {
		simplex_set.clear();
	}

	/**
	 * Constructs the singleton {a}
	 */
	Simplex(T a)
	{
		add_vertex(a);
	}

	/**
	 * Constructs the edge {a,b}
	 */
	Simplex(T a, T b)
	{
		add_vertex(a); add_vertex(b);
	}

	/**
	 * Constructs the triangle {a,b,c}
	 */
	Simplex(T a, T b, T c)
	{
		add_vertex(a); add_vertex(b); add_vertex(c);
	}

	/**
	 * Constructs the tetrahedron {a,b,c,d}
	 */
	Simplex(T a, T b, T c, T d)
	{
		add_vertex(a); add_vertex(b); add_vertex(c); add_vertex(d);
	}

	/**
	 * Initialize a simplex with a string such as {0,1,2}
	 */
	Simplex(string token){
		clear();
		if ((token[0] == '{')  && (token[token.size()-1] == '}' ) )
		{
			token.erase(0,1);
			token.erase(token.size()-1,1);
			while(token.size()!=0 ){
				int coma_position=token.find_first_of(',');
				//cout << "coma_position:"<<coma_position<<endl;
				string n = token.substr(0,coma_position);
				//cout << "token:"<<token<<endl;
				token.erase(0,n.size()+1);
				add_vertex((T)(atoi(n.c_str())));
			}
		}
	}


	/**
	 * Constructs the Simplex which is the same than sigma
	 */
	Simplex(const Simplex* sigma)
	{
		for (const_iterator i=sigma->begin(); i!=sigma->end(); ++i)
			simplex_set.insert(*i);
	}

	//@}

	/** @name Simplex manipulation
	 */
	//@{

	/**
	 * Add the vertex v to the simplex:
	 * \f$ (*this) \leftarrow (*this) \cup \{ v \} \f$
	 */
	void inline add_vertex(T v)
	{
		simplex_set.insert(v);
	}

	/**
	 * Remove the vertex v from the simplex:
	 * \f$ (*this) \leftarrow (*this) \setminus \{ v \} \f$
	 */
	void inline remove_vertex(T v)
	{
		simplex_set.erase(v);
	}

	/**
	 * Intersects the simplex with the simplex a:
	 * \f$ (*this) \leftarrow (*this) \cap a \f$
	 */
	void intersection(const Simplex & a){
		vector<T> v;
		v.reserve(std::min(simplex_set.size(), a.simplex_set.size()));

		set_intersection(simplex_set.begin(),simplex_set.end(),
				a.simplex_set.begin(),a.simplex_set.end(),
				std::back_inserter(v));
		clear();
		for (auto i:v)
			simplex_set.insert(i);
	}

	/**
	 * Substracts a from the simplex:
	 * \f$ (*this) \leftarrow (*this) \setminus a \f$
	 */
	void difference(const Simplex & a){
		vector<T> v;
		v.reserve(simplex_set.size());

		set_difference(simplex_set.begin(),simplex_set.end(),
				a.simplex_set.begin(),a.simplex_set.end(),
				std::back_inserter(v));

		clear();
		for (auto i:v)
			simplex_set.insert(i);
	}

	/**
	 * Add vertices of a to the simplex:
	 * \f$ (*this) \leftarrow (*this) \cup a \f$
	 */
	void union_vertices(const Simplex & a){
		vector<T> v;
		v.reserve(simplex_set.size() + a.simplex_set.size());

		set_union(simplex_set.begin(),simplex_set.end(),
				a.simplex_set.begin(),a.simplex_set.end(),
				std::back_inserter(v));

		clear();
		for (typename vector<T>::iterator i=v.begin(); i!=v.end(); ++i)
			simplex_set.insert(*i);
	}



	const_iterator begin() const{
		return const_iterator(simplex_set.cbegin());
	}

	const_iterator end() const{
		return const_iterator(simplex_set.cend());
	}


	//@}

	/** @name Queries
	 */
	//@{

	/**
	 * Returns the dimension of the simplex.
	 */
	int inline dimension() const
	{
		return (simplex_set.size() - 1);
	}

	/**
	 * Returns the first vertex of the (oriented) simplex.
	 *
	 * Be careful : assumes the simplex is non-empty.
	 */
	T inline first_vertex() const
	{
		return *(begin());
	}

	/**
	 * Returns the second vertex of the (oriented) simplex.
	 *
	 * Be careful : assumes the simplex has at least two vertices.
	 */
	int inline second_vertex() const
	{
		const_iterator it=(begin());
		++it;
		return *it;
	}

	/**
	 * Returns the last vertex of the (oriented) simplex.
	 *
	 * Be careful : assumes the simplex is non-empty.
	 */
	T inline last_vertex() const
	{
		assert(!simplex_set.empty());
		return *(simplex_set.rbegin());
	}
	/**
	 * @return true iff the simplex contains the simplex a, i.e. iff \f$ a \subset (*this) \f$.
	 */
	bool contains(const Simplex & a) const{
		return includes(simplex_set.cbegin(),simplex_set.cend(),a.simplex_set.cbegin(),a.simplex_set.cend());
	}

	/**
	 * @return true iff the simplex contains the vertex v, i.e. iff \f$ v \in (*this) \f$.
	 */
	bool contains(T v) const{
		return (simplex_set.find(v) != simplex_set.end());

	}

	/**
	 * @return \f$ (*this) \cap a = \emptyset \f$.
	 */
	bool disjoint(const Simplex& a) const{
		vector<T> v;
		v.reserve(std::min(simplex_set.size(), a.simplex_set.size()));

		set_intersection(simplex_set.cbegin(),simplex_set.cend(),
				a.simplex_set.cbegin(),a.simplex_set.cend(),
				std::back_inserter(v));

		return (v.size()==0);
	}


	bool operator==(const Simplex& other) const{
		return (this->simplex_set == other.simplex_set);
	}

	bool operator!=(const Simplex& other) const{
		return (this->simplex_set != other.simplex_set);
	}

	//@}



	/**
	 * Display a simplex
	 */
	friend ostream& operator << (ostream& o, const Simplex & sigma)
	{
		int first = 1;
		const_iterator i = sigma.begin();
		o << '{';
		while(i != sigma.end())
		{
			if (first) first = 0; else o << ',';
			o << *i;
			++i;
		}
		o << '}';
		return o;
	}


};




#endif


