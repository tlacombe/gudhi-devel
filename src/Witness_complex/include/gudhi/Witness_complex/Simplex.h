#ifndef SIMPLEX_H_
#define SIMPLEX_H_

#include <gudhi/Active_witness/Sib_vertex_pair.h>
#include <vector>

namespace Gudhi {
namespace witness_complex {  
  
  template <typename SimplicialComplexForWitness>
class Simplex {

private:
    // typedef typename SimplicialComplexForWitness::Vertex_handle Vertex_handle;
  typedef typename SimplicialComplexForWitness::Simplex_handle Simplex_handle;
  typedef typename SimplicialComplexForWitness::Siblings Siblings;
  typedef Sib_vertex_pair<SimplicialComplexForWitness, Vertex_handle> Simplex_key;
    typedef std::vector<std::size_t> Vertex_vector;
  typedef typename SimplicialComplexForWitness::Simplex_vertex_range Simplex_vertex_range;
  
  bool is_simplex_key_;
  SimplicialComplexForWitness& complex_;
  Simplex_key simplex_key_;
  Vertex_vector vertex_vector_;

public:  

  /* Simplex handle constructor */
  Simplex(Simplex_handle sh, SimplicialComplexForWitness& complex)
    : is_simplex_key_(true), complex_(complex)
  {
    Siblings* sib = complex.self_siblings(sh);
    Vertex_handle v = sh->first;
    simplex_key_ = Simplex_key(sib, v);
  }
  
  /* Vertex vector constructor */
  Simplex(Vertex_vector& input, SimplicialComplexForWitness& complex)
    : is_simplex_key_(false), complex_(complex)
  {
    vertex_vector_ = Vertex_vector(input);
  }


    // /* Simplex key constructor */
  // Simplex(Simplex_key& input,  SimplicialComplexForWitness& complex)
  //   : is_simplex_key_(true), complex_(complex)
  // {
  //   simplex_key_ = new Simplex_key(input);
  // }

  // /* Simplex handle constructor */
  // Simplex(Simplex_handle sh, SimplicialComplexForWitness& complex)
  //   : is_simplex_key_(true), complex_(complex)
  // {
  //   Siblings* sib = complex.self_siblings(sh);
  //   Vertex_handle v = sh->first;
  //   simplex_key_ = new Simplex_key(sib, v);
  // }
  
  // /* Vertex vector constructor */
  // Simplex(Vertex_vector& input, SimplicialComplexForWitness& complex)
  //   : is_simplex_key_(false), complex_(complex)
  // {
  //   vertex_vector_ = new Vertex_vector(input);
  // }

  // /* Copy constructor */
  // Simplex(const Simplex& s)
  //   : is_simplex_key_(s.is_simplex_key_), complex_(s.complex_)
  // {
  //   if (is_simplex_key_)
  //     simplex_key_ = new Simplex_key(*s.simplex_key_);
  //   else
  //     vertex_vector_ = new Vertex_vector(*s.vertex_vector_);
  // }

  // /* Move constructor */
  // Simplex(Simplex&& s) noexcept
  //   : is_simplex_key_(s.is_simplex_key_), complex_(s.complex_)
  // {
  //   if (is_simplex_key_) {
  //     simplex_key_ = new Simplex_key(*s.simplex_key_);
  //     s.simplex_key_ = nullptr;
  //   }
  //   else {
  //     vertex_vector_ = new Vertex_vector(*s.vertex_vector_);
  //     s.vertex_vector_ = nullptr;
  //   }
  // }
  
  // /* Copy assignment */
  // Simplex& operator= (const Simplex& s)
  // {
  //   Simplex tmp(s);
  //   *this = std::move(tmp);
  //   return *this;
  // }

  // /** Move assignment operator */
  // Simplex& operator= (Simplex&& s) noexcept
  // {
  //   if (is_simplex_key_) {
  //     delete simplex_key_;
  //     simplex_key_ = s.simplex_key_;
  //     complex_ = s.complex_;
  //     s.simplex_key_ = nullptr;
  //     return *this;
  //   }
  //   else {
  //     delete vertex_vector_;
  //     vertex_vector_ = new Vertex_vector(*s.vertex_vector_);
  //     complex_ = s.complex_;
  //     s.vertex_vector_ = nullptr;
  //     return *this;
  //   }
  // }  

  // /* Destructor */
  // ~Simplex() noexcept
  // {
  //   if (is_simplex_key_)
  //     delete simplex_key_;
  //   else
  //     delete vertex_vector_;
  // }

  // // Can't put const because of simplex_vertex_range :(
  // Vertex_vector vertex_range()
  // {
  //   if (is_simplex_key_) {
  //     Simplex_vertex_range sv_range = complex_.simplex_vertex_range(simplex_key_->simplex_handle());
  //     return Vertex_vector(sv_range.begin(), sv_range.end());
  //   }
  //   else
  //     return Vertex_vector(*vertex_vector_);
  // }

  // // Can't put const because of find :(
  // Simplex_handle simplex_handle()
  // {
  //   if (is_simplex_key_)
  //     return simplex_key_->simplex_handle();
  //   else
  //     return complex_.find(*vertex_vector_);
  // }

  // Can't put const because of simplex_vertex_range :(
  Vertex_vector vertex_range()
  {
    if (is_simplex_key_) {
      Simplex_vertex_range sv_range = complex_.simplex_vertex_range(simplex_key_.simplex_handle());
      return Vertex_vector(sv_range.begin(), sv_range.end());
    }
    else
      return Vertex_vector(vertex_vector_);
  }

  // Can't put const because of find :(
  Simplex_handle simplex_handle()
  {
    if (is_simplex_key_)
      return simplex_key_.simplex_handle();
    else
      return complex_.find(vertex_vector_);
  }
    
    
  /* A comparison operator */
  bool operator <(const Simplex& rhs) const
  {
    return (&simplex_key_ < &rhs.simplex_key_) || (&vertex_vector_ < &rhs.vertex_vector_);
  }

  bool operator ==(const Simplex& rhs) const
  {
    return (simplex_key_.first == rhs.simplex_key_.first &&
            simplex_key_.second == rhs.simplex_key_.second) ||
           (vertex_vector_ == rhs.vertex_vector_);
  }
};

}
}
#endif
