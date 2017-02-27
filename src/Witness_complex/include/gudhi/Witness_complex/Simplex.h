#ifndef SIMPLEX_H_
#define SIMPLEX_H_

#include <gudhi/Active_witness/Sib_vertex_pair.h>
#include <vector>

namespace Gudhi {
namespace witness_complex {  
  
template <typename SimplicialComplexForWitness>
class Simplex {

private:
  typedef typename SimplicialComplexForWitness::Vertex_handle Vertex_handle;
  typedef Sib_vertex_pair<SimplicialComplexForWitness, Vertex_handle> Simplex_key;
  typedef std::vector<Vertex_handle> Vertex_vector;
  typedef typename SimplicialComplexForWitness::Simplex_vertex_range Simplex_vertex_range;
  
  bool is_simplex_key_;
  const SimplicialComplexForWitness& complex_;
  Simplex_key* simplex_key_;
  Vertex_vector* vertex_vector_;

public:  
  /* Simplex key constructor */
  Simplex(const Simplex_key& input, const SimplicialComplexForWitness& complex)
    : is_simplex_key_(true), complex_(complex)
  {
    simplex_key_ = new Simplex_key(input);
  }

  /* Vertex vector constructor */
  Simplex(const Vertex_vector& input, const SimplicialComplexForWitness& complex)
    : is_simplex_key_(false), complex_(complex)
  {
    vertex_vector_ = new Vertex_vector(input);
  }

  /* Copy constructor */
  Simplex(const Simplex& s)
    : is_simplex_key_(s.is_simplex_key_)
  {
    if (is_simplex_key_)
      simplex_key_ = new Simplex_key(*s.simplex_key_);
    else
      vertex_vector_ = new Vertex_vector(*s.vertex_vector_);
  }

  /* Move constructor */
  Simplex(Simplex&& s) noexcept
    : is_simplex_key_(s.is_simplex_key_)
  {
    if (is_simplex_key_) {
      simplex_key_ = new Simplex_key(*s.simplex_key_);
      s.simplex_key_ = nullptr;
    }
    else {
      vertex_vector_ = new Vertex_vector(*s.vertex_vector_);
      s.vertex_vector_ = nullptr;
    }
  }
  
  /* Copy assignment */
  Simplex& operator= (const Simplex& s)
  {
    Simplex tmp(s);
    *this = std::move(tmp);
    return *this;
  }

  /** Move assignment operator */
  Simplex& operator= (Simplex&& s) noexcept
  {
    if (is_simplex_key_) {
      delete simplex_key_;
      simplex_key_ = s.simplex_key_;
      s.simplex_key_ = nullptr;
      return *this;
    }
    else {
      delete vertex_vector_;
      vertex_vector_ = new Vertex_vector(*s.vertex_vector_);
      s.vertex_vector_ = nullptr;
      return *this;
    }
  }  

  /* Destructor */
  ~Simplex() noexcept
  {
    if (is_simplex_key_)
      delete simplex_key_;
    else
      delete vertex_vector_;
  }

  Simplex_vertex_range vertex_range()
  {
    if (is_simplex_key_)
      return complex_.simplex_vertex_range(simplex_key_->simplex_handle());
    else
      return *vertex_vector_;
  }

  Simlex_handle simplex_handle()
  {
    if (is_simplex_key_)
      return simplex_key_->simplex_handle();
    else
      return complex_.find(*vertex_vector_);
  }
};

}
}
#endif
