#ifndef SRC_SKELETON_BLOCKER_INCLUDE_GUDHI_MAP_BLOCKERS_H_
#define SRC_SKELETON_BLOCKER_INCLUDE_GUDHI_MAP_BLOCKERS_H_

#include <map>
#include <boost/iterator/transform_iterator.hpp>
#include <boost/range/adaptor/map.hpp>
#include "boost/iterator/iterator_facade.hpp"

namespace Gudhi {

namespace skbl {

/**
 * @brief Iterator through the blockers of a vertex.
  */
// ReturnType = const Simplex* or Simplex*
// MapIteratorType = BlockerMapConstIterator or BlockerMapIterator
template<typename MapIteratorType, typename ReturnType>
class Blocker_map_iterator_internal : public boost::iterator_facade<
  Blocker_map_iterator_internal<MapIteratorType,ReturnType>,
  ReturnType,
  boost::forward_traversal_tag,
  ReturnType
  >{
private:
  MapIteratorType current_position;
  MapIteratorType end_of_map;
public:

  Blocker_map_iterator_internal():current_position(){}

  Blocker_map_iterator_internal(MapIteratorType position,MapIteratorType end_of_map_ ):
    current_position(position), end_of_map(end_of_map_)
  { }

  bool equal(const Blocker_map_iterator_internal& other) const{
    return current_position == other.current_position;
  }

  void increment(){ goto_next_blocker(); }

  ReturnType dereference() const  { return(current_position->second); }

private:
  /**
   * Let the current pair be (v,sigma) where v is a vertex and sigma is a blocker.
   * If v is not the first vertex of sigma then we already have seen sigma as a blocker
   * and we look for the next one.
   */
  void goto_next_blocker(){
    do {
      ++current_position;
    } while (!(current_position == end_of_map) && !first_time_blocker_is_seen());
  }

  bool first_time_blocker_is_seen() const{
    return current_position->first  == current_position->second->first_vertex();
  }
};



/**
 * @brief Iterator through the blockers of a vertex
 */
// ReturnType = const Simplex* or Simplex*
// MapIteratorType = BlockerMapConstIterator or BlockerMapIterator
template<typename MapIteratorType, typename ReturnType>
class Blocker_map_iterator_around_vertex_internal : public boost::iterator_facade<
  Blocker_map_iterator_around_vertex_internal<MapIteratorType,ReturnType>,
  ReturnType,
  boost::forward_traversal_tag,
  ReturnType
>{
private:
  MapIteratorType current_position_;
public:

  Blocker_map_iterator_around_vertex_internal():current_position_(){}

  Blocker_map_iterator_around_vertex_internal(MapIteratorType position):
    current_position_(position)
  {}

  Blocker_map_iterator_around_vertex_internal& operator=(Blocker_map_iterator_around_vertex_internal other){
    this->current_position_ = other.current_position_;
    return *this;
  }

  bool equal(const Blocker_map_iterator_around_vertex_internal& other) const{
    return current_position_ == other.current_position_;
  }

  void increment(){ current_position_++; }

  ReturnType dereference() const{ return(current_position_->second); }
  MapIteratorType current_position(){ return this->current_position_; }
};

template<typename Vertex_handle> class Blockers_map_ds {

public:
  typedef Skeleton_blocker_simplex<Vertex_handle> Simplex;
  typedef Simplex* Blocker_handle;
  typedef std::multimap<Vertex_handle, Simplex *> BlockerMap;
  typedef typename std::multimap<Vertex_handle, Simplex *>::value_type BlockerPair;
  typedef typename std::multimap<Vertex_handle, Simplex *>::iterator BlockerMapIterator;
  typedef typename std::multimap<Vertex_handle, Simplex *>::const_iterator BlockerMapConstIterator;

private:
  BlockerMap blocker_map_;
  int num_blockers_ = 0;

public:
  template <typename Complex>
  Blockers_map_ds(Complex& complex){
  }


  ~Blockers_map_ds() {
    clear();
  }
  
  void clear() {
    while (!blocker_map_.empty()) {
      delete_blocker(blocker_map_.begin()->second);
    }
    blocker_map_.clear();
  }

  int num_blockers() const {
    return num_blockers_;
  }

  Blocker_handle add_blocker(const Simplex& blocker) {
      Blocker_handle blocker_pt = new Simplex(blocker);
      for(auto v : *blocker_pt)
        blocker_map_.insert(BlockerPair(v, blocker_pt));
      ++num_blockers_;
      return blocker_pt;
  }
 
  void remove_blocker(const Blocker_handle sigma) {
    --num_blockers_;
  	for (auto vertex : *sigma)
        remove_blocker(sigma, vertex);
  }
 
private:
  /**
   * Removes sigma from the blocker map of vertex v
   */
  void remove_blocker(const Blocker_handle sigma, Vertex_handle v) {
    Blocker_around_vertex_iterator blocker;
    for (blocker = blocker_range(v).begin(); blocker != blocker_range(v).end();
         ++blocker) {
      if (*blocker == sigma)
        break;
    }
    if (*blocker != sigma) {
      std::cerr
          << "bug ((*blocker).second == sigma) ie try to remove a blocker not present\n";
      assert(false);
    } else {
      blocker_map_.erase(blocker.current_position());
    }
  }

public:
  void delete_blocker(Blocker_handle sigma) {
    remove_blocker(sigma);
    delete sigma;
  }

  bool is_included_in_one_blocker(Vertex_handle v) const {
    return blocker_map_.lower_bound(v) != blocker_map_.upper_bound(v);
  }

  /**
   * @return true iff s is a blocker of the simplicial complex
   */
  bool contains_blocker(const Simplex & s) const {
    if (s.dimension() < 2)
      return false;

    Vertex_handle a = s.first_vertex();

    for (auto blocker : const_blocker_range(a)) {
      if (s == *blocker)
        return true;
    }
    return false;
  }


  typedef Blocker_map_iterator_around_vertex_internal<
  typename std::multimap<Vertex_handle, Simplex *>::iterator,
  Blocker_handle>
  Blocker_around_vertex_iterator;

  /**
   * @brief Iterator over (constant) blockers adjacent to a vertex
   */
  typedef Blocker_map_iterator_around_vertex_internal<
  typename std::multimap<Vertex_handle, Simplex *>::const_iterator,
  const Blocker_handle>
  Const_blocker_around_vertex_iterator;

  typedef boost::iterator_range <Blocker_around_vertex_iterator> Blocker_around_vertex_range;
  typedef boost::iterator_range <Const_blocker_around_vertex_iterator> Const_blocker_around_vertex_range;

 public:
  /**
   * @brief Returns a range of the blockers of the complex passing through a vertex
   */
  Blocker_around_vertex_range blocker_range(Vertex_handle v) {
    auto begin = Blocker_around_vertex_iterator(blocker_map_.lower_bound(v));
    auto end = Blocker_around_vertex_iterator(blocker_map_.upper_bound(v));
    return Blocker_around_vertex_range(begin, end);
  }

  /**
   * @brief Returns a range of the blockers of the complex passing through a vertex
   */
  Const_blocker_around_vertex_range const_blocker_range(Vertex_handle v) const {
    auto begin = Const_blocker_around_vertex_iterator(blocker_map_.lower_bound(v));
    auto end = Const_blocker_around_vertex_iterator(blocker_map_.upper_bound(v));
    return Const_blocker_around_vertex_range(begin, end);
  }


  /**
   * @brief Iterator over the blockers.
   */
  typedef Blocker_map_iterator_internal<
  typename std::multimap<Vertex_handle, Simplex *>::iterator,
  Blocker_handle>
  Blocker_iterator;

  /**
   * @brief Iterator over the (constant) blockers.
   */
  typedef Blocker_map_iterator_internal<
  typename std::multimap<Vertex_handle, Simplex *>::const_iterator,
  const Blocker_handle>
  Const_blocker_iterator;

  typedef boost::iterator_range <Blocker_iterator> Blocker_range;
  typedef boost::iterator_range <Const_blocker_iterator> Const_blocker_range;

 public:
  /**
   * @brief Returns a range of the blockers of the complex
   */
  Blocker_range blocker_range() {
    auto begin = Blocker_iterator(blocker_map_.begin(), blocker_map_.end());
    auto end = Blocker_iterator(blocker_map_.end(), blocker_map_.end());
    return Blocker_range(begin, end);
  }

    /**
   * @brief Returns a range of the blockers of the complex
   */
  Const_blocker_range const_blocker_range() const{
    auto begin = Const_blocker_iterator(blocker_map_.begin(), blocker_map_.end());
    auto end = Const_blocker_iterator(blocker_map_.end(), blocker_map_.end());
    return Const_blocker_range(begin, end);
  }



};

}

}


#endif 