#ifndef SIMPLEX_TREE_H
#define SIMPLEX_TREE_H

#include <vector>
#include <unordered_map>
#include <queue>
#include <stack>
#include <algorithm>

namespace Gudhi {
namespace tmp_package_name {

class Simplex_tree
{
public:
    class Simplicial_tree_node;
    using vertex = double;
    using index = double;
    using simplex_base = std::vector<vertex>;
    using labels_dictionary = std::unordered_map<vertex, Simplicial_tree_node*>;

    Simplex_tree();
    ~Simplex_tree();

    class Simplicial_tree_node
    {
    public:
        Simplicial_tree_node(vertex label, index insertionIndex, int dim, Simplicial_tree_node *parent);
        ~Simplicial_tree_node();

        vertex get_label() const;
        index get_insertion_index() const;

        Simplicial_tree_node* get_parent() const;
        void set_parent(Simplicial_tree_node *value);
        labels_dictionary *get_children() const;

        Simplicial_tree_node* get_next() const;
        void set_next(Simplicial_tree_node *value);
        Simplicial_tree_node* get_prev() const;
        void set_prev(Simplicial_tree_node *value);

        int get_dim() const;

    private:
        vertex label_;
        index insertionIndex_;
        int dim_;
        Simplicial_tree_node *parent_;
        labels_dictionary *children_;
        Simplicial_tree_node *next_;
        Simplicial_tree_node *prev_;
    };

    /* Important for tower_converter.h */

    Simplicial_tree_node* insert_simplex(simplex_base *numVertices);
    void remove_simplex(simplex_base *simplex);
    void remove_simplex(simplex_base *simplex, std::vector<index> *removedIndices);
    vertex get_smallest_closed_star(vertex v, vertex u, std::vector<simplex_base*> *closedStar);   //returns vertex with smallest closed star; closedStar is ordered
    index get_boundary(simplex_base *simplex, std::vector<index> *boundary);    //ordered by increasing insertion numbers
    double get_size() const;
    double get_max_index() const;

    /* Other */

    //index contract_vertices(vertex v, vertex u, double timestamp, std::vector<std::vector<index>*> *boundaries, std::vector<index> *&inactiveInsertionNumbers);
    double get_max_size() const;
    int get_max_dimension() const;
    void get_cofaces(simplex_base *simplex, std::vector<simplex_base*> *cofaces);
    void get_cofacets(simplex_base *simplex, std::vector<simplex_base *> *cofacets);
    bool contains(simplex_base *simplex);

private:
    std::vector<labels_dictionary*> *dictionaries_; // circular list for each label at each height
    std::unordered_map<double, vertex> *verticesTranslation_;
    std::vector<double> *numberOfSimplicesPerDim_; //number of simplices in each dimension
    double numberOfSimplices_;
    double maxIndex_;
    double maxSize_;
    int maxDim_;

    Simplicial_tree_node* insert_simplex_in_tree(simplex_base *simplex);
    void insert_node_in_dictionary(Simplicial_tree_node *simplex);
    void delete_node_from_dictionary(Simplicial_tree_node *simplex);
    void remove_simplex(Simplicial_tree_node *simplex, std::vector<index> *removedIndices);
    Simplicial_tree_node* find(simplex_base *simplex);
    void get_cofaces(Simplicial_tree_node *simplex, std::vector<Simplicial_tree_node*> *cofaces);
    void get_cofacets(Simplicial_tree_node *simplex, std::vector<Simplicial_tree_node*> *cofacets);
    bool is_coface(Simplicial_tree_node *node, Simplicial_tree_node *simplex);
    void get_simplices_in_subtree(Simplicial_tree_node *node, std::vector<Simplicial_tree_node*> *simplices);
    vertex get_smallest_star(vertex v, vertex u, std::queue<Simplicial_tree_node*> *qv, std::queue<Simplicial_tree_node*> *qu);
    bool get_smallest_star_find_next_node(vertex v, vertex toAvoid, std::vector<labels_dictionary*>::size_type &currentHeight,
                                          Simplicial_tree_node *&currentRoot, Simplicial_tree_node *&startRoot, Simplicial_tree_node *&currentNode, bool &currentNodeIsRoot,
                                          labels_dictionary *&children, labels_dictionary::iterator &it, std::queue<Simplicial_tree_node*> *&tail);
    bool get_smallest_star_find_next_root(vertex v, vertex toAvoid, std::vector<labels_dictionary*>::size_type &currentHeight,
                                          Simplicial_tree_node *&currentRoot, Simplicial_tree_node *&startRoot, Simplicial_tree_node *&currentNode, bool &currentNodeIsRoot);
    simplex_base* node_to_vector(Simplicial_tree_node *node);
    Simplicial_tree_node* get_opposite_facet(Simplicial_tree_node *simplex, vertex v);
};

Simplex_tree::Simplex_tree() : numberOfSimplices_(0), maxIndex_(-1), maxSize_(0), maxDim_(0)
{
    verticesTranslation_ = new std::unordered_map<double, vertex>();
    dictionaries_ = new std::vector<labels_dictionary*>();
    numberOfSimplicesPerDim_ = new std::vector<double>();
}

Simplex_tree::~Simplex_tree()
{
    delete verticesTranslation_;
    for (std::vector<labels_dictionary*>::size_type i = 0; i < dictionaries_->size(); i++){
        for (labels_dictionary::iterator it = dictionaries_->at(i)->begin(); it != dictionaries_->at(i)->end(); it++) delete it->second;
        delete dictionaries_->at(i);
    }
    delete dictionaries_;
    delete numberOfSimplicesPerDim_;
}

Simplex_tree::Simplicial_tree_node* Simplex_tree::insert_simplex(simplex_base *numVertices)
{
    Simplicial_tree_node *simplexNode = insert_simplex_in_tree(numVertices);
    if (simplexNode == NULL) return simplexNode;

    int dim = numVertices->size() - 1;

    if ((int)numberOfSimplicesPerDim_->size() == dim) numberOfSimplicesPerDim_->push_back(1);
    else numberOfSimplicesPerDim_->at(dim)++;
    numberOfSimplices_++;
    maxIndex_++;
    if (maxSize_ < numberOfSimplices_) maxSize_ = numberOfSimplices_;
    if (maxDim_ < dim) maxDim_ = dim;

    return simplexNode;
}

inline void Simplex_tree::remove_simplex(simplex_base *simplex)
{
    Simplicial_tree_node *simplexNode = find(simplex);
    remove_simplex(simplexNode, NULL);
}

inline void Simplex_tree::remove_simplex(simplex_base *simplex, std::vector<index> *removedIndices)
{
    Simplicial_tree_node *simplexNode = find(simplex);
    remove_simplex(simplexNode, removedIndices);
}

Simplex_tree::vertex Simplex_tree::get_smallest_closed_star(vertex v, vertex u, std::vector<simplex_base *> *closedStar)
{
    std::queue<Simplicial_tree_node*> qv;
    std::queue<Simplicial_tree_node*> qu;
    Simplicial_tree_node *s;

    if (get_smallest_star(v, u, &qv, &qu) == v){
        while (!qv.empty()){
            s = qv.front();
            qv.pop();
            if (s->get_parent() != NULL){
                closedStar->push_back(node_to_vector(get_opposite_facet(s, v)));
            }
            closedStar->push_back(node_to_vector(s));
        }

        return v;
    } else {
        while (!qu.empty()){
            s = qu.front();
            qu.pop();
            if (s->get_parent() != NULL){
                closedStar->push_back(node_to_vector(get_opposite_facet(s, u)));
            }
            closedStar->push_back(node_to_vector(s));
        }

        return u;
    }
}

Simplex_tree::index Simplex_tree::get_boundary(simplex_base *simplex, std::vector<index> *boundary)
{
    Simplicial_tree_node *simplexNode = find(simplex);
    if (simplexNode->get_parent() == NULL) return simplexNode->get_insertion_index();

    simplex_base tail;
    tail.push_back(simplexNode->get_label());
    Simplicial_tree_node *it = simplexNode->get_parent();
    simplex_base::reverse_iterator tailIt;
    Simplicial_tree_node *itb;

    while (it != NULL){
        itb = it;
        tailIt = tail.rbegin();
        tailIt++;
        while (tailIt != tail.rend()){
            itb = itb->get_children()->at(*tailIt);
            tailIt++;
        }
        boundary->push_back(itb->get_insertion_index());
        tail.push_back(it->get_label());
        it = it->get_parent();
    }
    itb = dictionaries_->front()->at(tail.at(tail.size() - 2));
    tailIt = tail.rbegin();
    tailIt += 2;
    while (tailIt != tail.rend()){
        itb = itb->get_children()->at(*tailIt);
        tailIt++;
    }
    boundary->push_back(itb->get_insertion_index());

    std::sort(boundary->begin(), boundary->end());
    return simplexNode->get_insertion_index();
}

inline double Simplex_tree::get_size() const
{
    return numberOfSimplices_;
}

inline double Simplex_tree::get_max_index() const
{
    return maxIndex_;
}

inline double Simplex_tree::get_max_size() const
{
    return maxSize_;
}

inline int Simplex_tree::get_max_dimension() const
{
    return maxDim_;
}

void Simplex_tree::get_cofaces(simplex_base *simplex, std::vector<simplex_base *> *cofaces)
{
    Simplicial_tree_node *simplexNode = find(simplex);
    std::vector<Simplicial_tree_node*> cofaceNodes;
    get_cofaces(simplexNode, &cofaceNodes);

    for (std::vector<Simplicial_tree_node*>::size_type i = 0; i < cofaceNodes.size(); i++){
        cofaces->push_back(node_to_vector(cofaceNodes.at(i)));
    }
}

void Simplex_tree::get_cofacets(simplex_base *simplex, std::vector<simplex_base*> *cofacets)
{
    Simplicial_tree_node *simplexNode = find(simplex);
    std::vector<Simplicial_tree_node*> cofacetNodes;
    get_cofacets(simplexNode, &cofacetNodes);

    for (std::vector<Simplicial_tree_node*>::size_type i = 0; i < cofacetNodes.size(); i++){
        cofacets->push_back(node_to_vector(cofacetNodes.at(i)));
    }
}

inline bool Simplex_tree::contains(Simplex_tree::simplex_base *simplex)
{
    return find(simplex) != NULL;
}

Simplex_tree::Simplicial_tree_node* Simplex_tree::insert_simplex_in_tree(simplex_base *simplex)
{
    int dim = simplex->size() - 1;

    if (dim == 0){
        if (numberOfSimplices_ == 0) dictionaries_->push_back(new labels_dictionary());
        labels_dictionary *vertices = dictionaries_->front();
        if (vertices->find(simplex->at(0)) != vertices->end()) return NULL;
        Simplicial_tree_node *simplexNode = new Simplicial_tree_node(simplex->at(0), maxIndex_ + 1, 0, NULL);
        vertices->emplace(simplex->at(0), simplexNode);
        return simplexNode;
    }

    Simplicial_tree_node *it = dictionaries_->front()->at(simplex->at(0));
    for (int i = 1; i < dim; i++){
        it = it->get_children()->at(simplex->at(i));
    }

    if (it->get_children()->find(simplex->back()) != it->get_children()->end()) return NULL;

    if ((int)dictionaries_->size() == dim) dictionaries_->push_back(new labels_dictionary());
    Simplicial_tree_node *simplexNode = new Simplicial_tree_node(simplex->back(), maxIndex_ + 1, dim, it);
    it->get_children()->emplace(simplex->back(), simplexNode);
    insert_node_in_dictionary(simplexNode);

    return simplexNode;
}

void Simplex_tree::insert_node_in_dictionary(Simplicial_tree_node *simplex)
{
    labels_dictionary *dict = dictionaries_->at(simplex->get_dim());

    if (dict->find(simplex->get_label()) == dict->end()) {
        dict->emplace(simplex->get_label(), simplex);
        simplex->set_next(simplex);
        simplex->set_prev(simplex);
    } else {
        Simplicial_tree_node *head = dict->at(simplex->get_label());
        simplex->set_next(head->get_next());
        simplex->set_prev(head);
        head->get_next()->set_prev(simplex);
        head->set_next(simplex);
    }
}

void Simplex_tree::delete_node_from_dictionary(Simplicial_tree_node *simplex)
{
    labels_dictionary *dict = dictionaries_->at(simplex->get_dim());

    if (simplex->get_next() == simplex) dict->erase(simplex->get_label());
    else {
        simplex->get_next()->set_prev(simplex->get_prev());
        simplex->get_prev()->set_next(simplex->get_next());
        if (dict->at(simplex->get_label()) == simplex) dict->at(simplex->get_label()) = simplex->get_next();
    }
}

void Simplex_tree::remove_simplex(Simplicial_tree_node *simplex, std::vector<index> *removedIndices)
{
    std::vector<Simplicial_tree_node*> cofaces;
    get_cofaces(simplex, &cofaces);

    std::sort(cofaces.begin(), cofaces.end(), [](Simplicial_tree_node *n1, Simplicial_tree_node *n2){return n1->get_dim() > n2->get_dim();});

    for (std::vector<Simplicial_tree_node*>::size_type i = 0; i < cofaces.size(); i++){
        Simplicial_tree_node *toDelete = cofaces.at(i);
        if (removedIndices != NULL) removedIndices->push_back(toDelete->get_insertion_index());

        if (toDelete->get_dim() > 0) toDelete->get_parent()->get_children()->erase(toDelete->get_label());
        delete_node_from_dictionary(toDelete);

        numberOfSimplicesPerDim_->at(toDelete->get_dim())--;
        if (numberOfSimplicesPerDim_->at(toDelete->get_dim()) == 0) numberOfSimplicesPerDim_->pop_back();
        numberOfSimplices_--;

        delete toDelete;
    }
}

Simplex_tree::Simplicial_tree_node *Simplex_tree::find(simplex_base *simplex)
{
    if (simplex->empty() || dictionaries_->empty()) return NULL;

    vertex label = simplex->front();
    labels_dictionary *vertices = dictionaries_->front();

    if (vertices->find(label) == vertices->end()) return NULL;

    Simplicial_tree_node *it = vertices->at(label);
    for (simplex_base::size_type i = 1; i < simplex->size(); ++i) {
        label = simplex->at(i);
        labels_dictionary::iterator next = it->get_children()->find(label);

        if (next == it->get_children()->end()) return NULL;
        else it = next->second;
    }
    return it;
}

void Simplex_tree::get_cofaces(Simplicial_tree_node *simplex, std::vector<Simplicial_tree_node *> *cofaces)
{
    for (int d = dictionaries_->size() - 1; d >= simplex->get_dim(); d--){
        if (dictionaries_->at(d)->find(simplex->get_label()) != dictionaries_->at(d)->end()){
            Simplicial_tree_node *it = dictionaries_->at(d)->at(simplex->get_label());
            Simplicial_tree_node *startingNode = it;
            if (is_coface(it, simplex)) get_simplices_in_subtree(it, cofaces);
            it = it->get_next();
            while (it != startingNode){
                if (is_coface(it, simplex)) get_simplices_in_subtree(it, cofaces);
                it = it->get_next();
            }
        }
    }
}

void Simplex_tree::get_cofacets(Simplicial_tree_node *simplex, std::vector<Simplicial_tree_node*> *cofacets)
{
    Simplicial_tree_node *it = dictionaries_->at(simplex->get_dim() + 1)->at(simplex->get_label());
    Simplicial_tree_node *startingNode = it;
    if (is_coface(it, simplex)) cofacets->push_back(it);
    it = it->get_next();
    while (it != startingNode){
        if (is_coface(it, simplex)) cofacets->push_back(it);
        it = it->get_next();
    }

    for (auto childIt = simplex->get_children()->begin(); childIt != simplex->get_children()->end(); childIt++){
        cofacets->push_back(childIt->second);
    }
}

bool Simplex_tree::is_coface(Simplicial_tree_node *node, Simplicial_tree_node *simplex)
{
    Simplicial_tree_node *nodeIt = node;
    Simplicial_tree_node *simplexIt = simplex;

    while (simplexIt != NULL && nodeIt != NULL){
        if (simplexIt->get_label() == nodeIt->get_label()) simplexIt = simplexIt->get_parent();
        nodeIt = nodeIt->get_parent();
    }

    if (simplexIt == NULL) return true;
    else return false;
}

void Simplex_tree::get_simplices_in_subtree(Simplicial_tree_node *node, std::vector<Simplicial_tree_node *> *simplices)
{
    simplices->push_back(node);
    for (auto it = node->get_children()->begin(); it != node->get_children()->end(); it++){
        get_simplices_in_subtree(it->second, simplices);
    }
}

Simplex_tree::vertex Simplex_tree::get_smallest_star(vertex v, vertex u, std::queue<Simplicial_tree_node*> *qv, std::queue<Simplicial_tree_node*> *qu)
{
    labels_dictionary *childrenV;
    labels_dictionary *childrenU;
    labels_dictionary::iterator vit;
    labels_dictionary::iterator uit;
    std::queue<Simplicial_tree_node*> *tv = new std::queue<Simplicial_tree_node*>();
    std::queue<Simplicial_tree_node*> *tu = new std::queue<Simplicial_tree_node*>();
    std::vector<labels_dictionary*>::size_type currentHeightV = 0;
    std::vector<labels_dictionary*>::size_type currentHeightU = 0;
    Simplicial_tree_node *startRootV = dictionaries_->at(0)->at(v);
    Simplicial_tree_node *startRootU = dictionaries_->at(0)->at(u);
    Simplicial_tree_node *currentRootV = startRootV;
    Simplicial_tree_node *currentRootU = startRootU;
    Simplicial_tree_node *currentNodeV = startRootV;
    Simplicial_tree_node *currentNodeU = startRootU;
    bool currentNodeIsRootV = true;
    bool currentNodeIsRootU = true;
    bool continueV = true;
    bool continueU = true;

    while (continueV && continueU){
        qv->push(currentNodeV);
        if (!currentNodeIsRootV) tv->push(currentNodeV);
        qu->push(currentNodeU);
        if (!currentNodeIsRootU) tu->push(currentNodeU);

        continueV = get_smallest_star_find_next_node(v, u, currentHeightV, currentRootV, startRootV, currentNodeV, currentNodeIsRootV, childrenV, vit, tv);
        continueU = get_smallest_star_find_next_node(u, v, currentHeightU, currentRootU, startRootU, currentNodeU, currentNodeIsRootU, childrenU, uit, tu);
    }

    delete tv;
    delete tu;

    if (!continueU) return u;
    else return v;
}

bool Simplex_tree::get_smallest_star_find_next_node(vertex v, vertex toAvoid, std::vector<labels_dictionary*>::size_type &currentHeight, Simplicial_tree_node *&currentRoot,
                                                    Simplicial_tree_node *&startRoot, Simplicial_tree_node *&currentNode, bool &currentNodeIsRoot,
                                                    labels_dictionary *&children, labels_dictionary::iterator &it, std::queue<Simplicial_tree_node *> *&tail)
{
    if (currentNodeIsRoot){
        children = currentRoot->get_children();
        it = children->begin();
        currentNodeIsRoot = false;
    } else {
        it++;
    }

    if (it != children->end() && it->first == toAvoid) it++;
    while (it == children->end()){
        if (!tail->empty()){
            children = tail->front()->get_children();
            it = children->begin();
            tail->pop();
        } else {
            return get_smallest_star_find_next_root(v, toAvoid, currentHeight, currentRoot, startRoot, currentNode, currentNodeIsRoot);
        }
    }

    if (it->first == toAvoid){
        return get_smallest_star_find_next_node(v, toAvoid, currentHeight, currentRoot, startRoot, currentNode, currentNodeIsRoot, children, it, tail);
    }

    currentNode = it->second;
    return true;
}

bool Simplex_tree::get_smallest_star_find_next_root(vertex v, vertex toAvoid, std::vector<labels_dictionary*>::size_type &currentHeight, Simplicial_tree_node *&currentRoot,
                                                    Simplicial_tree_node *&startRoot, Simplicial_tree_node *&currentNode, bool &currentNodeIsRoot)
{
    if (currentRoot->get_next() == startRoot){
        currentHeight++;
        while (currentHeight < dictionaries_->size() && dictionaries_->at(currentHeight)->find(v) == dictionaries_->at(currentHeight)->end()) currentHeight++;
        if (currentHeight == dictionaries_->size()) return false;
        startRoot = dictionaries_->at(currentHeight)->at(v);
        currentRoot = startRoot;
    } else {
        currentRoot = currentRoot->get_next();
    }

    if (is_coface(currentRoot, dictionaries_->at(0)->at(toAvoid)))
        return get_smallest_star_find_next_root(v, toAvoid, currentHeight, currentRoot, startRoot, currentNode, currentNodeIsRoot);
    currentNode = currentRoot;
    currentNodeIsRoot = true;
    return true;
}

Simplex_tree::simplex_base *Simplex_tree::node_to_vector(Simplicial_tree_node *node)
{
    simplex_base *vector = new simplex_base();
    Simplicial_tree_node *nodeIt = node;
    while (nodeIt != NULL){
        vector->push_back(nodeIt->get_label());
        nodeIt = nodeIt->get_parent();
    }
    std::reverse(vector->begin(), vector->end());
    return vector;
}

Simplex_tree::Simplicial_tree_node *Simplex_tree::get_opposite_facet(Simplicial_tree_node *simplex, vertex v)
{
    if (simplex->get_parent() == NULL) return NULL;

    Simplicial_tree_node *trav = simplex;
    std::stack<vertex> tail;

    while (trav != NULL && trav->get_label() > v) {
        tail.push(trav->get_label());
        trav = trav->get_parent();
    }

    if (trav == NULL || trav->get_label() != v) return NULL;

    if (trav->get_parent() != NULL){
        trav = trav->get_parent();
    } else {
        trav = dictionaries_->at(0)->at(tail.top());
        tail.pop();
    }
    while (!tail.empty()){
        trav = trav->get_children()->at(tail.top());
        tail.pop();
    }

    return trav;
}

Simplex_tree::Simplicial_tree_node::Simplicial_tree_node(vertex label, index insertionIndex, int dim, Simplicial_tree_node *parent) :
    label_(label), insertionIndex_(insertionIndex), dim_(dim), parent_(parent), next_(this), prev_(this)
{
    children_ = new labels_dictionary();
}

Simplex_tree::Simplicial_tree_node::~Simplicial_tree_node()
{
    delete children_;
}

inline Simplex_tree::vertex Simplex_tree::Simplicial_tree_node::get_label() const
{
    return label_;
}

inline Simplex_tree::index Simplex_tree::Simplicial_tree_node::get_insertion_index() const
{
    return insertionIndex_;
}

inline Simplex_tree::Simplicial_tree_node *Simplex_tree::Simplicial_tree_node::get_parent() const
{
    return parent_;
}

inline void Simplex_tree::Simplicial_tree_node::set_parent(Simplicial_tree_node *value)
{
    parent_ = value;
}

inline Simplex_tree::labels_dictionary *Simplex_tree::Simplicial_tree_node::get_children() const
{
    return children_;
}

inline Simplex_tree::Simplicial_tree_node *Simplex_tree::Simplicial_tree_node::get_next() const
{
    return next_;
}

inline void Simplex_tree::Simplicial_tree_node::set_next(Simplicial_tree_node *value)
{
    next_ = value;
}

inline Simplex_tree::Simplicial_tree_node *Simplex_tree::Simplicial_tree_node::get_prev() const
{
    return prev_;
}

inline void Simplex_tree::Simplicial_tree_node::set_prev(Simplicial_tree_node *value)
{
    prev_ = value;
}

inline int Simplex_tree::Simplicial_tree_node::get_dim() const
{
    return dim_;
}

}
}

#endif // SIMPLEX_TREE_H
