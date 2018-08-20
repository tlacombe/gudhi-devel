#ifndef SIMPLEX_TREE_H
#define SIMPLEX_TREE_H

#include <vector>
#include <unordered_map>
#include <queue>
#include <stack>
#include <algorithm>

namespace Gudhi {
namespace tower_to_filtration {

class Simplex_tree
{
public:
    class Node;
    using vertex = double;
    using index = double;
    using simplex_base = std::vector<vertex>;
    using label_dictionary = std::unordered_map<vertex, Node*>;

    Simplex_tree();
    ~Simplex_tree();

    class Node
    {
    public:
        Node(vertex label, index insertionIndex, int dim, Node *parent);
        ~Node();

        vertex get_label() const;
        index get_insertion_index() const;

        Node* get_parent() const;
        void set_parent(Node *value);
        label_dictionary *get_children() const;

        Node* get_next() const;
        void set_next(Node *value);
        Node* get_prev() const;
        void set_prev(Node *value);

        int get_dim() const;

    private:
        vertex label_;
        index insertionIndex_;
        int dim_;
        Node *parent_;
        label_dictionary *children_;
        Node *next_;
        Node *prev_;
    };

    /* Important for tower_converter.h */

    Node* insert_simplex(simplex_base &numVertices);
    void remove_simplex(simplex_base &simplex);
    void remove_simplex(simplex_base &simplex, std::vector<index> *removedIndices);
    vertex get_smallest_closed_star(vertex v, vertex u, std::vector<simplex_base*> *closedStar);   //returns vertex with smallest closed star; closedStar is ordered
    index get_boundary(simplex_base &simplex, std::vector<index> *boundary);    //ordered by increasing insertion numbers
    double get_size() const;
    double get_max_index() const;

    /* Other */

    //index contract_vertices(vertex v, vertex u, double timestamp, std::vector<std::vector<index>*> *boundaries, std::vector<index> *&inactiveInsertionNumbers);
    double get_max_size() const;
    int get_max_dimension() const;
    void get_cofaces(simplex_base &simplex, std::vector<simplex_base*> *cofaces);
    void get_cofacets(simplex_base &simplex, std::vector<simplex_base*> *cofacets);
    bool contains(simplex_base &simplex);

private:
    std::vector<label_dictionary*> *dictionaries_; // circular list for each label at each height
    std::unordered_map<double, vertex> *verticesTranslation_;
    double numberOfSimplices_;
    double maxIndex_;
    double maxSize_;
    int maxDim_;

    Node* insert_simplex_in_tree(simplex_base &simplex);
    void insert_node_in_dictionary(Node *simplex);
    void delete_node_from_dictionary(Node *simplex);
    void remove_simplex(Node *simplex, std::vector<index> *removedIndices);
    Node* find(simplex_base &simplex);
    void get_cofaces(Node *simplex, std::vector<Node*> *cofaces);
    void get_cofacets(Node *simplex, std::vector<Node*> *cofacets);
    bool is_coface(Node *node, Node *simplex);
    void get_simplices_in_subtree(Node *node, std::vector<Node*> *simplices);
    vertex get_smallest_star(vertex v, vertex u, std::queue<Node*> *qv, std::queue<Node*> *qu);
    bool get_smallest_star_find_next_node(vertex v, vertex toAvoid, std::vector<label_dictionary*>::size_type &currentHeight,
                                          Node *&currentRoot, Node *&startRoot, Node *&currentNode, bool &currentNodeIsRoot,
                                          label_dictionary *&children, label_dictionary::iterator &it, std::queue<Node*> *&tail);
    bool get_smallest_star_find_next_root(vertex v, vertex toAvoid, std::vector<label_dictionary*>::size_type &currentHeight,
                                          Node *&currentRoot, Node *&startRoot, Node *&currentNode, bool &currentNodeIsRoot);
    simplex_base* node_to_vector(Node *node);
    Node* get_opposite_facet(Node *simplex, vertex v);
};

inline Simplex_tree::Simplex_tree() : numberOfSimplices_(0), maxIndex_(-1), maxSize_(0), maxDim_(0)
{
    verticesTranslation_ = new std::unordered_map<double, vertex>();
    dictionaries_ = new std::vector<label_dictionary*>();
}

inline Simplex_tree::~Simplex_tree()
{
    delete verticesTranslation_;
    for (std::vector<label_dictionary*>::size_type i = 0; i < dictionaries_->size(); i++){
        for (label_dictionary::iterator it = dictionaries_->at(i)->begin(); it != dictionaries_->at(i)->end(); it++) delete it->second;
        delete dictionaries_->at(i);
    }
    delete dictionaries_;
}

inline Simplex_tree::Node* Simplex_tree::insert_simplex(simplex_base &numVertices)
{
    Node *simplexNode = insert_simplex_in_tree(numVertices);
    if (simplexNode == NULL) return simplexNode;

    int dim = numVertices.size() - 1;

    numberOfSimplices_++;
    maxIndex_++;
    if (maxSize_ < numberOfSimplices_) maxSize_ = numberOfSimplices_;
    if (maxDim_ < dim) maxDim_ = dim;

    return simplexNode;
}

inline void Simplex_tree::remove_simplex(simplex_base &simplex)
{
    Node *simplexNode = find(simplex);
    remove_simplex(simplexNode, NULL);
}

inline void Simplex_tree::remove_simplex(simplex_base &simplex, std::vector<index> *removedIndices)
{
    Node *simplexNode = find(simplex);
    remove_simplex(simplexNode, removedIndices);
}

inline Simplex_tree::vertex Simplex_tree::get_smallest_closed_star(vertex v, vertex u, std::vector<simplex_base *> *closedStar)
{
    std::queue<Node*> qv;
    std::queue<Node*> qu;
    Node *s;

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

inline Simplex_tree::index Simplex_tree::get_boundary(simplex_base &simplex, std::vector<index> *boundary)
{
    Node *simplexNode = find(simplex);
    if (simplexNode->get_parent() == NULL) return simplexNode->get_insertion_index();

    simplex_base tail;
    tail.push_back(simplexNode->get_label());
    Node *it = simplexNode->get_parent();
    simplex_base::reverse_iterator tailIt;
    Node *itb;

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

inline void Simplex_tree::get_cofaces(simplex_base &simplex, std::vector<simplex_base *> *cofaces)
{
    Node *simplexNode = find(simplex);
    std::vector<Node*> cofaceNodes;
    get_cofaces(simplexNode, &cofaceNodes);

    for (std::vector<Node*>::size_type i = 0; i < cofaceNodes.size(); i++){
        cofaces->push_back(node_to_vector(cofaceNodes.at(i)));
    }
}

inline void Simplex_tree::get_cofacets(simplex_base &simplex, std::vector<simplex_base*> *cofacets)
{
    Node *simplexNode = find(simplex);
    std::vector<Node*> cofacetNodes;
    get_cofacets(simplexNode, &cofacetNodes);

    for (std::vector<Node*>::size_type i = 0; i < cofacetNodes.size(); i++){
        cofacets->push_back(node_to_vector(cofacetNodes.at(i)));
    }
}

inline bool Simplex_tree::contains(simplex_base &simplex)
{
    return find(simplex) != NULL;
}

inline Simplex_tree::Node* Simplex_tree::insert_simplex_in_tree(simplex_base &simplex)
{
    int dim = simplex.size() - 1;

    if (dim == 0){
        if (numberOfSimplices_ == 0) dictionaries_->push_back(new label_dictionary());
        label_dictionary *vertices = dictionaries_->front();
	if (vertices->find(simplex.at(0)) != vertices->end()) return NULL;
	Node *simplexNode = new Node(simplex.at(0), maxIndex_ + 1, 0, NULL);
	vertices->emplace(simplex.at(0), simplexNode);
        return simplexNode;
    }

    Node *it = dictionaries_->front()->at(simplex.at(0));
    for (int i = 1; i < dim; i++){
	it = it->get_children()->at(simplex.at(i));
    }

    if (it->get_children()->find(simplex.back()) != it->get_children()->end()) return NULL;

    if ((int)dictionaries_->size() == dim) dictionaries_->push_back(new label_dictionary());
    Node *simplexNode = new Node(simplex.back(), maxIndex_ + 1, dim, it);
    it->get_children()->emplace(simplex.back(), simplexNode);
    insert_node_in_dictionary(simplexNode);

    return simplexNode;
}

inline void Simplex_tree::insert_node_in_dictionary(Node *simplex)
{
    label_dictionary *dict = dictionaries_->at(simplex->get_dim());

    if (dict->find(simplex->get_label()) == dict->end()) {
        dict->emplace(simplex->get_label(), simplex);
        simplex->set_next(simplex);
        simplex->set_prev(simplex);
    } else {
        Node *head = dict->at(simplex->get_label());
        simplex->set_next(head->get_next());
        simplex->set_prev(head);
        head->get_next()->set_prev(simplex);
        head->set_next(simplex);
    }
}

inline void Simplex_tree::delete_node_from_dictionary(Node *simplex)
{
    label_dictionary *dict = dictionaries_->at(simplex->get_dim());

    if (simplex->get_next() == simplex) dict->erase(simplex->get_label());
    else {
        simplex->get_next()->set_prev(simplex->get_prev());
        simplex->get_prev()->set_next(simplex->get_next());
        if (dict->at(simplex->get_label()) == simplex) dict->at(simplex->get_label()) = simplex->get_next();
    }
}

inline void Simplex_tree::remove_simplex(Node *simplex, std::vector<index> *removedIndices)
{
    std::vector<Node*> cofaces;
    get_cofaces(simplex, &cofaces);

    std::sort(cofaces.begin(), cofaces.end(), [](Node *n1, Node *n2){return n1->get_dim() > n2->get_dim();});

    for (std::vector<Node*>::size_type i = 0; i < cofaces.size(); i++){
        Node *toDelete = cofaces.at(i);
        if (removedIndices != NULL) removedIndices->push_back(toDelete->get_insertion_index());

        if (toDelete->get_dim() > 0) toDelete->get_parent()->get_children()->erase(toDelete->get_label());
        delete_node_from_dictionary(toDelete);

        numberOfSimplices_--;

        delete toDelete;
    }
}

inline Simplex_tree::Node *Simplex_tree::find(simplex_base &simplex)
{
    if (simplex.empty() || dictionaries_->empty()) return NULL;

    vertex label = simplex.front();
    label_dictionary *vertices = dictionaries_->front();

    if (vertices->find(label) == vertices->end()) return NULL;

    Node *it = vertices->at(label);
    for (simplex_base::size_type i = 1; i < simplex.size(); ++i) {
	label = simplex.at(i);
        label_dictionary::iterator next = it->get_children()->find(label);

        if (next == it->get_children()->end()) return NULL;
        else it = next->second;
    }
    return it;
}

inline void Simplex_tree::get_cofaces(Node *simplex, std::vector<Node *> *cofaces)
{
    for (int d = dictionaries_->size() - 1; d >= simplex->get_dim(); d--){
        if (dictionaries_->at(d)->find(simplex->get_label()) != dictionaries_->at(d)->end()){
            Node *it = dictionaries_->at(d)->at(simplex->get_label());
            Node *startingNode = it;
            if (is_coface(it, simplex)) get_simplices_in_subtree(it, cofaces);
            it = it->get_next();
            while (it != startingNode){
                if (is_coface(it, simplex)) get_simplices_in_subtree(it, cofaces);
                it = it->get_next();
            }
        }
    }
}

inline void Simplex_tree::get_cofacets(Node *simplex, std::vector<Node*> *cofacets)
{
    Node *it = dictionaries_->at(simplex->get_dim() + 1)->at(simplex->get_label());
    Node *startingNode = it;
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

inline bool Simplex_tree::is_coface(Node *node, Node *simplex)
{
    Node *nodeIt = node;
    Node *simplexIt = simplex;

    while (simplexIt != NULL && nodeIt != NULL){
        if (simplexIt->get_label() == nodeIt->get_label()) simplexIt = simplexIt->get_parent();
        nodeIt = nodeIt->get_parent();
    }

    if (simplexIt == NULL) return true;
    else return false;
}

inline void Simplex_tree::get_simplices_in_subtree(Node *node, std::vector<Node *> *simplices)
{
    simplices->push_back(node);
    for (auto it = node->get_children()->begin(); it != node->get_children()->end(); it++){
        get_simplices_in_subtree(it->second, simplices);
    }
}

inline Simplex_tree::vertex Simplex_tree::get_smallest_star(vertex v, vertex u, std::queue<Node*> *qv, std::queue<Node*> *qu)
{
    label_dictionary *childrenV;
    label_dictionary *childrenU;
    label_dictionary::iterator vit;
    label_dictionary::iterator uit;
    std::queue<Node*> *tv = new std::queue<Node*>();
    std::queue<Node*> *tu = new std::queue<Node*>();
    std::vector<label_dictionary*>::size_type currentHeightV = 0;
    std::vector<label_dictionary*>::size_type currentHeightU = 0;
    Node *startRootV = dictionaries_->at(0)->at(v);
    Node *startRootU = dictionaries_->at(0)->at(u);
    Node *currentRootV = startRootV;
    Node *currentRootU = startRootU;
    Node *currentNodeV = startRootV;
    Node *currentNodeU = startRootU;
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

inline bool Simplex_tree::get_smallest_star_find_next_node(vertex v, vertex toAvoid, std::vector<label_dictionary*>::size_type &currentHeight, Node *&currentRoot,
                                                    Node *&startRoot, Node *&currentNode, bool &currentNodeIsRoot,
                                                    label_dictionary *&children, label_dictionary::iterator &it, std::queue<Node*> *&tail)
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

inline bool Simplex_tree::get_smallest_star_find_next_root(vertex v, vertex toAvoid, std::vector<label_dictionary*>::size_type &currentHeight, Node *&currentRoot,
                                                    Node *&startRoot, Node *&currentNode, bool &currentNodeIsRoot)
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

inline Simplex_tree::simplex_base *Simplex_tree::node_to_vector(Node *node)
{
    simplex_base *vector = new simplex_base();
    Node *nodeIt = node;
    while (nodeIt != NULL){
        vector->push_back(nodeIt->get_label());
        nodeIt = nodeIt->get_parent();
    }
    std::reverse(vector->begin(), vector->end());
    return vector;
}

inline Simplex_tree::Node *Simplex_tree::get_opposite_facet(Node *simplex, vertex v)
{
    if (simplex->get_parent() == NULL) return NULL;

    Node *trav = simplex;
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

inline Simplex_tree::Node::Node(vertex label, index insertionIndex, int dim, Node *parent) :
    label_(label), insertionIndex_(insertionIndex), dim_(dim), parent_(parent), next_(this), prev_(this)
{
    children_ = new label_dictionary();
}

inline Simplex_tree::Node::~Node()
{
    delete children_;
}

inline Simplex_tree::vertex Simplex_tree::Node::get_label() const
{
    return label_;
}

inline Simplex_tree::index Simplex_tree::Node::get_insertion_index() const
{
    return insertionIndex_;
}

inline Simplex_tree::Node *Simplex_tree::Node::get_parent() const
{
    return parent_;
}

inline void Simplex_tree::Node::set_parent(Node *value)
{
    parent_ = value;
}

inline Simplex_tree::label_dictionary *Simplex_tree::Node::get_children() const
{
    return children_;
}

inline Simplex_tree::Node *Simplex_tree::Node::get_next() const
{
    return next_;
}

inline void Simplex_tree::Node::set_next(Node *value)
{
    next_ = value;
}

inline Simplex_tree::Node *Simplex_tree::Node::get_prev() const
{
    return prev_;
}

inline void Simplex_tree::Node::set_prev(Node *value)
{
    prev_ = value;
}

inline int Simplex_tree::Node::get_dim() const
{
    return dim_;
}

}
}

#endif // SIMPLEX_TREE_H
