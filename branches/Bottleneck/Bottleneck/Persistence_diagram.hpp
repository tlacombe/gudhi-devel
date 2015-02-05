#pragma once
#include <vector>
#include <utility>
#include <cmath>

// Diagram_point is the type of the persistence diagram's points
typedef typename std::pair<double,double> Diagram_point;

// Intern_container is the private type used to store the persistence diagram's points
typedef std::vector<Diagram_point> Intern_container;

// Diagram_point_iterator is the type used to iterate over the persistence diagram's points
typedef typename Intern_container::const_iterator Diagram_point_iterator;

int null_point_index();

double distance(Diagram_point p1, Diagram_point p2);

// Persistence_diagram is the interface beetwen any external representation of persistence diagrams and
// the bottleneck distance computation. An interface is necessary to ensure basic functions complexity.
class Persistence_diagram{

private:
    Intern_container points;

public:
    // Container is the type of any external representation of persistence diagrams. It have to have an iterator over points,
    // which have to have fields first (for birth) and second (for death).
    template<typename Container>
    Persistence_diagram(Container& c, double e = 0.);
    Persistence_diagram();
    Diagram_point get_point(int i) const;
    int size() const;
    Diagram_point_iterator begin() const;
    Diagram_point_iterator end() const;
};



inline int null_point_index()
{
    return -1;
}

inline double distance(Diagram_point p1, Diagram_point p2)
{
    return std::max(std::fabs(p1.first - p2.first), std::fabs(p1.second - p2.second));
}

template<typename Container>
Persistence_diagram::Persistence_diagram(Container& c, double e)
 : points()
{
    for(auto it = c.begin(); it != c.end(); ++it)
        if(it->second - it->first > e)
            points.emplace_back(it->first, it->second);
}

Persistence_diagram::Persistence_diagram()
    : points()
{}

inline Diagram_point Persistence_diagram::get_point(int i) const
{
    return points.at(i);
}

inline int Persistence_diagram::size() const
{
    return points.size();
}

inline Diagram_point_iterator Persistence_diagram::begin() const
{
    return points.cbegin();
}

inline Diagram_point_iterator Persistence_diagram::end() const
{
    return points.cend();
}
