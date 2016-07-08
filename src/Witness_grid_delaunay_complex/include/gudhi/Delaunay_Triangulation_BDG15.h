#ifndef Delaunay_Triangulation_BDG15_H_
#define Delaunay_Triangulation_BDG15_H_
#define CGAL_EIGEN3_ENABLED
#include <gudhi/Simplex_tree.h>
#include <gudhi/distance_functions.h>
#include <gudhi/Point.h>
#include <gudhi/reader_utils.h>
#include <boost/math/constants/constants.hpp>
#include <algorithm>
#include <utility>
#include <vector>
#include <list>
#include <set>
#include <queue>
#include <limits>
#include <ctime>
#include <cmath>
#include <iostream>
#include <CGAL/point_generators_d.h>
#include <CGAL/Cartesian_d.h>
#include <CGAL/Epick_d.h>
#include <CGAL/constructions_d.h>
#include <CGAL/Kernel_d/Hyperplane_d.h>
#include <CGAL/Origin.h>
#include <Good_links.h>
namespace Gudhi
{
namespace Delaunay_Triangulation_BDG15
{
template <class Simplicial_complex>
class Thick_relaxed_witness_complex
{
private:
    typedef CGAL::Cartesian_d<double> Kd;
    typedef Kd::Point_d Point;
    typedef Kd::Direction_d Direction;
    typedef Kd::Vector_d Vector;
    typedef Kd::Hyperplane_d Hyperplane;
    typedef std::vector<Point> Point_Vector;
    typedef std::vector< Vertex_handle > typeVectorVertex;
    typedef std::pair< Simplex_tree<>::Simplex_handle, bool > typePairSimplexBool;
    void subgen(std::vector<std::vector<int> > &C,std::vector<int> const &N,int const dim,int const index)//finds all possible d subsets of a "index" sized vector
    {
        if(dim>index)return;
        if(dim==1)
        {
            for(int i=0; i<index; i++)
            {
                std::vector<int> t;
                t.push_back(N[i]);
                C.push_back(t);
            }
            return;
        }
        else if(dim==index)
        {
            std::vector<int> t;
            for(int i=0; i<index; i++)
            {
                t.push_back(N[i]);
            }
            C.push_back(t);
            return;
        }
        else
        {
            std::vector<std::vector<int> > temp,temp2;
            subgen(temp,N,dim-1,index-1);
            for(int i=0; i<temp.size(); i++)
            {
                temp[i].push_back(N[index-1]);
                C.push_back(temp[i]);
            }
            subgen(temp2,N,dim,index-1);
            for(int i=0; i<temp2.size(); i++)
            {
                C.push_back(temp2[i]);
            }
            return;
        }
    }
    void full_cell_counter(std::vector<Point> &full_leaves,int &cnt,double const B, std::vector<Hyperplane> const &H,Point center, double len, double const eps)
    {
        int dim=center.dimension();
        if(cnt>B)return;
        else
        {
            bool all=1;//all is true if all hyperplanes are cutting the cell
            for(int pp=0; pp<H.size(); pp++)
            {
                Point far=center;
                Vector vfar;
                if(H[pp].has_on_positive_side(center)==0)//forming vfar vector
                {
                    vfar=H[pp].orthogonal_vector();
                }
                else
                {
                    vfar=H[pp].orthogonal_vector();
                    vfar*=-1;
                }
                std::vector<double> temp;
                for(int ppp=0; ppp<dim; ppp++) //forming far point
                {
                    if(vfar[ppp]>=0)
                    {
                        temp.push_back(len/2);
                    }
                    else
                    {
                        temp.push_back(-1*len/2);
                    }
                }
                vfar=Vector(dim,temp.begin(),temp.end());
                far=far+vfar;
                if(H[pp].has_on_positive_side(center)==H[pp].has_on_positive_side(far))
                {
                    all=0;
                    break;
                }
            }
            if(all)
            {
                cnt++;   //if it is full cell increase count
                if(len<=eps/sqrt(2))
                {
                    full_leaves.push_back(center);
                    return;
                }
                else
                {
                    std::vector<int> b;//first child
                    for(int pp=0; pp<dim; pp++)
                    {
                        b.push_back(-1);
                    }
                    do//loop over children
                    {
                        Vector vchild(dim,b.begin(),b.end());
                        vchild*=(len/4);
                        Point child=center+vchild;
                        full_cell_counter(full_leaves,cnt,B,H,child,len/2,eps);//call on its children
                    }
                    while(next_child(b,b.size())==1);
                    return;
                }
            }
            else
            {
                return;
            }
        }
    }
    bool next_child(std::vector<int> &b,int index)
    {
        if(index==1)
        {
            if(b[0]==-1)
            {
                b[0]=1;
                return 1;
            }
            else
            {
                return 0;
            }
        }
        else
        {
            if(b[index-1]==-1)
            {
                b[index-1]=1;
                return 1;
            }
            else
            {
                b[index-1]=-1;
                return next_child(b,index-1);
            }
        }
    }
    int factorial(int x)
    {
        const int X=x;
        for(int i=1; i<X; i++)
        {
            x=x*(X-i);
        }
        return x;
    }
    class sorter
    {
    public:
        Point_Vector center;
        sorter(std::vector<double> cent)
        {
            center.push_back(Point(cent.size(),cent.begin(),cent.end()));
        }
        bool operator()(Point a,Point b)
        {
            return (squared_distance(center[0],a)<squared_distance(center[0],b));
        }
    };
    const double pi = boost::math::constants::pi<double>();
public:
    //constructor
    Thick_relaxed_witness_complex(Point_Vector L,int dim,Simplicial_complex &sc,Simplicial_complex &prev,double lambda,double eps,bool check,bool fullleaf_deltacheck,std::vector<bool> &boundary,std::vector<bool> &perturbed, std::map<std::vector<int>,std::vector<Point> > &leaves_of_simplex)
    {
        std::map<std::vector<int>,std::vector<Point> > leaves_of_simplex_curr;//map to store full leaves based on simlex vertex vector as key
        double spar;
        double U_d=pow(pi,dim/2);//volume of unit d-ball
        if(dim%2==0)U_d=U_d/factorial(dim/2);
        else U_d=U_d*pow(2,dim/2+1)/factorial(dim);
        int nP=L.size()-2*dim;
        std::vector<double> org(dim,0);
        Point origin(dim,org.begin(),org.end());//origin point
        double far=0;
        for(int k=0; k<nP; k++)
        {
            if(squared_distance(L[k],origin)>far)far=squared_distance(L[k],origin);
        }
        double lambda_big=2*sqrt(dim)-far;//this lambda is not accurate
        double lambda_small=lambda;
        lambda=lambda_big;
        spar=squared_distance(L[0],L[1]);
        for(int k=0; k<L.size(); k++) //calculating sparsity
        {
            for(int kk=k+1; kk<L.size(); kk++)
            {
                if(squared_distance(L[k],L[kk])<spar)
                {
                    spar=squared_distance(L[k],L[kk]);
                }
            }
        }
        spar=sqrt(spar);
        double sparR=spar/lambda;//sparsity ratio
        const double delta=24*dim*eps/sparR;//parameter
        const double deltaR=delta/lambda;
        const double theta0=deltaR*sparR/(24*dim);
        const double B=U_d*(pow(4*dim/(theta0*sparR),dim))*(log(5*sqrt(dim)*lambda/eps));//bound on fullcells//to be drastically improved
        std::vector<std::vector<int> > neighbors(nP+2*dim);
        for(int i=0; i<L.size(); i++)
        {
            if(i<nP&&lambda_small<lambda_big)
            {
                lambda=lambda_small;
            }
            else
            {
                lambda=lambda_big;
            }
            double R=2*lambda+4*eps;
            double mini;
            int mindex;
            if(i==0)
            {
                mini=squared_distance(L[i],L[1]);
                mindex=1;
            }
            else
            {
                mini=squared_distance(L[i],L[0]);
                mindex=0;
            }
            //std::cout<<"Neighbors of "<<i;
            for(int j=0; j<L.size(); j++)//creating neighbors N(p)
            {
                if(squared_distance(L[i],L[j])<R*R&&i!=j)
                {
                    neighbors[i].push_back(j);
                    //std::cout<<j<<" ";
                }
                if(squared_distance(L[i],L[j])<mini&&i!=j)
                {
                    mini=squared_distance(L[i],L[j]);
                    mindex=j;
                }
            }
            if(neighbors[i].size()==0)
            {
                //std::cout<<mindex<<"' ";
                neighbors[i].push_back(mindex);
                neighbors[mindex].push_back(i);
                sort(neighbors[mindex].begin(),neighbors[mindex].end());
            }
            //std::cout<<std::endl;
        }
        for(int i=0; i<L.size(); i++)//i represents point p...looping over all points as in algorithm 3 in the paper
        {
            std::vector<std::vector<int> > C;//assume nbP>dim, C represents C_d(p)
            subgen(C,neighbors[i],dim,neighbors[i].size());
            for(int k=0; k<C.size(); k++)
            {
                C[k].push_back(i);
            }
            for(int sim=0; sim<C.size(); sim++) //checks if the simplices of C are to be included are not
            {
                bool fullcell_check=0;
                std::vector<Point> full_leaves;
                sort(C[sim].begin(),C[sim].end());
                bool ptb=0;
                for(int f=0; f<dim+1; f++)
                {
                    if(perturbed[C[sim][f]]==1)
                    {
                        ptb=1;
                        break;
                    }
                }
                if(ptb==0&&prev.num_simplices()!=0)//the simplex with same points i.e with same coordinates is in previous tree
                {
                    fullcell_check=1;
                    full_leaves=leaves_of_simplex[C[sim]];
                }
                if(fullcell_check==0)
                {
                    std::vector<double> cent;
                    for(int d=0; d<dim; d++)
                    {
                        double mini=L[C[sim][0]][d];
                        double maxi=L[C[sim][0]][d];
                        for(int s=1; s<dim+1; s++)
                        {
                            if(L[C[sim][s]][d]<mini)
                            {
                                mini=L[C[sim][s]][d];
                            }
                            if(L[C[sim][s]][d]>maxi)
                            {
                                maxi=L[C[sim][s]][d];
                            }
                        }
                        double c=(mini+maxi)/2;
                        cent.push_back(c);
                    }//center is obtained
                    //to create the hyperplanes
                    std::vector<Hyperplane> H;
                    for(int i1=0; i1<dim+1; i1++)
                    {
                        for(int j1=i1+1; j1<dim+1; j1++)
                        {
                            std::vector<double> midp;
                            Direction dir=(L[C[sim][i1]]-L[C[sim][j1]]).direction();

                            for(int d=0; d<dim; d++)
                            {
                                midp.push_back((L[C[sim][i1]][d]+L[C[sim][j1]][d])/2);
                            }
                            Point p(dim,midp.begin(),midp.end());
                            H.push_back(Hyperplane(p,dir));
                        }
                    }//to recurse over the box
                    int cnt=0;
                    Point center(dim,cent.begin(),cent.end());
                    full_cell_counter(full_leaves,cnt,B,H,center,4*lambda+8*eps,eps);
                    if(cnt>B||full_leaves.empty())continue;
                    else
                    {
                        fullcell_check=1;
                    }
                }
                if(fullcell_check==1)
                {
                    if(check==1) //check if the diameter of the full leaves is bounded and also the delta protection
                    {
                        double const dia_bound=16*sqrt(dim)*eps/(theta0*sparR);//to be drastically improved
                        double maxc,minc,dia;
                        for(int ss=0; ss<dim; ss++)
                        {
                            maxc=full_leaves[0][ss];
                            maxc=minc;
                            for(int s=1; s<full_leaves.size(); s++)
                            {
                                if(full_leaves[s][ss]>maxc)
                                {
                                    maxc=full_leaves[s][ss];
                                }
                                if(full_leaves[s][ss]<minc)
                                {
                                    minc=full_leaves[s][ss];
                                }
                            }
                            if(maxc-minc>dia)
                            {
                                dia=maxc-minc;
                            }
                        }
                        if(dia>dia_bound)
                        {
                            continue;
                        }
                        bool prot=1;
                        for(int f=0; prot==1&&f<full_leaves.size(); f++) //all full leaves
                        {
                            int s;
                            for(s=0; prot==1&&s<L.size(); s++) //over all landmarks
                            {
                                int ss;
                                for(ss=0; ss<dim+1; ss++)
                                {
                                    if(s==C[sim][ss])
                                    {
                                        break;
                                    }
                                }
                                if(ss!=dim+1)continue;//found a point in the simplex that is the same as point s
                                for(ss=0; prot==1&&ss<dim+1; ss++)
                                {
                                    if(sqrt(squared_distance(L[C[sim][ss]],full_leaves[f]))+delta-2*eps>=sqrt(squared_distance(L[s],full_leaves[f])))
                                    {
                                        prot=0;
                                        break;//found a violation of protection
                                    }
                                }
                            }
                        }
                        if(prot==0)
                        {
                            continue;
                        }
                    }

bool bndry=0;
                    for(int f=0; f<dim+1; f++)
                    {
                        if(C[sim][f]>=nP)
                        {
                            bndry=1;
                            break;
                        }
                    }
                    if(fullleaf_deltacheck==1&&bndry==0)
                    {
                        std::vector<double> cent;
                        for(int d=0; d<dim; d++)
                        {
                            double mini=L[C[sim][0]][d];
                            double maxi=L[C[sim][0]][d];
                            for(int s=1; s<dim+1; s++)
                            {
                                if(L[C[sim][s]][d]<mini)
                                {
                                    mini=L[C[sim][s]][d];
                                }
                                if(L[C[sim][s]][d]>maxi)
                                {
                                    maxi=L[C[sim][s]][d];
                                }
                            }
                            double c=(mini+maxi)/2;
                            cent.push_back(c);
                        }
                        Point_Vector sortL=L;
                        sort(sortL.begin(),sortL.end(),sorter(cent));
                        bool flop=0;
                        double maxc,minc,dia=0;
                        for(int ss=0; ss<dim; ss++)
                        {
                            maxc=full_leaves[0][ss];
                            maxc=minc;
                            for(int s=1; s<full_leaves.size(); s++)
                            {
                                if(full_leaves[s][ss]>maxc)
                                {
                                    maxc=full_leaves[s][ss];
                                }
                                if(full_leaves[s][ss]<minc)
                                {
                                    minc=full_leaves[s][ss];
                                }
                            }
                            if(maxc-minc>dia)
                            {
                                dia=maxc-minc;
                            }
                        }
                        double deltaW=2*dia;
                        for(int f=0; flop==0&&f<full_leaves.size(); f++)
                        {
                            for(int q=0; flop==0&&q<sortL.size(); q++)
                            {
                                int p;
                                for(p=0; p<dim+1; p++)
                                {
                                    if(sortL[q]==L[C[sim][p]])
                                    {
                                        break;
                                    }
                                }
                                if(p!=dim+1)continue;
                                for( p=0; flop==0&&p<dim+1; p++)
                                {
                                    if(squared_distance(sortL[q],full_leaves[f])+deltaW<squared_distance(L[C[sim][p]],full_leaves[f]))
                                    {
                                        flop=1;
                                    }
                                }
                            }
                        }
                        if(flop==1)continue;

                    }
                    //add the simplex to the tree
                    typeVectorVertex simplex=C[sim];
                    typePairSimplexBool returnValue =sc.insert_simplex_and_subfaces(simplex);
                    if (returnValue.second == true)
                    {
                        leaves_of_simplex_curr[C[sim]]=full_leaves;
                    }
                    //marking the boundary simplices
                    bndry=0;
                    for(int f=0; f<dim+1; f++)
                    {
                        if(C[sim][f]>=nP)
                        {
                            bndry=1;
                            break;
                        }
                    }
                    if(bndry==1)
                    {
                        for(int f=0; f<dim+1; f++)
                        {
                            boundary[C[sim][f]]=1;
                            //std::cout<<"boundary"<<C[sim][f]<<" ";
                        }
                        //std::cout<<std::endl;
                    }
                }
            }
        }
        leaves_of_simplex.clear();
        leaves_of_simplex=leaves_of_simplex_curr;
    }
};
//--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template<class Simplicial_complex>
class Delaunay_Triangulation_BDG15
{
private:
    typedef CGAL::Cartesian_d<double> Kd;
    typedef Kd::Point_d Point;
    typedef Kd::Vector_d Vector;
    typedef std::vector<Point> Point_Vector;
    void perturb(const Point_Vector &L, Point_Vector &L1, int index,int dim,double rho)
    {
        CGAL:: Random_points_in_ball_d<Vector> rp(dim,rho);
        L1[index]=(L[index]+(*rp++));
    }
public:
    //constructor
    Delaunay_Triangulation_BDG15(Simplex_tree<> &delt,Point_Vector & L, int const dim,double lambda, const double eps,double const rho)
    {
        int nP=L.size();
        std::vector<std::vector<Vertex_handle> > I(nP);
        std::map<std::vector<int>,std::vector<Point> > leaves_of_simplex;
        Point_Vector L1=L;
        const double lambda1=lambda+rho;//lambda1 is lambda'
        double spar=squared_distance(L[0],L[1]);
        for(int k=0; k<L.size(); k++) //calculating sparsity
        {

            for(int kk=k+1; kk<L.size(); kk++)
            {
                if(squared_distance(L[k],L[kk])<spar)
                {
                    spar=squared_distance(L[k],L[kk]);
                }
            }
        }
        spar=sqrt(spar);
        double sparR=spar/lambda;
        double sparR1=sparR/3;
        const double R=(5+(3*sparR)/2)*lambda;
        for(int i=0; i<nP; i++)// creating I(p) set for every point
        {
            for(int j=0; j<nP; j++)
            {
                if(squared_distance(L[i],L[j])<R*R&&i!=j)
                {
                    I[i].push_back(j);
                }
            }
        }
        for(int k=0; k<dim; k++)//inserting points on dual of cube
        {
            std::vector<double> A_p,A_n;//inserting points on dual of cube
            double const dual=2*sqrt(dim);
            for(int i=0; i<dim; i++)
            {
                if(i==k)
                {
                    A_p.push_back(dual);
                }
                else
                {
                    A_p.push_back(0);
                }
            }
            L1.push_back(Point(dim,A_p.begin(),A_p.end()));
            for(int i=0; i<dim; i++)
            {
                if(i==k)
                {
                    A_n.push_back(-1*dual);
                }
                else
                {
                    A_n.push_back(0);
                }
            }
            L1.push_back(Point(dim,A_n.begin(),A_n.end()));
        }//point set P U Q created.
        Simplex_tree<> temp_tree;
        bool bad=0;
        int badcnt=0;
        do
        {
            Simplex_tree<> dt;
            badcnt=0;
            bad=0;
            std::vector<bool> boundary(nP+2*dim,0);
            std::vector<bool> perturbed(nP,0);

            Thick_relaxed_witness_complex<Simplicial_complex> T(L1,dim,dt,temp_tree,lambda,eps,true,true,boundary,perturbed,leaves_of_simplex);
            Simplex_tree<> temp_tree=dt;
            Good_links<Simplicial_complex> temp(temp_tree);
            for(int i=0; i<nP; i++)
            {
                //std::cout<<i<<" :"<<boundary[i]<<"  ";
                if(!boundary[i])
                {

                    if(!temp.has_good_link(i))//perturb
                    {
                        perturb(L,L1,i,dim,rho);
                        perturbed[i]=1;
                        for(int j=0; j<I[i].size(); j++)
                        {
                            perturb(L,L1,I[i][j],dim,rho);
                            perturbed[I[i][j]]=1;
                        }
                        bad=1;
                        badcnt++;
                    }
                }
                //std::cout<<std::endl;
            }
            lambda=lambda1;
            sparR=sparR1;
            std::cout<<badcnt<<"badcnt"<<std::endl;
            if(bad==0)
            {
                for (auto simplex : dt.complex_simplex_range())
                {
                    std::vector<int> simp;
                    for (auto vertex : dt.simplex_vertex_range(simplex))
                    {
                        simp.push_back(vertex);
                    }
                    delt.insert_simplex_and_subfaces(simp);
                }
		for(int l=0;l<L.size();l++)
		{
			L[l]=L1[l];		
		}
            }
        }
        while(badcnt>0);

    }
};
}
}
#endif // Delaunay_Triangulation_BDG15_H_
