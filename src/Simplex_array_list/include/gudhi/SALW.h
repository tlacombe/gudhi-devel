#ifndef SALW_H
#define SALW_H

#include <gudhi/LSAL.h>

namespace Gudhi {




class SALW {

private :
    SAL K;

public :
    typedef Simplex Simplex_handle;

    typedef Vertex Vertex_handle;

    typedef std::list<Simplex> Complex_simplex_range; //Iterator over the simplices of the complex, in an arbitrary order. More...

    typedef Simplex Simplex_vertex_range; //Iterator over vertices of a simplex. More...

    typedef int Filtration_value;

     Simplex_handle null_simplex(){
         return Simplex();
     }

   /*  Complex_simplex_range complex_simplex_range(){
         return K.max_cofaces(Simplex());
     }
*/
     Simplex_vertex_range simplex_vertex_range (Simplex_handle const &simplex){
        return simplex;
     }


     // template <typename Input_vertex_range>
     // void insert_simplex (Input_vertex_range const &vertex_range){
     //     Simplex s(vertex_range.begin(),vertex_range.end());
     //     K.insert_max(s);
     // }

     template <typename Input_vertex_range>
     void insert_simplex (Input_vertex_range const &vertex_range, double filtration_value = 0){
         Simplex s(vertex_range.begin(),vertex_range.end());
         K.insert_max(s);
     }

  
     template<typename Input_vertex_range>
     Simplex_handle find (Input_vertex_range const &vertex_range){
         Simplex s(vertex_range.begin(),vertex_range.end());
         return K.membership(s) ? s : Simplex();
     }

     std::size_t num_simplices(){
        return K.size();
     }

     std::size_t num_vertices(){
        return K.num_vertices();
     }

     void set_dimension(int k){
     }

     double filtration(Simplex_handle sh){
       return 0;
     }
};
    

} //namespace Gudhi

#endif /* SALW_H */
