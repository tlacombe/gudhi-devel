// Code zigzag


Only coefficients in Z2

birth_vector_   //maintains the b and d orders.

// The Matrix:

std::list< matrix_chain >   matrix_; // 0 ... m-1
class matrix_chain 
{
	...
  Column            * column_      ;
  Row               * row_         ;
  matrix_chain      * paired_col_  ; //\in F -> NULL, \in H -> g, \in G -> h
  SimplexKey          birth_       ; //\in F -> b, \in H -> -2 \in G -> -1
  SimplexKey          lowest_idx_  ; //to update...
}

// todo.  deal with the maximal dimension
cpx_->dimension(zzsh) < max_dim_+12

//todo check the new implementation of the pairing in surjective diamond. See article

//todo represent directly col_ in matrix_chain, instead of column_->col_. Can we do that?

//todo check the order of simplices and the order enforced by col_bsh, and in the col_fp in surjective diamond

//todo turn the public access, that is everywhere, into private.

//todo template the SimplexKey in matrix_ ?

//todo Documentation of compute_zigzag_persistence + syntax like compute_persistence in GUDHI.

//use a pool ?

//todo interval compare and interval display policies (see Persistent_cohomology.h)