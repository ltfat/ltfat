function cadj=ref_spreadadj(coef);
%REF_SPREADADJ  Symbol of adjoint preading function.
%   Usage: cadj=ref_spreadadj(c);
%
%   cadj=SPREADADJ(c) will compute the symbol cadj of the spreading
%   operator that is the adjoint of the spreading operator with symbol c.
%
%   The algorithm converts the symbol to the matrix representation,
%   adjoints its, and finds its spreading function.
  

L=size(coef,1);

T=tfmat('spread',coef);
cadj=spreadfun(T');



