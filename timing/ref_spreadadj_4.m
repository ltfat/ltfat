function cadj=ref_spreadadj_4(coef)
%REF_SPREADADJ_4  Symbol of adjoint spreading function.
%   Usage: cadj=ref_spreadadj_4(c,number);
%
%   Development version by FJ for comparison of different implementations
%   cadj=SPREADADJ(c) will compute the symbol cadj of the spreading 
%   operator that is the adjoint of the spreading operator with symbol c. 
%
%   Improved implementation of the direct formula with higher memory
%   needs (but better speed)
%
%   This implementation uses the same improvements as case 3, but 
%   also avoid the use of loop by using matrix pointwise 
%   multiplication, which improves the computation time, but also 
%   requires more memory due to the contruction of the (L-1)x(L-1) 
%   matrix temp
%

L=size(coef,1);

cadj=zeros(L);
        
% Proceesing for ii==0 or jj==0
cadj(1,1)=conj(coef(1,1));
cadj(2:end,1)=conj(coef(end:-1:2,1));
cadj(1,2:end,1)=conj(coef(1,end:-1:2));

% Processing for ii~=0 and jj~=0

% Precomputation for exponential term

% Optimization note : As said in note of case 3 ,we are computing 
% the Lth root of unity which have special properties and symetries 
% that could be exploited to highly reduce this computation
temp2=exp((-i*2*pi/L)*(0:L-1));

% Optimization note : Here we are computing indexes for all the
% exponential terms, which leads to a highly structured matrix
% which strcture can be formalized (notably it is symetric) and
% used to reduce the computation cost
temp=mod((1:L-1)'*(1:L-1),L)+1;


% Optimization note : Finaly we construct the matrix containing all
% the needed exponential terms. 
% This matrix is known as the DFT matrix and appears in the matrix 
% formulation of the dft, but in our case it is used for pointwise 
% multiplication instead of matrix mulitplication.
% There might be optimized algorithms for the computation of this
% matrix described in the context of fft.
temp=temp2(temp);

cadj(2:L,2:L)=conj(coef(L:-1:2,L:-1:2)).*temp;




