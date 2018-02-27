function R = UpdateCholesky(R,P,Pt,index,m)
% UpdateCholesky: Updates a cholesky decomposition matrix R'R=A 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Usage
% R = UpdateCholesky(R,P,Pt,index);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input
%   Mandantory:
%               R       upper triangular matrix, such that R'*R =
%                       PT(index(1:end-1))*P((1:end-1))
%               P       is a function handle (type "help function_format" 
%                       for more information)
%               Pt      is a function handle
%               index is the set of all non-zero indices, must be 1 larger
%                       than dimension of R.
%               m       dimension of space from which P maps
%
% Outputs :     R is traingular matrix such that R'*R = PT(index)*P(index)

[n n]=size(R);
li=length(index);
if ~(n+1==li)
    error('Incorect index length or size of R.')
end
linsolve_options_transpose.UT = true;
linsolve_options_transpose.TRANSA = true;

mask = zeros(m,1);
mask(index(end))=1;
new_vector = P(mask); 

if li==1
    R=sqrt(new_vector'*new_vector);
else
    Pt_new_vector   = Pt(new_vector);
    new_col         = linsolve(R, Pt_new_vector(index(1:end-1)),linsolve_options_transpose);
    R_ii            = sqrt(new_vector'*new_vector - new_col'*new_col);
    R               = [R new_col; zeros(1, size(R,2)) R_ii];
end
    