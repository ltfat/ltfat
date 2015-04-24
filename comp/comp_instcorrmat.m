function R = comp_instcorrmat(f, g, Ls);
%COMP_INSTCORRMAT Compute instantaneous correlation matrix
%   Usage R = comp_instcorrmat(f, g, Ls);
%
%   Input parameters:
%         f,g    : Input vectors
%	  Ls   	 : length of vectors.
%
%   Output parameters:
%         R      : Instantaneous correlation matrix.
%
%

% AUTHOR: Jordy van Velthoven

if ~all(size(f)==size(g))
  error('%s: f and g must have the same size.', upper(mfilename));
end

if ~all(mod(Ls,2) == 0)
 f = postpad(f, Ls+1);
 g = postpad(g, Ls+1);
end
 
	
R = zeros(Ls,Ls);

for l = 0 : Ls-1;
   m = -min([Ls-l, l, round(Ls/2)-1]) : min([Ls-l, l, round(Ls/2)-1]);
   R(mod(Ls+m,Ls)+1, l+1) =  f(mod(l+m, Ls)+1).*conj(g(mod(l-m, Ls)+1));
end

R = R(1:Ls, 1:Ls);
