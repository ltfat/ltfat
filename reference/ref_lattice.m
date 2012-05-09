function [lat]=ref_lattice(V,L);
%REF_LATTICE  List of lattice points.
%   Usage:  lat=ref_lattice(V,L)  
%
%   Returns the lattice given by av and bv.
%   The output format is a 2xMxN matrix, where
%   each column is a point on the lattice.
%
%   The lattice must be in lower triangular Hermite normal form.
  
a=V(1,1);
b=V(2,2);
s=V(2,1);
M=abs(L/b);
N=abs(L/a);

% Create lattice.
lattice=zeros(2,M*N);
for n=0:N-1
  for m=0:M-1
    soffset=mod(s*n,b);
    % Determine gridpoint in rectangular coordinates.
    %lat(:,m+n*M+1) = V(:,1)*n+V(:,2)*m;
    lat(1,m+n*M+1) = n*a;
    lat(2,m+n*M+1) = m*b+soffset;
  end;
end;

% Mod' the lattice.
lat=mod(lat,L);


