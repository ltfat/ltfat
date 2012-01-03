function f=idft(f,N,dim)
%IDFT  Inverse DFT
%   Usage: f=idft(f);
%          f=idft(f,N,dim);
%
%   This function computes a normalized inverse discrete Fourier transform.
%   This is nothing but a scaled version of the output from IFFT. The
%   function takes exactly the same arguments as IFFT. See the help on IFFT
%   for a thorough description.
%
%   See also:  dft

%   AUTHOR : Peter Soendergaard
%   TESTING: OK
%   REFERENCE: OK

error(nargchk(1,3,nargin));

if nargin<3
  dim=[];  
end;

if nargin<2
  N=[];
end;

[f,N,Ls,W,dim,permutedsize,order]=assert_sigreshape_pre(f,N,dim,'IDFT');

% Force IFFT along dimension 1, since we have permuted the dimensions
% manually
f=ifft(f,N,1)*sqrt(N);

f=assert_sigreshape_post(f,dim,permutedsize,order);

%OLDFORMAT
