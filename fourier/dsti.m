function c=dsti(f,L,dim)
%DSTI  Discrete Sine Transform type I
%   Usage:  c=dsti(f);
%           c=dsti(f,L);
%           c=dsti(f,[],dim);
%           c=dsti(f,L,dim);
%
%   `dsti(f)` computes the discrete sine transform of type I of the
%   input signal *f*. If *f* is multi-dimensional, the transformation is
%   applied along the first non-singleton dimension.
%
%   `dsti(f,L)` zero-pads or truncates *f* to length *L* before doing the
%   transformation.
%
%   `dsti(f,[],dim)` or `dsti(f,L,dim)` applies the transformation along
%   dimension *dim*.
%
%   The transform is real (output is real if input is real) and orthonormal.
%
%   This transform is its own inverse.
%
%   Let f be a signal of length *L* and let `c=dsti(f)`. Then 
%
%   ..                         L-1
%     c(n+1) = sqrt(2/(L+1)) * sum sin(pi*(n+1)*(m+1)/(L+1)) 
%                              m=0 
%   .. math:: c\left(n+1\right)=\sqrt{\frac{2}{L+1}}\sum_{m=0}^{L-1}f\left(m+1\right)\sin\left(\frac{\pi \left(n+1\right)\left(m+1\right)}{L+1}\right)
%
%   The implementation of this functions uses a simple algorithm that requires
%   an FFT of length $2N+2$, which might potentially be the product of a large
%   prime number. This may cause the function to sometimes execute slowly.
%   If guaranteed high speed is a concern, please consider using one of the
%   other DST transforms.
%
%   Examples:
%   ---------
%
%   The following figures show the first 4 basis functions of the DSTI of
%   length 20:::
%
%     % The dsti is its own adjoint.
%     F=dsti(eye(20));
%
%     for ii=1:4
%       subplot(4,1,ii);
%       stem(F(:,ii));
%     end;
%
%   See also:  dcti, dstiii, dstiv
%
%   References: rayi90 wi94

%   AUTHOR: Peter L. Søndergaard
%   TESTING: TEST_PUREFREQ
%   REFERENCE: REF_DSTI

complainif_argnonotinrange(nargin,1,3,mfilename);

if nargin<3
  dim=[];
end;

if nargin<2
  L=[];
end;

[f,L,Ls,W,dim,permutedsize,order]=assert_sigreshape_pre(f,L,dim,'DSTI');

if ~isempty(L)
  f=postpad(f,L);
end;

if L==1
  c=f;
 
else

  c = comp_dst(f,1);

end;

c=assert_sigreshape_post(c,dim,permutedsize,order);

