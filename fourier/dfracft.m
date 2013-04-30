function frf=dfracft(f,a,varargin)
%DFRACFT  Discrete Fractional Fourier transform
%   Usage:  V=dfracft(f,a,p);
%           V=dfracft(f,a);
%
%   `dfracft(f,a)` computes the discrete fractional Fourier Transform of the
%   signal *f* to the power *a*. For *a=1* it corresponds to the ordinary
%   discrete Fourier Transform. If *f* is multi-dimensional, the
%   transformation is applied along the first non-singleton dimension.
%
%   `dfracft(f,a,dim)` does the same along dimension *dim*.   
%
%   `dfracft(f,a,[],p)` or `dfracft(f,a,dim,p)` allows to choose the order
%   of approximation of the second difference operator (default: *p=2*).
%
%   See also:  ffracft, dft, hermbasis, pherm
%
%   References: ozzaku01,buma04

%   AUTHOR : Christoph Wiesmeyr 
%   TESTING: TEST_HERMBASIS
%   REFERENCE: OK

if nargin<2
  error('%s: Too few input parameters.',upper(mfilename));
end;

definput.keyvals.p   = 2;
definput.keyvals.dim = [];
[flags,keyvals,dim,p]=ltfatarghelper({'dim','p'},definput,varargin);

[f,L,Ls,W,dim,permutedsize,order]=assert_sigreshape_pre(f,[],dim,upper(mfilename));

H = hermbasis(L,p);

% set up the eigenvalues
k=0:L-1;
lam = exp(-1i*k*a*pi/2);
lam=lam(:);

% correction for even signal lengths
if ~rem(L,2)
    lam(end)=exp(-1i*L*a*pi/2);
end

% shuffle the eigenvalues in the right order
even=~mod(L,2);
cor=2*floor(L/4)+1;
for k=(cor+1):2:(L-even)
    lam([k,k+1])=lam([k+1,k]);
end

frf =H*(bsxfun(@times,lam,H'*f));

frf=assert_sigreshape_post(frf,dim,permutedsize,order);



