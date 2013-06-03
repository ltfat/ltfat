function [c,Ls,g]=dwilt(f,g,M,varargin)
%DWILT  Discrete Wilson transform
%   Usage:  c=dwilt(f,g,M);
%           c=dwilt(f,g,M,L);
%           [c,Ls]=dwilt(...);
%
%   Input parameters:
%         f     : Input data
%         g     : Window function.
%         M     : Number of bands.
%         L     : Length of transform to do.
%   Output parameters:
%         c     : $2M \times N$ array of coefficients.
%         Ls    : Length of input signal.
%
%   `dwilt(f,g,M)` computes a discrete Wilson transform with *M* bands and
%   window *g*.
%
%   The length of the transform will be the smallest possible that is
%   larger than the signal. *f* will be zero-extended to the length of the 
%   transform. If *f* is a matrix, the transformation is applied to each column.
%
%   The window *g* may be a vector of numerical values, a text string or a
%   cell array. See the help of |wilwin| for more details.
%
%   `dwilt(f,g,M,L)` computes the Wilson transform as above, but does a
%   transform of length *L*. *f* will be cut or zero-extended to length *L*
%   before the transform is done.
%
%   `[c,Ls]=dwilt(f,g,M)` or `[c,Ls]=dwilt(f,g,M,L)` additionally return the
%   length of the input signal *f*. This is handy for reconstruction::
%
%     [c,Ls]=dwilt(f,g,M);
%     fr=idwilt(c,gd,M,Ls);
%
%   will reconstruct the signal *f* no matter what the length of *f* is, provided
%   that *gd* is a dual Wilson window of *g*.
%
%   `[c,Ls,g]=dwilt(...)` additionally outputs the window used in the
%   transform. This is useful if the window was generated from a description
%   in a string or cell array.
%
%   A Wilson transform is also known as a maximally decimated, even-stacked
%   cosine modulated filter bank.
%
%   Use the function |wil2rect| to visualize the coefficients or to work
%   with the coefficients in the TF-plane.
%
%   Assume that the following code has been executed for a column vector *f*::
%
%     c=dwilt(f,g,M);  % Compute a Wilson transform of f.
%     N=size(c,2)*2;   % Number of translation coefficients.
%
%   The following holds for $m=0,\ldots,M-1$ and $n=0,\ldots,N/2-1$:
%
%   If $m=0$:
%
%   ..               L-1 
%     c(m+1,n+1)   = sum f(l+1)*g(l-2*n*M+1)
%                    l=0  
%
%   .. math:: c\left(1,n+1\right) = \sum_{l=0}^{L-1}f(l+1)g\left(l-2nM+1\right)
%
%
%   If $m$ is odd and less than $M$
%
%   ..               L-1 
%     c(m+1,n+1)   = sum f(l+1)*sqrt(2)*sin(pi*m/M*l)*g(k-2*n*M+1)
%                    l=0  
% 
%                    L-1 
%     c(m+M+1,n+1) = sum f(l+1)*sqrt(2)*cos(pi*m/M*l)*g(k-(2*n+1)*M+1)
%                    l=0  
%
%   .. math:: c\left(m+1,n+1\right) & = & \sqrt{2}\sum_{l=0}^{L-1}f(l+1)\sin(\pi\frac{m}{M}l)g(l-2nM+1)\\
%             c\left(m+M+1,n+1\right) & = &
%             \sqrt{2}\sum_{l=0}^{L-1}f(l+1)\cos(\pi\frac{m}{M}l)g\left(l-\left(2n+1\right)M+1\right)
%
%   If $m$ is even and less than $M$
%
%   ..               L-1 
%     c(m+1,n+1)   = sum f(l+1)*sqrt(2)*cos(pi*m/M*l)*g(l-2*n*M+1)
%                    l=0  
% 
%                    L-1 
%     c(m+M+1,n+1) = sum f(l+1)*sqrt(2)*sin(pi*m/M*l)*g(l-(2*n+1)*M+1)
%                    l=0  
%
%   .. math:: c\left(m+1,n+1\right) & = & \sqrt{2}\sum_{l=0}^{L-1}f(l+1)\cos(\pi\frac{m}{M}l)g(l-2nM+1)\\
%             c\left(m+M+1,n+1\right) & = &
%             \sqrt{2}\sum_{l=0}^{L-1}f(l+1)\sin(\pi\frac{m}{M}l)g\left(l-\left(2n+1\right)M+1\right)
%
%   if $m=M$ and $M$ is even:
%
%   ..               L-1 
%     c(m+1,n+1)   = sum f(l+1)*(-1)^(l)*g(l-2*n*M+1)
%                    l=0
%
%   .. math:: c\left(M+1,n+1\right) = \sum_{l=0}^{L-1}f(l+1)(-1)^{l}g(l-2nM+1)
%
%   else if $m=M$ and $M$ is odd
%
%   ..               L-1 
%     c(m+1,n+1)   = sum f(l+1)*(-1)^l*g(l-(2*n+1)*M+1)
%                    l=0
%
%   .. math:: c\left(M+1,n+1\right) = \sum_{k=0}^{L-1}f(l+1)(-1)^{l}g\left(l-\left(2n+1\right)M+1\right)
%
%   See also:  idwilt, wilwin, wil2rect, dgt, wmdct, wilorth
%
%   References: bofegrhl96-1 liva95 dajajo91

%   AUTHOR : Peter L. Søndergaard.
%   TESTING: TEST_DWILT
%   REFERENCE: REF_DWILT

if nargin<3
  error('%s: Too few input parameters.',upper(mfilename));
end;

definput.keyvals.L=[];
definput.keyvals.dim=[];
[flags,kv,L]=ltfatarghelper({'L'},definput,varargin);


%% ----- step 1 : Verify f and determine its length -------
% Change f to correct shape.
[f,dummy,Ls,W,dim,permutedsize,order]=assert_sigreshape_pre(f,[],kv.dim,upper(mfilename));

%% ------ step 2: Verify a, M and L
if isempty(L)

    % ----- step 2b : Verify a, M and get L from the signal length f----------
    L=dwiltlength(Ls,M);

else

    % ----- step 2a : Verify a, M and get L
    Luser=dwiltlength(L,M);
    if Luser~=L
        error(['%s: Incorrect transform length L=%i specified. Next valid length ' ...
               'is L=%i. See the help of DWILTLENGTH for the requirements.'],...
              upper(mfilename),L,Luser);
    end;

end;

%% ----- step 3 : Determine the window 

[g,info]=wilwin(g,M,L,upper(mfilename));

if L<info.gl
  error('%s: Window is too long.',upper(mfilename));
end;

%% ----- step 4: final cleanup ---------------

f=postpad(f,L);

% If the signal is single precision, make the window single precision as
% well to avoid mismatches.
if isa(f,'single')
  g=single(g);
end;

%% ----- Call the computational subroutines.
c=comp_dwilt(f,g,M);

%% ----- reorder coefficients to correct final layout
order=assert_groworder(order);
permutedsize=[2*M,L/(2*M),permutedsize(2:end)];

c=assert_sigreshape_post(c,dim,permutedsize,order);
