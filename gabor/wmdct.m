function [c,Ls,g]=wmdct(f,g,M,varargin)
%WMDCT  Windowed MDCT transform
%   Usage:  c=wmdct(f,g,M);
%           c=wmdct(f,g,M,L);
%           [c,Ls]=wmdct(...);
%
%   Input parameters:
%         f     : Input data
%         g     : Window function.
%         M     : Number of bands.
%         L     : Length of transform to do.
%   Output parameters:
%         c     : $M \times N$ array of coefficients.
%         Ls    : Length of input signal.
%
%   `wmdct(f,g,M)` computes a Windowed Modified Discrete Cosine Transform with
%   *M* bands and window *g*.
%
%   The length of the transform will be the smallest possible that is
%   larger than the signal. *f* will be zero-extended to the length of the 
%   transform. If *f* is a matrix, the transformation is applied to each column.
%   *g* must be whole-point even.
%
%   The window *g* may be a vector of numerical values, a text string or a
%   cell array. See the help of |wilwin| for more details.
%
%   `wmdct(f,g,M,L)` computes the MDCT transform as above, but does
%   a transform of length *L*. *f* will be cut or zero-extended to length *L*
%   before the transform is done.
%
%   `[c,Ls]=wmdct(f,g,M)` or [c,Ls]=wmdct(f,g,M,L)` additionally returns the
%   length of the input signal *f*. This is handy for reconstruction::
%
%     [c,Ls]=wmdct(f,g,M);
%     fr=iwmdct(c,gd,M,Ls);
%
%   will reconstruct the signal *f* no matter what the length of *f* is, provided
%   that *gd* is a dual Wilson window of *g*.
%
%   `[c,Ls,g]=wmdct(...)` additionally outputs the window used in the
%   transform. This is useful if the window was generated from a description
%   in a string or cell array.
%
%   The WMDCT is sometimes known as an odd-stacked cosine modulated filter
%   bank. The WMDCT defined by this routine is slightly different from the
%   most common definition of the WMDCT, in order to be able to use the
%   same window functions as the Wilson transform.
%
%   Assume that the following code has been executed for a column vector f
%   of length L::
%
%     c=wmdct(f,g,M);  % Compute the WMDCT of f.
%     N=size(c,2);     % Number of translation coefficients.
%
%   The following holds for $m=0,\ldots,M-1$ and $n=0,\ldots,N-1$:
%
%   If $m+n$ is even:
%
%   ..               L-1
%       c(m+1,n+1) = sum f(l+1)*cos(pi*(m+.5)*l/M+pi/4)*g(l-n*M+1)
%                    l=0
%
%   .. math:: c\left(m+1,n+1\right) = \sqrt{2}\sum_{l=0}^{L-1}f(l+1)\cos\left(\frac{\pi}{M}\left(m+\frac{1}{2}\right)l+\frac{\pi}{4}\right)g(l-nM+1)
%
%   If $m+n$ is odd:
%
%   ..              L-1
%      c(m+1,n+1) = sum f(l+1)*sin(pi*(m+.5)*l/M+pi/4)*g(l-n*M+1)
%                   l=0
%
%   .. math:: c\left(m+1,n+1\right) =
%      \sqrt{2}\sum_{l=0}^{L-1}f(l+1)\sin\left(\frac{\pi}{M}\left(m+\frac{1}{2}\right)l+\frac{\pi}{4}\right)g(l-nM+1)
%
%   Examples:
%   ---------
%
%   The following example shows the WMDCT coefficients (128 channels) of the
%   |greasy| test signal:::
%
%     fs=16000; % Sampling rate
%     c=wmdct(greasy,{'hann',0.02*fs},128);
%     plotwmdct(c,fs,90);
%
%   Compare the visual difference with the redundant expansion of the
%   same signal given in the example of the |dgtreal| function.
%
%   See also:  iwmdct, wilwin, dwilt, wildual, wilorth
%
%   References: prbr86 prjobr87 ma92 bohl96-1 

%   AUTHOR:    Peter L. Søndergaard
%   TESTING:   TEST_WMDCT
%   REFERENCE: REF_WMDCT

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
c  = comp_dwiltiii(f,g,M);

%% ----- reorder coefficients to correct final layout
order=assert_groworder(order);
permutedsize=[M,L/M,permutedsize(2:end)];

c=assert_sigreshape_post(c,dim,permutedsize,order);

