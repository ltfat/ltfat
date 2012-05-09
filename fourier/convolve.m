function h=convolve(f,g,varargin)
%CONVOLVE  Convolution
%   Usage:  h=convolve(f,g);
%
%   `convolve(f,g)` will convolve two vectors *f* and *g*. The
%   convolution is not periodic, so the output will have a length of ::
%
%     output_len = length(f)+length(g)-1;
%
%   This function works entirely similar to the Matlab routine `conv`
%   using instead a fast FFT algorithm, making it much faster if one or
%   more of the signals are long.
%
%   See also: pconv
  
if nargin<2
  error('%s: Too few input parameters.',upper(mfilename));
end;

definput.keyvals.L=[];
definput.keyvals.dim=[];

[flags,kv,L,dim]=ltfatarghelper({'L','dim'},definput,varargin);

[f,L1,Lf,Wf,dimout,permutedsize_f,order_f]=assert_sigreshape_pre(f,L,dim,'CONVOLVE');
[g,L2,Lg,Wg,dimout,permutedsize_g,order_g]=assert_sigreshape_pre(g,L,dim,'CONVOLVE');

Lh=Lf+Lg-1;

if (Wf>1) && (Wg>1)
  error('%s: Only one of the inputs can be multi-dimensional.',upper(mfilename));
end;

W=max(Wf,Wg);
if Wf<W
  f=repmat(f,1,W);
end;

if Wg<W
  g=repmat(g,1,W);
end;

if isreal(f) && isreal(g)
  h=comp_ifftreal(comp_fftreal(postpad(f,Lh)).*...
                  comp_fftreal(postpad(g,Lh)),Lh);
else
  h=ifft(fft(postpad(f,Lh)).*...
         fft(postpad(g,Lh)));
  
end;
