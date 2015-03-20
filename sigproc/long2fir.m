function g=long2fir(g,varargin);
%LONG2FIR   Cut LONG window to FIR
%   Usage:  g=long2fir(g,L);
%
%   `long2fir(g,L)` will cut the LONG window *g* to a length *L* FIR window by
%   cutting out the middle part. Note that this is a slightly different
%   behaviour than |middlepad|.
%
%   `long2fir(g,L,'wp')` or `long2fir(g,L,'hp')` does the same assuming the
%   input window is a whole-point even or half-point even window,
%   respectively.
%
%   See also:  fir2long, middlepad

if nargin<1
  error('%s: Too few input parameters.',upper(mfilename));
end;

definput.flags.centering = {'unsymmetric','wp','hp'};
definput.keyvals.L      = [];
definput.keyvals.cutrel = [];

[flags,kv,L]=ltfatarghelper({'L'},definput,varargin);

W=length(g);

if W<L
  error('L must be smaller than length of window.');
end;

if ~isempty(kv.cutrel)
  maxval=max(abs(g));
  mask=abs(g)>maxval*kv.cutrel;
  L=W-2*min(abs(find(mask)-L/2));
end;

if isempty(L)
    error(['%s: You must specify a way to shorten the window, either by ' ...
           'specifying the length or through a flag.'],upper(mfilename));
end;

if flags.do_unsymmetric
  % No assumption on the symmetry of the window.

  if rem(L,2)==0
    % HPE middlepad works the same way as the FIR cutting (e.g. just
    % removing middle points) for even values of L.
    g=middlepad(g,L,'hp');
  else
    % WPE middlepad works the same way as the FIR cutting (e.g. just
    % removing middle points) for odd values of L.
    g=middlepad(g,L);
  end;
  
else
  if flags.do_wp
    g=middlepad(g,L);
    if rem(L,2)==0
      g(L/2+1)=0;
    end;
  else
    g=middlepad(g,L,'hp');
    if rem(L,2)==1
      g(ceil(L/2))=0;
    end;
  end;
end;

