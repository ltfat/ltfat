function f = ifwt(c,g,J,varargin)
%IFWT   Inverse Fast Wavelet Transform 
%   Usage:  f = ifwt(c,g,J,)
%           f = ifwt(c,g,J,Ls)
%
%   Input parameters:
%         c     : Coefficients stored in J+1 cell-array.
%         g     : Synthesis wavelet filters.
%         J     : Number of filterbank iterations.
%         Ls    : Length of signal.
%   Output parameters:
%         f     : Reconstructed data.
%
%
%   See also:
%
%   Demos:
%
%   References:

%   AUTHOR : Zdenek Prusa.
%   TESTING: TEST_IFWT
%   REFERENCE: REF_IFWT


if nargin<3
  error('%s: Too few input parameters.',upper(mfilename));
end;

% works with coefficients obtained from input signals of lengts a*2^J
definput.keyvals.Ls=[];
definput.flags.type = {'dec','undec'};
[flags,kv,Ls]=ltfatarghelper({'Ls'},definput,varargin);

if ~isnumeric(J) || ~isscalar(J)
  error('%s: "J" must be a scalar.',upper(mfilename));
end;

if(J<1)
   error('%s: J must be a positive integer.',upper(callfun)); 
end

if(~iscell(g) || length(g)<2)
   error('%s: h is expected to be a cell array containing two wavelet filters.',upper(callfun)); 
end

[sigHalfLen,W] = size(c{end});

if(flags.do_dec)
  f = comp_ifwt(c,g,J);
elseif(flags.do_undec)
  f = comp_ifwt_undec(c,g,J);  
end


% Cut or extend f to the correct length, if desired.
if ~isempty(Ls)
  f=postpad(f,Ls);
end;

f=comp_sigreshape_post(f,Ls,0,[0; W]);