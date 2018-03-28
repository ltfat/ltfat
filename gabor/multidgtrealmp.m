function [c, frec, info] = multidgtrealmp(f,dicts,varargin)
%FRANAMP  Frame Analysis by Matching Pursuit
%   Usage:  c = multidgtrealmp(f,dicts)
%           c = multidgtrealmp(f,dicts,errdb,maxit)
%           [c,frec,info] = multidgtrealmp(...)
%
%   Input parameters:
%       f        : Input signal
%       dicts    : Dictionaries. Format {g1,a1,M1,g2,a2,M2,...}
%       errdb    : Target normalized approximation error in dB
%       maxit    : Maximum number of iterations.
%   Output parameters:
%       c        : Sparse representation
%       frec     : Reconstructed signal
%       info     : Struct with additional output paramerets
%
%   `multidgtrealmp(f,{g1,a1,M1,g2,a2,M2,...,gW,aW,MW})` returns sparse 
%   representation of a signal in `W` Gabor dictionaries defined using the 
%   fast matching pursuit algorithm. `gw` is a Gabor window defined as in
%   as in |dgtreal|, `aw` is a hop factor, `Mw` is the number of frequency
%   channels. All `aw` and `Mw` must be divisible by `min(a1,...,aW)`.
%
%   `multidgtrealmp(f,dicts,errdb,maxit)` tries to reach normalized 
%   approximation error *errdb* dB in at most *maxit* iterations. 
%
%   `[c,frec,info] = multidgtrealmp(...)` in addition returns the aproximated
%   signal *frec* and a struct `info` with the following fields:
%
%     .iter    Number of iterations done.
%
%     .relres  Vector of length `.iter` with approximation error progress. 
%
%   The normalized approximation error is computed as 
%   `err=norm(f-frec)/norm(f)`. The relationship between the output 
%   coefficients and the approximation is 
%   `frec = sum(frsyn(F,c))`.
%
%   The function takes the following optional parameters at the end of
%   the line of input arguments:
%
%     'print'    Display the progress.
%
%     'printstep',p
%                If 'print' is specified, then print every p'th
%                iteration. Default value is 10;
%
%   Algorithms
%   ----------
%
%     'fast'    Fast MP using approximate update in the coefficient domain.
%               This is the default. 
%               Requires MEX files to be compiled.
%
%     'slow'    By-the-book implementation
%
%   Examples
%   --------
%
%   The following example show the development of the approx. error for the
%   MP and OMP algorithms. :::
%
%
%   See also: dgtreal
%
%   References: mazh93 ltfatnote052

%AUTHOR: Zdenek Prusa

thismfile = upper(mfilename);
complainif_notenoughargs(nargin,2,thismfile);
    
% Define initial value for flags and key/value pairs.
definput.keyvals.errdb=-40;
definput.keyvals.maxit=[];
definput.keyvals.printstep=100;
definput.keyvals.kernthr = 1e-4;
definput.keyvals.dictmask = [];
definput.flags.print={'quiet','print'};
definput.flags.algorithm={'fast','slow'};
definput.flags.search={'plainsearch','pedanticsearch'};
definput.flags.phaseconv={'freqinv','timeinv'};
[flags,kv]=ltfatarghelper({'errdb','maxit'},definput,varargin);

%% ----- step 1 : Verify f and determine its length -------
% Change f to correct shape.
[f,~,Ls,W,dim,permutedsize,order]=assert_sigreshape_pre(f,[],[],upper(mfilename));

if W>1
    error('%s: Input signal can be single channel only.',upper(mfilename));
end

if kv.errdb > 0
    error('%s: Target error must be lower than 0 dB.',upper(mfilename));
end

if ~iscell(dicts), error('%s: dicts must be cell',thismfile); end
if rem(numel(dicts),3) ~= 0 || ~all(cellfun(@(x)isscalar(x), dicts([2:3:end,3:3:end])))
    error('%s: bad format of dicts. Check {g1,a1,M1,g2,a2,M2,...,gW,aW,MW}',...
        thismfile); 
end

dictno = numel(dicts)/3;

if isempty(kv.dictmask) 
    kv.dictmask = ones(dictno,1);
else
    if ~isvector(kv.dictmask) || numel(kv.dictmask) ~= dictno
        error('%s: dictmask length must be equal to the number of dictionaries.',thismfile);
    end
end

gin = dicts(1:3:end);
a = cell2mat(dicts(2:3:end));
M = cell2mat(dicts(3:3:end));

amin = min(a);
Mmax = max(M);

if any(rem([a(:);M(:)],amin) ~= 0)
    error('%s: all a and M must be divisible by min(a1,...,aW)',thismfile);
end

L = dgtlength(Ls,amin,Mmax);

g = cell(dictno,1);
for dIdx = 1:dictno
    g{dIdx} = normalize(gabwin(gin{dIdx},a(dIdx),M(dIdx),L),'2');
end

if isempty(kv.maxit), kv.maxit = L; end
info.iter = 0;
info.relres = [];

fpad = postpad(f,L);
fnorm = norm(fpad);

c = comp_multidgtrealmp(fpad,g,a,M,flags.do_timeinv,kv.kernthr,kv.errdb,kv.maxit,kv.maxit,flags.do_pedanticsearch);

if nargout>1
  frec = zeros(size(fpad,1),dictno,class(fpad));

  for dIdx = 1:dictno
        frec(:,dIdx) = idgtreal(c{dIdx},g{dIdx},a(dIdx),M(dIdx),flags.phaseconv); 
  end
  frec = frec(1:Ls,:);

  permutedsize2 = permutedsize; permutedsize2(2) = dictno;
  info.frecdict = assert_sigreshape_post(frec,dim,permutedsize2,order);
  frec = sum(frec,2);
  info.relres = norm(frec-f)/fnorm;
  frec = assert_sigreshape_post(frec,dim,permutedsize,order);
end
