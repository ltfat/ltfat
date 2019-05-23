function [c, frec, info] = multidgtrealmp(f,dicts,varargin)
%MULTIDGTREALMP  Matching Pursuit Decomposition with Multi-Gabor Dictionary
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
%   representation of a signal in `W` Gabor dictionaries using the 
%   fast matching pursuit algorithm. `gw` is a Gabor window defined
%   as in |dgt| and |dgtreal|, `aw` is a hop factor, `Mw` is the number of 
%   frequency channels. All `aw` and `Mw` must be divisible by `min(a1,...,aW)`.
%   The algorithm will try to reach -40 dB relative approximation error
%   in at most `numel(f)` iterations.
%   The function returns a cell-array with elements storing coefficients
%   for individual Gabor systems such that they can be directly used in
%   |idgtreal|. 
%
%   `multidgtrealmp(f,dicts,errdb,maxit)` tries to reach normalized 
%   approximation error *errdb* dB in at most *maxit* iterations.
%
%   `[c,frec,info] = multidgtrealmp(...)` in addition returns the
%   aproximation *frec* and a struct `info` with the following fields:
%
%     .iter     Number of iterations done.
%
%     .atoms    Number of atoms selected.
%
%     .relres   Final approximation error. 
%
%     .g        Cell array of numeric windows used in the multi-dictionary
%
%     .a        Array of hop factors for indivitual dictionaries
%
%     .M        Array of numbers of channels for individual dictionaries
%
%     .synthetize  Anonymous function which can be used to synthetize from
%                  the (modified) coefficients as 
%                  `frec = sum(info.synthetize(c),dim)`
%                  where `dim=2` if the input *f* was a column vector and
%                  `dim=1` if it was a row vector. 
%
%   The normalized approximation error is computed as 
%   `err=norm(f-frec)/norm(f)`.
%
%   The function takes the following optional parameters at the end of
%   the line of input arguments:
%
%     'kenrnthr',kt    Kernel threshold. Must be in range ]0,1]. Default is 1e-4.
%
%     'timeinv'        Use the time invariant phase convention. The
%                      default is 'freqinv'.
%
%     'pedanticsearch' Be pedantic about the energy of pairs of conjugated
%                      atoms in the selection step. Disbaled by default.
%
%     'algorithm',alg  Algorithm to use. Available: 
%                      'mp'(default),'selfprojmp','cyclicmp'
%
%     'resetit',it     Reset the decomposition every `it` iteration.
%                      The reset means a complete synthesis and re-analysis.
%                      The default is 0 (no reset)
%
%   The computational routine is only available in C. Use |ltfatmex| to
%   to compile it.
%
%   Algorithms
%   ----------
%
%   By default, the function uses the fast MP using approximate update 
%   in the coefficient domain as described in::
%   
%   "Z. Prusa, N. Holighaus, P. Balazs: Fast Matching Pursuit with Multi-Gabor Dictionaries"
%   
%   The kernel threshold limits the minimum approximation error which
%   can be reached. For example the default threshold 1e-4 in general 
%   allows achieving at least -40 dB.
%   
%   Examples
%   --------
%
%   The following example shows the decomposition in 3 dictionaries and
%   plots contributions from the individual dictionaries and the residual.:::
%
%       [f,fs] = gspi;
%       [c, frec, info] = multidgtrealmp(f,...
%       {'blackman',128,512,'blackman',512,2048,'blackman',2048,8192});
%       frecd = info.synthetize(c);
%       figure(1); 
%       xvals = (0:numel(f)-1)/fs;
%       subplot(4,1,1); plot(xvals,frecd(:,1));ylim([-0.5,0.5]);
%       subplot(4,1,2); plot(xvals,frecd(:,2));ylim([-0.5,0.5]);
%       subplot(4,1,3); plot(xvals,frecd(:,3));ylim([-0.5,0.5]);
%       subplot(4,1,4); plot(xvals,f-frec);ylim([-0.5,0.5]); xlabel('Time (s)');
%
%   See also: dgtreal idgtreal
%
%   References: ltfatnote052 mazh93 stchr10 rerose17

%AUTHOR: Zdenek Prusa

thismfile = upper(mfilename);
complainif_notenoughargs(nargin,2,thismfile);
% Define initial value for flags and key/value pairs.
definput.keyvals.errdb=-40;
definput.keyvals.maxit=[];
definput.keyvals.resetit=0;
%definput.keyvals.iterstep=[];
definput.keyvals.kernthr = 1e-4;
%definput.flags.print={'quiet','print'};
definput.flags.algversion={'fast','slow'};
definput.flags.algorithm={'mp','selfprojmp','cyclicmp'};
definput.flags.search={'plainsearch','pedanticsearch'};
definput.flags.phaseconv={'freqinv','timeinv'};
[flags,kv]=ltfatarghelper({'errdb','maxit'},definput,varargin);

if exist('comp_multidgtrealmp','file') ~= 3 && flags.do_fast
    error(['%s: MEX/OCT file is missing. Either compile the MEX/OCT ',...
           'interfaces or re-run the function with ''slow'''], thismfile);
end

if flags.do_slow
    error('%s: ''slow'' is not supported yet.',thismfile)
end

%% ----- step 1 : Verify f and determine its length -------
% Change f to correct shape.
[f,~,Ls,W,dim,permutedsize,order]=assert_sigreshape_pre(f,[],[],upper(mfilename));

if W>1
    error('%s: Input signal can be single channel only.',upper(mfilename));
end

if kv.errdb > 0
    error('%s: Target error must be lower than 0 dB.',upper(mfilename));
end

if ~(kv.kernthr > 0 && kv.kernthr <= 1)
    error('%s: Kenel threshold must be in range ]0,1].',upper(mfilename));
end

if ~iscell(dicts), error('%s: dicts must be cell',thismfile); end
if rem(numel(dicts),3) ~= 0 || ~all(cellfun(@(x)isscalar(x), dicts([2:3:end,3:3:end])))
    error('%s: bad format of dicts. Check {g1,a1,M1,g2,a2,M2,...,gW,aW,MW}',...
        thismfile);
end

dictno = numel(dicts)/3;
gin = dicts(1:3:end);
a = cell2mat(dicts(2:3:end));
M = cell2mat(dicts(3:3:end));

if any(rem(M,a) ~= 0) || any(M./a<2)
    error(['%s: Only integer oversampling greater than 1 is allowed ',...
           'i.e. M/a must be an integer>=2.'],...
    upper(mfilename));
end

if dictno > 1
    asort  = sort(a);
    Msort  = sort(M);
    if any(rem(asort(2:end),asort(1:end-1)) ~= 0)
        error('%s: all au and av must be divisible by min(au,av)',thismfile);
    end
    if any(rem(Msort(2:end),Msort(1:end-1)) ~= 0)
        error('%s: all Mu and Mv must be divisible by min(Mu,Mv)',thismfile);
    end
end

info.a = a;
info.M = M;
info.iter = 0;
info.relres = [];
fnorm = norm(f);

L = filterbanklength(Ls,[a(:);M(:)]);
if isempty(kv.maxit), kv.maxit = L; end

info.g = cell(dictno,1);
for dIdx = 1:dictno
    info.g{dIdx} = normalize(gabwin(gin{dIdx},a(dIdx),M(dIdx),L),'2');
end

for dIdx = 1:dictno
    condnum = gabframebounds(info.g{dIdx},a(dIdx),M(dIdx));
    if condnum > 1e3
        error('%s: Dictionary %d is badly conditioned.',dIdx,upper(mfilename));
    end
end

% Initial residuum
fpad = postpad(f,L);

if kv.resetit == 0
    [c,info.atoms,info.iter,info.status] = ...
        comp_multidgtrealmp(fpad,info.g,a,M,flags.do_timeinv,...
                            kv.kernthr,kv.errdb,kv.maxit,kv.maxit,...
                            flags.do_pedanticsearch, flags.algorithm );
else
    c = cell(1,dictno);
    for dIdx = 1:dictno
        c{dIdx} = zeros(floor(M(dIdx)/2) + 1, L/a(dIdx),class(f));
    end

    remainingit = kv.maxit;
    currerrdb = kv.errdb;

    while remainingit > 0 
        currit = min([remainingit,kv.resetit]);

        [ctmp,info.atoms,info.iter,info.status] = ...
            comp_multidgtrealmp(fpad,info.g,a,M,flags.do_timeinv,...
                                kv.kernthr,currerrdb,currit,currit,...
                                flags.do_pedanticsearch, flags.algorithm );

        c = cellfun(@(cEl,ctmpEl) cEl + ctmpEl, c(:)',ctmp(:)','UniformOutput', 0);

        remainingit = remainingit - info.iter;

        %The following means that some other stopping criterion has been reached
        if currit ~= info.iter, break; end

        if remainingit > 0
            % Recompute residuum
            frec = sum(cell2mat(cellfun(@(cEl,gEl,aEl,MEl) idgtreal(cEl,gEl,aEl,MEl,flags.phaseconv),...
                   c(:)',info.g(:)',num2cell(a(:))',num2cell(M(:))','UniformOutput',0)),2);
            
            fpad = postpad(f,L) - frec;
            % Since we are changing the residual, we must adjust the target error
            currerrdb = kv.errdb - 20*log10( norm(fpad(:)) / fnorm );
        end
    end

    %Fix the info
    info.atoms = sum(cellfun(@(cEl) numel(find( abs(cEl) > 0 )), c));
    info.iter = kv.maxit - remainingit;
end


if nargout>1
  permutedsize2 = permutedsize; permutedsize2(2) = dictno;
  info.synthetize = @(c) ...
      assert_sigreshape_post(...
      postpad(cell2mat(cellfun(@(cEl,gEl,aEl,MEl) idgtreal(cEl,gEl,aEl,MEl,flags.phaseconv),...
      c(:)',info.g(:)',num2cell(a(:))',num2cell(M(:))','UniformOutput',0)),Ls),...
      dim,permutedsize2,order);
  dim2 = 2;
  if dim == 2, dim2 = 1; end
  frec = sum(info.synthetize(c),dim2);
  if fnorm == 0
    info.relres = 0;
  else
    info.relres = norm(frec(:)-f(:))/fnorm;
  end

  status_str = {...
  'Target error reached',...
  'Maximum number of atoms reached',...
  'Maximum number of iterations reached',...
  'Stalled (abs. norm. error estimate became negative)',...
  'Selected coefficient tolerance reached',...
  'All zeros'...
  };
  info.message = status_str{1 + info.status};
end
