function [c,newphase,usedmask,tgrad,fgrad]=filterbankconstphase(s,a,tfr,fc,varargin)
%% FILTERBANKCONSTPHASE Construct phase for FILTERBANK and UFILTERBANK
%   Usage:  c=filterbankconstphase(s,a,tfr,fc);
%           [c,newphase,usedmask,tgrad,fgrad] = filterbankconstphase(...);
%
%   Input parameters:
%         s        : Initial coefficients.
%         a        : Downsampling factor(s).
%         tfr      : Time-frequency rations of the filters
%         fc       : Center frequencies of the filters
%         mask     : Mask for selecting known phase.
%         usephase : Explicit known phase.
%   Output parameters:
%         c        : Coefficients with the constructed phase.
%         newphase : Just the (unwrapped) phase.
%         usedmask : Mask for selecting coefficients with the new phase.
%         tgrad    : Relative time phase derivative.
%         fgrad    : Relative frequency phase derivative.
%
%   `filterbankconstphase(s,a,tfr,fc)` will construct a suitable phase for 
%   the positive valued coefficients `abs(s)` of a filterbank with center
%   frequencies *fc* and time-frequency ratios *tfr* and subsampling
%   factors *a*. It is assumed that the coefficients are
%
%   If *s* contains the absolute values of the Gabor coefficients of a signal
%   obtained using the window *g*, time-shift *a* and number of channels 
%   *M*, i.e.:
%
%     c=dgtreal(f,g,a,M);
%     s=abs(c);
%
%   `filterbankconstphase(s,g,a,M,tol)` does as above, but sets the phase of
%   coefficients less than *tol* to random values.
%   By default, *tol* has the value 1e-10. 
%
%   `filterbankconstphase(c,g,a,M,tol,mask)` accepts real or complex valued
%   *c* and real valued *mask* of the same size. Values in *mask* which can
%   be converted to logical true (anything other than 0) determine
%   coefficients with known phase which is used in the output. Only the
%   phase of remaining coefficients (for which mask==0) is computed.
%
%   `filterbankconstphase(c,g,a,M,tol,mask,usephase)` does the same as before
%   but uses the known phase values from *usephase* rather than from *c*.
%
%   The function accepts the following additional paramaters:
%
%       'tol',tol 
%   In addition, *tol* can be a vector containing decreasing values. In 
%   that case, the algorithm is run `numel(tol)` times, initialized with
%   the result from the previous step in the 2nd and the further steps.
%
%   This function requires a computational subroutine that is only
%   available in C. Use |ltfatmex| to compile it.
%
%   See also:  ltfatmex
%
%   References: ltfatnote051
%

% AUTHOR: Nicki Holighaus, Zdenek Prusa

thismfilename = upper(mfilename);
complainif_notenoughargs(nargin,4,thismfilename);

if ~(isnumeric(s) || iscell(s)) || isempty(s)
    error('%s: *s* must be numeric or cell.',thismfilename);
 end

if ~isnumeric(a) || isempty(a)
    error('%s: a must be non-empty numeric.',upper(mfilename));
end

definput.keyvals.tol=[1e-1,1e-10];
definput.keyvals.gderivweight=1/2;
definput.keyvals.mask=[];
definput.keyvals.usephase=[];
definput.flags.real = {'real','complex'};
[flags,kv,mask,usephase]=ltfatarghelper({'mask','usephase'},definput,varargin);
tol = kv.tol;
  
if ~isnumeric(tol) || ~isequal(tol,sort(tol,'descend'))
    error(['%s: *tol* must be a scalar or a vector sorted in a ',...
          'descending manner.'],thismfilename);
end

if ~isempty(usephase) && isempty(mask)
    error('%s: Both mask and usephase must be used at the same time.',...
          upper(mfilename));
end

if ~isempty(usephase)
    complainif_notequalsize('s',s,'usephase',usephase,thismfilename);
end

if ~isempty(mask)
    complainif_notequalsize('s',s,'mask',mask,thismfilename);
end

do_uniform = 1;
wasCell = 0;

if iscell(s)
    M = numel(s);
    N = cellfun(@(sEl) size(sEl,1),s);
    W = size(s{1},2);

    asan = comp_filterbank_a(a,M);
    a = asan(:,1)./asan(:,2);
    
    wasCell = 1;
    
    if any( N ~= N(1)) && any( a ~= a(1) )
        do_uniform = 0;
        s = cell2mat(s);
        if ~isempty(mask)
            mask = cell2mat(s);
        end
        if ~isempty(usephase)
            usephase = cell2mat(usephase);
        end
    else
        a = a(1);
        smat = zeros(N(1),M,W);
        for m=1:M    
            smat(:,m,:)=(s{m});
        end
        s = smat;
    end
else
    [N,M,W] = size(s);
    
    asan = comp_filterbank_a(a,M);
    a = asan(:,1)./asan(:,2);
    a = a(1);
end

abss = abs(s);

if isempty(usephase)
    usephase = angle(s);
else
    if ~isreal(usephase)
        error('%s: usephase must be real.',thismfilename);
    end
end

if ~isempty(mask) 
    if ~isreal(mask)
        error('%s: mask must be real.',thismfilename);
    end
    mask = cast(mask,'double');
    mask(mask~=0) = 1;
end

if ~isnumeric(tfr) || isempty(tfr) || numel(tfr) ~= M
  error('%s: tfr must be non-empty numeric.',upper(mfilename));
end;

if ~isnumeric(fc) || isempty(fc) || numel(fc) ~= M
  error('%s: fc must be non-empty numeric.',upper(mfilename));
end

if do_uniform
    [tgrad,fgrad,logs] = ...
        comp_ufilterbankphasegradfrommag(abss,N(1),a,M,sqrt(tfr),fc,flags.do_real);    
    [newphase,usedmask] = ...
        comp_ufilterbankconstphase(abss,tgrad,fgrad,fc,mask,usephase,a,tol,flags.do_real);
else
    NEIGH = comp_filterbankneighbors(a,M,N,flags.do_real);
    chanStart = [0;cumsum(N)];

    posInfo = zeros(chanStart(end),2);
    for kk = 1:M
        posInfo(chanStart(kk)+(1:N(kk)),:) = [(kk-1)*ones(N(kk),1),(0:N(kk)-1)'.*a(kk)];
    end
    posInfo = posInfo.';  

    NEIGH = NEIGH-1;
    [tgrad,fgrad,logs] = ...
        comp_filterbankphasegradfrommag(abss,N,a,M,sqrt(tfr),fc,NEIGH,posInfo,kv.gderivweight);
    [newphase,usedmask] = ...
        comp_filterbankconstphase(abss,tgrad,fgrad,NEIGH,posInfo,fc,mask,usephase,a,M,N,tol);
end

c = abss.*exp(1i*newphase);
if wasCell
    % Apply the phase and convert back to cell array 
    c = mat2cell(c(:),N,W);
    newphase = mat2cell(newphase(:),N,W);
    usedmask = mat2cell(usedmask(:),N,W);
    tgrad = mat2cell(tgrad(:),N,W);
    fgrad = mat2cell(fgrad(:),N,W);
end

function complainif_notequalsize(aname,a,bname,b,thismfilename)

if ~isequal(class(a), class(b))
    error('%s: %s and %s must be of the same type.',...
          aname, bname, thismfilename);
end

if ~isequal(size(a), size(b))
    error('%s: %s and %s must have equal sizes.',...
          aname, bname,thismfilename);
end

if iscell(a)
    if ~all(cellfun(@(aEl,bEl) isequal(aEl,bEl),a,b))
        error('%s: %s and %s must have equal sizes.',...
              aname, bname,thismfilename);
    end
end




