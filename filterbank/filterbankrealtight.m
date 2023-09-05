function gtout=filterbankrealtight(g,a,varargin)
%FILTERBANKREALTIGHT  Tight filters of filterbank for real signals only 
%   Usage:  gt=filterbankrealtight(g,a,L);
%           gt=filterbankrealtight(g,a);
%
%   `filterabankrealtight(g,a,L)` computes the canonical tight filters of 
%   *g* for a channel subsampling rate of *a* (hop-size) and a system 
%   length *L*. *L* must be compatible with subsampling rate *a* as 
%   `L==filterbanklength(L,a)`. The tight filters work only for real-valued
%   signals. Use this function on the common construction where the filters
%   in *g* only covers the positive frequencies.
%
%   `filterabankrealtight(g,a)` does the same, but the filters must be FIR
%   filters, as the transform length is unspecified. *L* will be set to 
%   next suitable length equal or bigger than the longest impulse response.  
%
%   The format of the filters *g* are described in the help of |filterbank|.
%
%   In addition, the function recognizes a 'forcepainless' flag which
%   forces treating the filterbank *g* and *a* as a painless case
%   filterbank.  
%
%   REMARK: The resulting system is tight for length *L*. In some cases, 
%   using tight system calculated for shorter *L* might work but check the
%   reconstruction error.
%
%   See also: filterbank, ufilterbank, ifilterbank

complainif_notenoughargs(nargin,2,'FILTERBANKREALTIGHT');

definput.import={'filterbankdual'};
definput.flags.outformat = {'fir','full','econ','asfreqfilter'};
definput.keyvals.efsuppthr = 10^(-5);
%definput.keyvals.bwthr = 10^(-3/10);

[flags,kv,L]=ltfatarghelper({'L'},definput,varargin,'filterbankrealtight');

[g,asan,info]=filterbankwin(g,a,L,'normal');
if isempty(L) 
    if info.isfir
        % Pick shortest possible length for FIR filterbank
        L = filterbanklength(info.longestfilter,asan);
    else
        % Just thow an error, nothing reasonable can be done without L
        error(['%s: L must be specified when not working with FIR ',...'
               'filterbanks.'], upper(mfilename));
    end
end
M=info.M;

% Force usage of the painless algorithm 
if flags.do_forcepainless
    info.ispainless = 1;
end

% Check user defined L
if L~=filterbanklength(L,a)
     error(['%s: Specified length L is incompatible with the length of ' ...
            'the time shifts.'],upper(mfilename));
end;

% Prioritize painless over uniform algorithm
if info.isuniform && info.ispainless
    info.isuniform = 0;
end

if info.isuniform
  % Uniform filterbank, use polyphase representation
  a=a(1);
  
  % Transfer functions of individual filters as cols
  G = filterbankfreqz(g,a,L);
  thisclass = class(G);
  
  N=L/a;

  gt=zeros(M,N,thisclass);
  
  for w=0:N-1
    idx_a = mod(w-(0:a-1)*N,L)+1;
    idx_b = mod((0:a-1)*N-w,L)+1;
    Ha = G(idx_a,:);
    Hb = conj(G(idx_b,:));
    
    Ha=sqrtm(Ha*Ha'+Hb*Hb')\Ha;
    
    gt(:,idx_a)=Ha.';
  end;
  % gt was created transposed because the indexing gt(:,idx_a)
  % is much faster than gt(idx_a,:)
  gt =  gt.';
  
  switch flags.outformat
      case 'fir'
          gt=ifft(gt)*sqrt(a);
          
          % Matrix cols to cell elements + cast
          gtout = cellfun(@(gtEl) cast(gtEl,thisclass), num2cell(gt,1),...
              'UniformOutput',0);
          
          % All filters in gdout will be treated as FIR of length L. Convert them
          % to a struct with .h and .offset format.
          gtout = filterbankwin(gtout,a);
      case 'full'
          
          gtout = gt*sqrt(a);
      case 'econ'
          % Shorten filters to essential support
          gt = gt*sqrt(a);
          gtout=economize_filters(gt,'efsuppthr',kv.efsuppthr);
          
      case 'asfreqfilter'
          gt = gt*sqrt(a);
          % All filters in gdout will be treated as (numeric) freqfilter format.
          % Manually convert them to a struct with .H and .foff.
          template = struct('H',[],'foff',0,'realonly',0,'delay',0,'L',L);
          gtout = cell(1,M);
          gtout(:) = {template};
          
          [H,foff,~]=economize_filters(gt,'efsuppthr',kv.efsuppthr);
          for kk = 1:M
              gtout{kk} = setfield(gtout{kk},'H',H{kk});
              gtout{kk} = setfield(gtout{kk},'foff',foff(kk));
          end
      otherwise
          error('%s: Unknown filter format.', upper(mfilename));
  end
  
elseif info.ispainless
    gtout = comp_painlessfilterbank(g,asan,L,'tight',1);

else
    error(['%s: The canonical dual frame of this system is not a ' ...
               'filterbank. You must call an iterative ' ...
               'method to perform the desired inverstion. Please see ' ...
               'FRANAITER or FRSYNITER.'],upper(mfilename));        

end;
