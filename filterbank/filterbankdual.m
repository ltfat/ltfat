function gdout=filterbankdual(g,a,varargin)
%FILTERBANKDUAL  Dual filters
%   Usage:  gd=filterbankdual(g,a,L);
%           gd=filterbankdual(g,a);
%           
%
%   `filterbankdual(g,a,L)` computes the canonical dual filters of *g* for a
%   channel subsampling rate of *a* (hop-size) and system length *L*.
%   *L* must be compatible with subsampling rate *a* as 
%   `L==filterbanklength(L,a)`.
%
%   `filterabankrealdual(g,a)` does the same, but the filters must be FIR
%   filters, as the transform length is unspecified. *L* will be set to 
%   next suitable length equal or bigger than the longest impulse response.
%
%   The input and output format of the filters *g* are described in the
%   help of |filterbank|.
%
%   In addition, the funtion recognizes a 'forcepainless' flag which
%   forces treating the filterbank *g* and *a* as a painless case
%   filterbank.  
%
%   To actually invert the output of a filterbank, use the dual filters
%   together with the |ifilterbank| function.
%
%   REMARK: Perfect reconstruction can be obtained for signals of length
%   *L*. In some cases, using dual system calculated for shorter *L* might
%   work but check the reconstruction error.
%
%   See also: filterbank, ufilterbank, ifilterbank

complainif_notenoughargs(nargin,2,'FILTERBANKDUAL');

definput.import={'filterbankdual'};
[flags,~,L]=ltfatarghelper({'L'},definput,varargin);

if isempty(L)
    if ~all(cellfun(@(gEl) isfield(gEl,'H'),g))
        % All filters are FIR, therefore filterbankwin can be called without L
        [~,info]=filterbankwin(g,a);
        if ~info.isfir
            % Just a sanity check
            error('%s: Internal error. Filterbank should be FIR. ',...
                  upper(mfilename));
        end
        % Use next suitable length
        L = filterbanklength(info.longestfilter,a);
    else
        error(['%s: L must be specified when working with filters defined ',...
           ' in frequency.'], upper(mfilename));
   end
end

[g,info]=filterbankwin(g,a,L,'normal');
M=info.M;

% Force usage of the painless algorithm 
if flags.do_forcepainless
    info.ispainless = 1;
end

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
  
  % G1 is done this way just so that we can determine the data type.
  G1=comp_transferfunction(g{1},L);
  thisclass=assert_classname(G1);
  G=zeros(L,M,thisclass);
  G(:,1)=G1;
  for ii=2:M
    G(:,ii)=comp_transferfunction(g{ii},L);
  end;
  
  N=L/a;
  
  gd=zeros(M,N,thisclass);
  
  for w=0:N-1
    idx = mod(w-(0:a-1)*N,L)+1;
    H = G(idx,:);
    
    H=pinv(H)';
    
    gd(:,idx)=H.';
  end;
  % gd was created transposed because the indexing gd(:,idx_a)
  % is much faster than gd(idx_a,:)
  gd =  gd.';
  
  gd=ifft(gd)*a;
  
  % Matrix cols to cell elements + cast
  gdout = cellfun(@(gdEl) cast(gdEl,thisclass), num2cell(gd,1),...
                  'UniformOutput',0);
  % All filters in gdout will be treated as FIR of length L. Convert them
  % to a struct with .h and .offset format.
  gdout = filterbankwin(gdout,a); 
  
elseif info.ispainless
   gdout = comp_painlessfilterbank(g,info.a,L,'dual',0);
else
        error(['%s: The canonical dual frame of this system is not a ' ...
               'filterbank. You must either call an iterative ' ...
               'method to perform the desired inverstion or transform ',...
               'or transform the filterbank to uniform one. Please see ' ...
               'FRANAITER or FRSYNITER for the former and ',...
               'NONU2UFILTERBANK for the latter case.'],upper(mfilename));        

    
end;
