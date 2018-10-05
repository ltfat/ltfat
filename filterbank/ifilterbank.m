function [f,Ls]=ifilterbank(c,g,a,varargin)
%IFILTERBANK  Filter bank inversion
%   Usage:  f=ifilterbank(c,g,a);
%
%   `ifilterbank(c,g,a)` synthesizes a signal *f* from the coefficients *c*
%   using the filters stored in *g* for a channel subsampling rate of *a* (the
%   hop-size). The coefficients has to be in the format returned by
%   either |filterbank| or |ufilterbank|.
%
%   The filter format for *g* is the same as for |filterbank|.
%
%   If perfect reconstruction is desired, the filters must be the duals
%   of the filters used to generate the coefficients. See the help on
%   |filterbankdual|.
%
%   Additional parameters
%   ---------------------
%
%   'complex' (default), 'real'
%       The 'real' flag indicates that the filters *g* cover only the positive
%       frequencies and does `2*real(f)` to effectivelly mirror the filters to
%       cover also the negative frequencies.
%       
%   See also: filterbank, ufilterbank, filterbankdual
%
%   References: bohlfe02

if nargin<3
  error('%s: Too few input parameters.',upper(mfilename));
end;

definput.import={'pfilt'};
definput.keyvals.Ls=[];
definput.flags.complex={'complex','real'};
[flags,kv,Ls]=ltfatarghelper({'Ls'},definput,varargin);

L=filterbanklengthcoef(c,a);


if iscell(c)
  M=numel(c);
else
  M=size(c,2);
end;

[g,asan]=filterbankwin(g,a,L,'normal');

if numel(g)~=M
  error(['%s: Number of filters must be equal to the number of channels ' ...
            'of coefficients.'],upper(mfilename));
end

 if size(a,1)>1 
   if  size(a,1)~=M
     error(['%s: The number of entries in "a" must match the number of ' ...
            'filters.'],upper(mfilename));
   end;
 end;

g = comp_filterbank_pre(g,asan,L,kv.crossover);

% Handle ufilterbank output format here
if isnumeric(c)
   ctmp = c;
   c = cell(M,1);
   for m=1:M
      c{m}=squeeze(ctmp(:,m,:));
   end;
end


f = comp_ifilterbank(c,g,asan,L);

if flags.do_real
    f = 2*real(f);
end

% Cut or extend f to the correct length, if desired.
if ~isempty(Ls)
  f=postpad(f,Ls);
else
  Ls=L;
end;
