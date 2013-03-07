function [C] = plotufwt(c,varargin)
%PLOTFWT  Plot wavelet coefficients
%   Usage:  plotufwt(c,g,fs) 
%
%   `plotufwt(c,g)` plots the wavelet coefficients *c* which were obtained
%   from |ufwt|_ with wavelet filters definition *g*. 
%
%   For possible formats of the parameter *g* see |fwt|_ function.
%
%   `plotufwt(c,g,fs)` does the same plot assuming a sampling rate of *fs* Hz
%   of the original signal.
%
%   `plotufwt(c,g,fs,'dynrange',dynrange)` additionally limits the dynamic range.
%
%   `C=plotufwt(...)` returns the processed image data used in the
%   plotting. Inputting this data directly to `imagesc` or similar functions
%   will create the plot. This is usefull for custom post-processing of the
%   image data.
%
%   `plotfwt` supports optional parameters of |tfplot|_. Please see
%   the help of |tfplot|_ for an exhaustive list.
%
%   See also: fwt, tfplot

if nargin<1
  error('%s: Too few input parameters.',upper(mfilename));
end;

definput.import={'tfplot'};
definput.flags.fwtplottype = {'tfplot','stem'};
definput.keyvals.g = [];
definput.keyvals.fs = [];
definput.keyvals.dynrange = [];
[flags,kv]=ltfatarghelper({'g','fs'},definput,varargin);

% Change to the cell format
if(iscell(c) && size(c,2)==2 && size(c,1)==1)
    c = wavpack2cell(c{1},c{2});
end

% Only one channel signals can be plotted.
if(ndims(c)>2)
   error('%s: Multichannel input.',upper(mfilename));
end

subbNo = size(c,2);
a = ones(subbNo,1);
if(~isempty(kv.g))
   % Determine number of levels *J* and subsampling factors *a* for
   % subbands. For correct y axis description.
   g = fwtinit(kv.g,'syn');
   filtNo = numel(g.filts);
   J = (subbNo-1)/(filtNo-1);
else
   % general values
   filtNo = 2;
   J = subbNo-1;
   a = ones(subbNo,1);
end


% Use plotfilterbank
C=plotfilterbank(c,a,[],kv.fs,kv.dynrange,flags.plottype,...
  flags.log,flags.colorbar,flags.display,'fontsize',kv.fontsize,'clim',kv.clim);

% Redo the yticks and ylabel
yTickLabels = cell(1,subbNo);
yTickLabels{1} = sprintf('a%d',J);
Jtmp = ones(filtNo-1,1)*(J:-1:1);
for ii=1:subbNo-1
   yTickLabels{ii+1} = sprintf('d%d',Jtmp(ii));
end

ylabel('Subbands','fontsize',kv.fontsize);
set(gca,'ytick',1:subbNo);
set(gca,'ytickLabel',yTickLabels,'fontsize',kv.fontsize);


