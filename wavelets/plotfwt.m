function [C] = plotfwt(c,info,varargin)
%PLOTFWT  Plot wavelet coefficients
%   Usage:  plotfwt(c,info,fs) 
%           plotfwt(c,info,fs,'dynrange',dynrange,...)
%
%   `plotfwt(c,g)` plots the wavelet coefficients *c* which were obtained
%   from |fwt|_ with wavelet filters definition *g*. The input cell-array
%   *c* can have two formats:
% 
%      1) Two-element cell array {c,Lc}. The first element *c* is vector of
%      coefficients in a packed format, second one *Lc* is vector of subband
%      lengths.
%    
%      2) Coefficients in cell array format *c*. Each entry of the array is
%      a single subband.
%
%   The formats are interchangable using |wavpack2cell|_ and |wavcell2pack|
%   functions.
%
%   For possible formats of the parameter *g* see |fwt|_ function.
%
%   `plotfwt(c,g,fs)` does the same plot assuming a sampling rate of *fs* Hz
%   of the original signal.
%
%   `plotfwt(c,g,fs,'dynrange',dynrange)` additionally limits the dynamic range.
%
%   `C=plotfwt(...)` returns the processed image data used in the
%   plotting. Inputting this data directly to `imagesc` or similar functions
%   will create the plot. This is usefull for custom post-processing of the
%   image data.
%
%   `plotfwt` supports optional parameters of |tfplot|_. Please see
%   the help of |tfplot|_ for an exhaustive list.
%
%   See also: fwt, tfplot

if nargin<2
  error('%s: Too few input parameters.',upper(mfilename));
end;

definput.import={'tfplot'};
definput.flags.fwtplottype = {'tfplot','stem'};
definput.keyvals.fs = [];
definput.keyvals.dynrange = [];
[flags,kv]=ltfatarghelper({'fs'},definput,varargin);

if(flags.do_stem)
   error('%s: Flag %s not supported yet.',upper(mfilename),flags.fwtplottype);
end

draw_ticks = 1;
if(strcmpi(info.fname,'fwt'))
   % Change to the cell format
   if(isnumeric(c))
       c = wavpack2cell(c,info.Lc,info.dim);
   end

   % Only one channel signals can be plotted.
   if(size(c{1},2)>1)
      error('%s: Multichannel input not supported.',upper(mfilename));
   end

   subbNo = numel(c);
   g = fwtinit(info.fwtstruct,'syn');
   aBase = g.a;
   filtNo = numel(g.filts);
   J = info.J;
   a = [aBase(1).^J, reshape(aBase(2:end)*aBase(1).^(J-1:-1:0),1,[])]';
elseif(strcmpi(info.fname,'ufwt'))
   % Only one channel signals can be plotted.
   if(ndims(c)>2)
      error('%s: Multichannel not supported.',upper(mfilename));
   end

   subbNo = size(c,2);
   a = ones(subbNo,1);

   % Determine number of levels *J* and subsampling factors *a* for
   % subbands. For correct y axis description.
   g = fwtinit(info.fwtstruct,'syn');
   filtNo = numel(g.filts);
   J = info.J; 
elseif(strcmpi(info.fname,'wfbt'))
   % Only one channel signals can be plotted.
   if(size(c{1},2)>1)
      error('%s: Multichannel input not supported.',upper(mfilename));
   end
   a = treeSub(info.wt);
   subbNo = numel(c);
   J = subbNo - 1;
   filtNo = 2;
   draw_ticks = 0;
end

% Use plotfilterbank
C=plotfilterbank(c,a,[],kv.fs,kv.dynrange,flags.plottype,...
  flags.log,flags.colorbar,flags.display,'fontsize',kv.fontsize,'clim',kv.clim);

if(draw_ticks)
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
end


