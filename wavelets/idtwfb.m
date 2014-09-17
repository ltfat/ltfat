function f=idtwfb(c,par,varargin)
%IDTWFB Inverse Dual-tree Filterbank
%   Usage:  f=idtwfb(c,info);
%           f=idtwfb(c,dualwt,Ls);
%
%   Input parameters:
%         c       : Input coefficients.
%         info    : Transform params. struct
%         dualwt  : Dual-tree Wavelet Filterbank definition
%         Ls      : Length of the reconstructed signal.
%
%   Output parameters:
%         f     : Reconstructed data.
%
%   `f = idtwfb(c,info)` reconstructs signal *f* from the coefficients *c*
%   using parameters from `info` struct. both returned by |dtwfb| function.
%
%   `f = idtwfb(c,dualwt,Ls)` reconstructs signal *f* from the coefficients
%   *c* using dual-tree filterbank defined by `dualwt`. Plese see |dtwfb| 
%   for supported formats. The *Ls* parameter is mandatory due to the
%   ambiguity of reconstruction lengths introduced by the subsampling
%   operation. 
%   Note that the same flag as in the |dtwfb| function have to be used, 
%   otherwise perfect reconstruction cannot be obtained. Please see help 
%   for |dtwfb| for description of the flags.
%
%   Examples:
%   ---------
%
%   A simple example showing perfect reconstruction using |idtwfb|:::
%
%      f = gspi;
%      J = 7;
%      wtdef = {'qshift3',J};
%      c = dtwfb(f,wtdef);
%      fhat = idtwfb(c,wtdef,length(f));
%      % The following should give (almost) zero
%      norm(f-fhat)
%
%   See also: dtwfb dtwfbinit
%
%


complainif_notenoughargs(nargin,2,'IDTWFB');

if ~iscell(c)
   error('%s: Unrecognized coefficient format.',upper(mfilename));
end

if(isstruct(par)&&isfield(par,'fname'))
   complainif_toomanyargs(nargin,2,'IDTWFB');
   
   if ~strcmpi(par.fname,'dtwfb')
      error(['%s: Wrong func name in info struct. ',...
             ' The info parameter was created by %s.'],...
             upper(mfilename),par.fname);
   end

   dtw = dtwfbinit({'dual',par.wt},par.fOrder);
   Ls = par.Ls;
   ext = 'per';
   L = wfbtlength(Ls,dtw,ext);
else
   complainif_notenoughargs(nargin,3,'IDTWFB');

   %% PARSE INPUT
   definput.keyvals.Ls=[];
   definput.import = {'wfbtcommon'};

   [flags,kv,Ls]=ltfatarghelper({'Ls'},definput,varargin);
   complainif_notposint(Ls,'Ls');

   ext = 'per';
   % Initialize the wavelet tree structure
   dtw = dtwfbinit(par,flags.forder);

   [Lc,L]=wfbtclength(Ls,dtw,ext);
   Lc = [Lc;Lc(end:-1:1)];

   % Do a sanity check
   if ~isequal(Lc,cellfun(@(cEl) size(cEl,1),c))
      error(['%s: The coefficient subband lengths do not comply with the'...
             ' signal length *Ls*.'],upper(mfilename));
   end
end

%% ----- step 3 : Run computation
[nodesBF, rangeLoc, rangeOut] = treeBFranges(dtw,'rev');
outLengths = nodesInLen(nodesBF,L,strcmpi(ext,'per'),dtw);
outLengths(end) = L;

f = comp_idtwfb(c,dtw.nodes(nodesBF),dtw.dualnodes(nodesBF),outLengths,...
                rangeLoc,rangeOut,ext,1);
f = postpad(f,Ls);
