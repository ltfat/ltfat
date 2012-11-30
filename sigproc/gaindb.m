function inoutsig = gaindb(inoutsig,gn,varargin)
%GAINDB  Increase/decrease level of signal
%   Usage:  outsig = gaindb(insig,gn);
%
%   `gaindb(insig,gn)` increases the energy level of the signal by *gn*
%   dB.
%
%   If *gn* is a scalar, the whole input signal is scaled.
%
%   If *gn* is a vector, each column is scaled by the entries in
%   *gn*. The length of *gn* must match the number of columns.
%
%   `gaindb(insig,gn,dim)` scales the signal along dimension *dim*.
%
%   See also: rms

%   AUTHOR: Peter L. Søndergaard, 2009

% ------ Checking of input parameters ---------
  
if nargin<2
  error('%s: Too few input parameters.',upper(mfilename));
end;

if ~isnumeric(inoutsig)
  error('%s: insig must be numeric.',upper(mfilename));
end;

if ~isnumeric(gn) 
  error('%s: gn must be numeric.',upper(mfilename));
end;

definput.keyvals.dim=[];
[flags,kv]=ltfatarghelper({'dim'},definput,varargin);


% ------ Computation --------------------------

if isscalar(gn)
  inoutsig = inoutsig*10^(gn/20);
else
  if isvector(gn)
    M=length(gn);
        
    [inoutsig,L,Ls,W,dim,permutedsize,order]=...
        assert_sigreshape_pre(inoutsig,[],kv.dim,upper(mfilename));
      
    if M~=W
      error('%s: Length of gn and signal size must match.',upper(mfilename));
    end;

    for ii=1:W
      inoutsig(:,ii)=inoutsig(:,ii)*10^(gn(ii)/20);
    end;
    
    inoutsig=assert_sigreshape_post(inoutsig,kv.dim,permutedsize,order);     
    
  else
    if ~isnumeric(gn) 
      error('%s: gn must be a scalar or vector.',upper(mfilename));
    end;
  end;
end;
