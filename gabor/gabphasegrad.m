function [tgrad,fgrad,c]=gabphasegrad(method,varargin)
%GABPHASEGRAD   Phase gradient of the DGT
%   Usage:  [tgrad,fgrad,c] = gabphasegrad('dgt',f,g,a,M);
%           [tgrad,fgrad]   = gabphasegrad('phase',cphase,a);
%           [tgrad,fgrad]   = gabphasegrad('abs',s,g,a);
%
%   `[tgrad,fgrad]=gabphasegrad(method,...)` computes the time-frequency
%   gradient of the phase of the |dgt| of a signal. The derivative in time
%   *tgrad* is the instantaneous frequency while the frequency derivative
%   *fgrad* is the local group delay.
%
%   *tgrad* and *fgrad* measure the deviation from the current time and
%   frequency, so a value of zero means that the instantaneous frequency is
%   equal to the center frequency of the considered channel.
%
%   *tgrad* is scaled such that distances are measured in samples. Similarly,
%   *fgrad* is scaled such that the Nyquest frequency (this highest possible
%   frequency) corresponds to a value of L/2.
%
%   The computation of *tgrad* and *fgrad* is inaccurate when the absolute
%   value of the Gabor coefficients is low. This is due to the fact the the
%   phase of complex numbers close to the machine precision is almost
%   random. Therefore, *tgrad* and *fgrad* may attain very large random values
%   when `abs(c)` is close to zero.
%
%   The computation can be done using three different methods.
%
%     'dgt'    Directly from the signal. This is the default method.
% 
%     'phase'  From the phase of a DGT of the signal. This is the
%              classic method used in the phase vocoder.
%
%     'abs'    From the absolute value of the DGT. Currently this
%              method works only for Gaussian windows.
%
%   `[tgrad,fgrad]=gabphasegrad('dgt',f,g,a,M)` computes the time-frequency
%   gradient using a DGT of the signal *f*. The DGT is computed using the
%   window *g* on the lattice specified by the time shift *a* and the number
%   of channels *M*. The algorithm used to perform this calculation computes
%   several DGTs, and therefore this routine takes the exact same input
%   parameters as |dgt|.
%
%   The window *g* may be specified as in |dgt|. If the window used is
%   'gauss', the computation will be done by a faster algorithm.
%
%   `[tgrad,fgrad,c]=gabphasegrad('dgt',f,g,a,M)` additionally returns the
%   Gabor coefficients *c*, as they are always computed as a byproduct of the
%   algorithm.
%
%   `[tgrad,fgrad]=gabphasegrad('phase',cphase,a)` computes the phase
%   gradient from the phase *cphase* of a DGT of the signal. The original DGT
%   from which the phase is obtained must have been computed using a
%   time-shift of *a*.
%
%   `[tgrad,fgrad]=gabphasegrad('abs',s,g,a)` computes the phase gradient
%   from the spectrogram *s*. The spectrogram must have been computed using
%   the window *g* and time-shift *a*.
%
%   `[tgrad,fgrad]=gabphasegrad('abs',s,g,a,difforder)` uses a centered finite
%   diffence scheme of order *difforder* to perform the needed numerical
%   differentiation. Default is to use a 4'th order scheme.
%
%   Currently the `'abs'` method only works if the window *g* is a Gaussian
%   window specified as a string or cell array.
%
%   See also: resgram, gabreassign, dgt
%
%   References: aufl95 cmdaaufl97 fl65

  
% AUTHOR: Peter L. Søndergaard, 2008.
  
%error(nargchk(4,6,nargin));

if ~ischar(method)
    error(['First argument must be the method name, "dgt", "phase" or ' ...
           '"abs".']);
end;
  
switch(lower(method))
 case 'dgt'
  % ---------------------------  DGT method ------------------------

  [f,g,a,M]=deal(varargin{1:4});
  
  definput.keyvals.L=[];
  definput.keyvals.minlvl=eps;
  definput.keyvals.lt=[0 1];
  [flags,kv,L,minlvl]=ltfatarghelper({'L','minlvl'},definput,varargin(5:end));
  
  
  %% ----- step 1 : Verify f and determine its length -------
  % Change f to correct shape.
  [f,Ls,W,wasrow,remembershape]=comp_sigreshape_pre(f,upper(mfilename),0);
  
  %% ------ step 2: Verify a, M and L
  if isempty(L)
      
      % ----- step 2b : Verify a, M and get L from the signal length f----------
      L=dgtlength(Ls,a,M,kv.lt);
      
  else
      
      % ----- step 2a : Verify a, M and get L
      Luser=dgtlength(L,a,M,kv.lt);
      if Luser~=L
          error(['%s: Incorrect transform length L=%i specified. Next valid length ' ...
                 'is L=%i. See the help of DGTLENGTH for the requirements.'],...
                upper(mfilename),L,Luser);
      end;
      
  end;
  
  %% ----- step 3 : Determine the window 
  
  [g,info]=gabwin(g,a,M,L,kv.lt,'callfun',upper(mfilename));
  
  if L<info.gl
      error('%s: Window is too long.',upper(mfilename));
  end;
  
  %% ----- step 4: final cleanup ---------------
  
  f=postpad(f,L);
  
  %% ------ algorithm starts --------------------
  
  % Compute the time weighted version of the window.
  hg=fftindex(L).*g;
  
  % The computation done this way is insensitive to whether the dgt is
  % phaselocked or not.
  c   = comp_dgt(f,g,a,M,kv.lt,0,0,0);
  c_h = comp_dgt(f,hg,a,M,kv.lt,0,0,0);
  
  c_s = abs(c).^2;
  
  % Remove small values because we need to divide by c_s
  c_s = max(c_s,minlvl*max(c_s(:)));
  
  % Compute the group delay
  fgrad=real(c_h.*conj(c)./c_s);
  
  if info.gauss
    % The method used below only works for the Gaussian window, because the
    % time derivative and the time multiplicative of the Gaussian are identical.
    tgrad=imag(c_h.*conj(c)./c_s)/info.tfr;
  else
    
    % The code below works for any window, and not just the Gaussian
    
    dg  = pderiv(g,[],Inf)/(2*pi);
    c_d = comp_dgt(f,dg,a,M,kv.lt,0,0,0);
    c_d = reshape(c_d,M,N,W);
    
    % Compute the instantaneous frequency
    tgrad=-imag(c_d.*conj(c)./c_s);
  end;
  
  
 case 'phase'
  % ---------------------------  phase method ------------------------

  [cphase,a]=deal(varargin{1:2});
  
  if ~isreal(cphase)
    error(['Input phase must be real valued. Use the "angle" function to ' ...
           'compute the argument of complex numbers.']);
  end;
  
  % --- linear method ---
  [M,N,W]=size(cphase);
  L=N*a;
  b=L/M;
  
  if 0
    
    % This is the classic phase vocoder algorithm by Flanagan.
    
    tgrad = cphase-circshift(cphase,[0,-1]);
    tgrad = tgrad- 2*pi*round(tgrad/(2*pi));
    tgrad = -tgrad/(2*pi)*L;
    
    % Phase-lock the angles.
    TimeInd = (0:(N-1))*a;
    FreqInd = (0:(M-1))/M;
    
    phl = FreqInd'*TimeInd;
    cphase = cphase+2*pi.*phl;
    
    fgrad = cphase-circshift(cphase,[1,0]);
    fgrad = fgrad- 2*pi*round(fgrad/(2*pi));
    fgrad = -fgrad/(2*pi)*L;
    
  end;
  
  
  if 1
    % This is the classic phase vocoder algorithm by Flanagan modified to
    % yield a second order centered difference approximation.
    
    % Forward approximation
    tgrad_1 = cphase-circshift(cphase,[0,-1]);
    tgrad_1 = tgrad_1 - 2*pi*round(tgrad_1/(2*pi));
    % Backward approximation
    tgrad_2 = circshift(cphase,[0,1])-cphase;
    tgrad_2 = tgrad_2 - 2*pi*round(tgrad_2/(2*pi));
    % Average
    tgrad = (tgrad_1+tgrad_2)/2;
    
    tgrad = -tgrad/(2*pi*a)*L;
    
    % Phase-lock the angles.
    TimeInd = (0:(N-1))*a;
    FreqInd = (0:(M-1))/M;
    
    phl = FreqInd'*TimeInd;
    cphase = cphase+2*pi.*phl;
    
    % Forward approximation
    fgrad_1 = cphase-circshift(cphase,[-1,0]);
    fgrad_1 = fgrad_1 - 2*pi*round(fgrad_1/(2*pi));
    % Backward approximation
    fgrad_2 = circshift(cphase,[1,0])-cphase;
    fgrad_2 = fgrad_2 - 2*pi*round(fgrad_2/(2*pi));
    % Average
    fgrad = (fgrad_1+fgrad_2)/2;
    
    fgrad = fgrad/(2*pi*b)*L;
    
  end;
  
  
 case 'abs'
  % ---------------------------  abs method ------------------------

  [s,g,a]=deal(varargin{1:3});
  if numel(varargin)>3
    difforder=varargin{4};
  else
    difforder=4;
  end;

  if ~(all(s(:)>=0))
    error('First input argument must be positive or zero.');
  end;
    
  [M,N,W]=size(s);
  
  L=N*a;
  tfr=1;

  g
  
  [g,info]=gabwin(g,a,M,L,'callfun','GABPHASEGRAD');

  info
  
  if ~info.gauss
    error(['The window must be a Gaussian window (specified as a string or ' ...
           'as a cell arrray).']);
  end;
  
  L=N*a;
  b=L/M;
  
  % We must avoid taking the log of zero.
  % Therefore we add the smallest possible
  % number
  logs=log(s+realmin);
  
  % XXX REMOVE Add a small constant to limit the dynamic range. This should
  % lessen the problem of errors in the differentation for points close to
  % (but not exactly) zeros points.
  maxmax=max(logs(:));
  tt=-11;
  logs(logs<maxmax+tt)=tt;
  
  fgrad=pderiv(logs,2,difforder)/(2*pi)*info.tfr;
  tgrad=pderiv(logs,1,difforder)/(2*pi*info.tfr);
  
  
 otherwise
  error(['First argument must be the method name, "dgt", "phase" or ' ...
         '"abs".']);
end;  
