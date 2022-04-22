function outsig=frana(F,insig);
%FRANA  Frame analysis operator
%   Usage: c=frana(F,f);
%
%   `c=frana(F,f)` computes the frame coefficients *c* of the input
%   signal *f* using the frame *F*. The frame object *F* must have been
%   created using |frame| or |framepair|.
%
%   If *f* is a matrix, the transform will be applied along the columns
%   of *f*. If *f* is an N-D array, the transform will be applied along
%   the first non-singleton dimension.
%
%   The output coefficients are stored as columns. This is usually
%   **not** the same format as the 'native' format of the frame. As an
%   examples, the output from |frana| for a gabor frame cannot be
%   passed to |idgt| without a reshape.
%
%   Examples:
%   ---------
%
%   In the following example the signal *bat* is analyzed through a wavelet 
%   frame. The result are the frame coefficients associated with the input  
%   signal *bat* and the analysis frame `'fwt'`:::
%
%      f = bat;
%      w = 'sym8';
%      J = 7;
%      F = frame('fwt', w, J); 
%      c = frana(F, f);
%      % A plot of the frame coefficients
%      plotframe(F, c, 'dynrange', 100);
%
%   See also: frame, framepair, frsyn, plotframe

complainif_notenoughargs(nargin,2,'FRANA');
complainif_notvalidframeobj(F,'FRANA');

if size(insig,1) == 1
%    error('%s: Currently only column vectors are supported. See bug #73.',...
%          upper(mfilename)); 
  insig = insig(:);
end


%% ----- step 1 : Verify f and determine its length -------
% Change f to correct shape.
[insig,~,Ls,W,dim,permutedsize,order]=assert_sigreshape_pre(insig,[],[],upper(mfilename));
 
F=frameaccel(F,Ls);

insig=postpad(insig,F.L);

%% ----- do the computation ----

outsig=F.frana(insig);

%% --- cleanup -----

permutedsize=[size(outsig,1),permutedsize(2:end)];

outsig=assert_sigreshape_post(outsig,dim,permutedsize,order);

  
