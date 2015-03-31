function f=iunsdgtreal(c,g,a,M,Ls)
%IUNSDGTREAL  Inverse uniform non-stationary discrete Gabor transform
%   Usage:  f=iunsdgtreal(c,g,a,M,Ls);
%
%   Input parameters:
%         c     : Cell array of coefficients.
%         g     : Cell array of window functions.
%         a     : Vector of time positions of windows.
%         M     : Numbers of frequency channels.
%         Ls    : Length of input signal.
%   Output parameters:
%         f     : Signal.
%
%   IUNSDGTREAL(c,g,a,M,Ls) computes the inverse uniform non-stationary Gabor
%   expansion of the input coefficients c.
%
%   IUNSDGTREAL is used to invert the function UNSDGTREAL. Read the help of
%   UNSDGTREAL for details of variables format and usage.
%
%   For perfect reconstruction, the windows used must be dual windows of 
%   the ones used to generate the coefficients. The windows can be
%   generated unsing NSGABDUAL.
%
%   See also:  unsdgt, nsgabdual, nsgabtight
%
%   Demos:  demo_nsdgt
%
%   References: ltfatnote018

%   AUTHOR : Florent Jaillet and Peter L. Søndergaard
%   TESTING: TEST_NSDGT
%   REFERENCE: OK


warning(['LTFAT: IUNSDGTREAL has been deprecated, use INSDGTREAL instead.']);
  
f=insdgtreal(varargin{:});
