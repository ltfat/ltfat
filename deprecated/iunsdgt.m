function f=iunsdgt(c,g,a,varargin)
%IUNSDGT  Inverse uniform non-stationary discrete Gabor transform
%   Usage:  f=iunsdgt(c,g,a,Ls);
%
%   Input parameters:
%         c     : Cell array of coefficients.
%         g     : Cell array of window functions.
%         a     : Vector of time positions of windows.
%         Ls    : Length of input signal.
%   Output parameters:
%         f     : Signal.
%
%   IUNSDGT(c,g,a,Ls) computes the non-stationary Gabor expansion of the 
%   input coefficients c.
%
%   IUNSDGT is used to invert the function NSDGT. Read the help of NSDGT
%   for details of variables format and usage.
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

%   AUTHOR : Florent Jaillet
%   TESTING: TEST_NSDGT
%   REFERENCE: 
%   Last changed 2009-05

warning(['LTFAT: IUNSDGT has been deprecated, use INSDGT instead.']);  

f=insdgt(varargin{:});

  
