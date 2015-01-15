function h = frameoperator(F, f);
%FRAMEOPERATOR Frame Operator
%   Usage:  o=frameoperator(F, f);
%
%   Input parameters:
%          F    : frame
%          f    : input vector
%
%   Output parameter: 
%          h    : output vector
%     
%   `h=frameoperator(F,f)` applies the frame operator associated with the frame 
%   *F* to the input *f*.
%
%   If the frame *F* is a tight frame, then *h* equals *f* up to the constant 
%   $\frac{1}{A}$ where $A$ is the lower frame bound of *F*. If the frame *F*
%   is an orthonormal basis, or more general a Parseval frame, then *h* equals 
%   *f*. 
%


% AUTHOR: Jordy van Velthoven

complainif_notenoughargs(nargin, 2, 'FRAMEOPERATOR');


h = frsyn(F, (frana(F,f)));
