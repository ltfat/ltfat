function o = frgramian(c, Fa, Fs)
%FRGRAMIAN Frame Gramian operator
%   Usage:  o=frgramian(c, F);
%           o=frgramian(c, Fa, Fs);
%
%   Input parameters:
%          c    : Input coefficients
%          Fa   : Analysis frame
%          Fs   : Synthesis frame
%
%   Output parameters: 
%          o    : Output coefficients
%     
%   `o=frgramian(c,F)` applies the Gramian operator or Gram matrix of the 
%   frame *F*. The entries of the Gram matrix are the inner products of the 
%   frame elements of *F*. The frame must have been created using |frame|.
%   If the frame *F* is a Parseval frame, the Gramian operator is a projection 
%   onto the range of the frame analysis operator.
%
%   `o=frgramian(c, Fa, Fs)` applies the (cross) Gramian operator with the 
%   frames *Fa* and *Fs*. Here *Fs* is the frame associated with the frame
%   synthesis operator and *Fa* the frame that is associated with the 
%   frame analysis operator. The entries of the matrix that is constructed
%   through the Gramian operator are the inner products of the frame 
%   elements of *Fa* and *Fs*.
%   If *Fa* and *Fs* are canonical dual frames, the Gramian operator is a 
%   projection onto the range of the frame analysis operator.
%

% AUTHOR: Jordy van Velthoven

complainif_notenoughargs(nargin, 2, 'FRGRAMIAN');
complainif_notvalidframeobj(Fa,'FRGRAMIAN');

if (nargin == 2)
   Fs = Fa;
else
   complainif_notvalidframeobj(Fs,'FRGRAMIAN'); 
end;

o = frana(Fa, frsyn(Fs, c));
