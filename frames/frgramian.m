function o = frgramian(c, Fa, Fs);
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
%   `o=frgramian(c,F)` computes the Gramian operator or Gram matrix of the 
%   frame *F*. The entries of the Gram matrix are the inner products of the 
%   frame elements of *F*. The frame must have been created using |frame|.
%   If the frame *F* is a Parseval frame, the Gramian operator is a projection 
%   onto the range of the frame analysis operator.
%
%   `o=frgramian(o, Fa, Fs)` computes the Gramian operator with the frames *Fa* 
%   and *Fs*. Here is *Fs* the frame associated with the frame synthesis operator 
%   and *Fa* the frame that is associated with the frame analysis operator. The 
%   entries of the matrix that is constructed through the Gramian operator are the 
%   inner products of the frame elements of *Fa* and *Fs*. The frames *Fa* and 
%   *Fs* must have been created using |framepair|.
%   If *Fa* and *Fs* are canonical dual frames, the Gramian operator is a projection
%   onto the range of the frame analysis operator.
%
% AUTHOR: Jordy van Velthoven

complainif_notenoughargs(nargin, 2, 'FRGRAMIAN');

if (nargin == 2)
  Fs = Fa;
end;

o = frana(Fa, frsyn(Fs, c));
