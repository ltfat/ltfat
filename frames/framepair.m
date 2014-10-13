function [F1,F2]=framepair(ftype,g1,g2,varargin)
%FRAMEPAIR  Construct a new frame
%   Usage: [F1,F2]=framepair(ftype,g1,g2,...);
%
%   `[F1,F2]=framepair(ftype,g1,g2,...)` constructs two new frame objects 
%   *F1* and *F2* of the same type *ftype* using the windows *g1* and *g2*.
%   The windows are specific to choosen frame type. See the help on *frame*
%   for the windows and arguments. 
%
%   This function makes it easy to create a pair of canonical dual frames:
%   simply specify `'dual'` as window if one frame should be the dual of the
%   other.
%
%   This is most easily explained through some examples. The following
%   example creates a Gabor frame for real-valued signals with a Gaussian
%   analysis window and its canonical dual frame as the synthesis frame:::
%
%      f=greasy;
%      [Fa,Fs]=framepair('dgtreal','gauss','dual',20,294);
%      c=frana(Fa,f);
%      r=frsyn(Fs,c);
%      norm(f-r)
%
%   The following example creates a Wilson basis with a Gaussian
%   synthesis window, and its canonical dual frame as the analysis
%   frame::
% 
%     [Fa,Fs]=framepair('dwilt','dual','gauss',20);
%
%   See also: frame, framedual, frametight

  
complainif_notenoughargs(nargin,3,'FRAMEPAIR');

ftype=lower(ftype);

if ~strcmp(g1,'dual')
    F1=frame(ftype,g1,varargin{:});
end;

if ~strcmp(g2,'dual')
    F2=frame(ftype,g2,varargin{:});
end;

if strcmp(g1,'dual')
    F1=framedual(F2);
end;

if strcmp(g2,'dual')
    F2=framedual(F1);
end;
