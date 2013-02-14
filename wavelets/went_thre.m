function E = went_thre(x,th)
%WENT_THRE Threshold entropy
%   Usage:  E = went_compn(x,th)
%
%   `E = went_thre(x,th)` calculates entropy of the input *x* using formula
%   XXX
%
E = sum(x(abs(x)>th));