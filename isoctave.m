function t=isoctave()
%ISOCTAVE  True if the operating environment is octave.
%   Usage: t=isoctave();
%
%   ISOCTAVE returns 1 if the operating environment is Octave, otherwise
%   0 (Matlab)

%   AUTHOR : Peter Soendergaard.  
%   TESTING: NA
%   REFERENCE: NA
persistent inout;

if isempty(inout),
  inout = exist('OCTAVE_VERSION','builtin') ~= 0;
end;
t = inout;

%OLDFORMAT
