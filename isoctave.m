function t=isoctave()
%ISOCTAVE  True if the operating environment is octave
%   Usage: t=isoctave();
%
%   `isoctave` returns 1 if the operating environment is Octave, otherwise
%   it returns 0 (Matlab)

%   AUTHOR : Peter L. Søndergaard.  
%   TESTING: NA
%   REFERENCE: NA
persistent inout;

if isempty(inout),
  inout = exist('OCTAVE_VERSION','builtin') ~= 0;
end;
t = inout;
