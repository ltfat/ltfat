function f=comp_sigreshape_post(f,fl,wasrow,remembershape)
%COMP_SIGRESHAPE_POST
%

%   AUTHOR : Peter L. Søndergaard.
%   TESTING: OK
%   REFERENCE: OK

% Get original dimensionality
fd=length(remembershape);

if fd>2
  f=reshape(f,[fl,remembershape(2:fd)]);
else
  if wasrow
    f=f.';
  end;
end;


