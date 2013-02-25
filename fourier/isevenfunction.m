function t=isevenfunction(f,varargin);
%ISEVENFUNCTION  True if function is even
%   Usage:  t=isevenfunction(f);
%           t=isevenfunction(f,tol);
%
%   `isevenfunction(f)` returns 1 if f is whole point even. Otherwise it
%   returns 0.
%
%   `isevenfunction(f,tol)` does the same, using the tolerance *tol* to measure
%   how large the error between the two parts of the vector can be. Default
%   is 1e-10.
%
%   Adding the flag `'hp'` as the last argument does the same for half point
%   even functions.
%  
%   See also: middlepad, peven
  
%   AUTHOR : Peter L. Søndergaard
%   TESTING: OK
%   REFERENCE: OK

if nargin<1
  error('Too few input parameters.');
end;

if size(f,2)>1
  if size(f,1)>1
    error('f must be a vector');
  else
    % f was a row vector.
    f=f(:);
  end;
end;

% Define initial values for flags
definput.flags.centering = {'wp','hp'};
definput.keyvals.tol     = 1e-10; 

[flags,keyvals,tol]=ltfatarghelper({'tol'},definput,varargin);

L=size(f,1);

if flags.do_wp
  % Determine middle point of sequence.
  if rem(L,2)==0
    middle=L/2;
  else
    middle=(L+1)/2;
  end;
  
  % Relative norm of difference between the parts of the signal.
  d=norm(f(2:middle)-conj(flipud(f(L-middle+2:L))))/norm(f);
else
  
  middle=floor(L/2);
  
  d=norm(f(1:middle)-conj(flipud(f(L-middle+1:L))))/norm(f);

end;

% Return true if d less than tolerance.
t=d<=tol;
