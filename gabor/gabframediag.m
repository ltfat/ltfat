function d=gabframediag(g,a,M,L,varargin)
%GABFRAMEDIAG  Diagonal of Gabor frame operator
%   Usage:  d=gabframediag(g,a,M,L);
%           d=gabframediag(g,a,M,L,'lt',lt);
%
%   Input parameters:
%         g     : Window function.
%         a     : Length of time shift.
%         M     : Number of channels.
%         L     : Length of transform to do.
%         lt    : Lattice type (for non-separable lattices).
%   Output parameters:
%         d     : Diagonal stored as a column vector
%
%   `gabframediag(g,a,M,L)` computes the diagonal of the Gabor frame operator
%   with respect to the window *g* and parameters *a* and *M*. The
%   diagonal is stored a as column vector of length *L*.
%
%   The diagonal of the frame operator can for instance be used as a
%   preconditioner.
%
%   See also: dgt

if nargin<4
  error('%s: Too few input parameters.',upper(mfilename));
end;

definput.keyvals.lt=[0 1];
definput.flags.phase={'freqinv','timeinv'}; 
[flags,kv]=ltfatarghelper({},definput,varargin);

% ----- step 2a : Verify a, M and get L

Luser=dgtlength(L,a,M,kv.lt);
if Luser~=L
    error(['%s: Incorrect transform length L=%i specified. Next valid length ' ...
           'is L=%i. See the help of DGTLENGTH for the requirements.'],...
          upper(mfilename),L,Luser);
end;


%% ----- step 3 : Determine the window 

[g,info]=gabwin(g,a,M,L,kv.lt,'callfun',upper(mfilename));

if L<info.gl
  error('%s: Window is too long.',upper(mfilename));
end;

%%  compute the diagonal 

glong2=abs(fir2long(g,L)).^2;
N=L/a;

d=zeros(L,1,assert_classname(glong2));

% The diagonal is a-periodic, so compute a single period by summing up
% glong2 in slices. 
d=repmat(sum(reshape(glong2,a,N),2),N,1)*M;


