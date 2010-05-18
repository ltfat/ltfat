function [xo,N]=thresh(ttype,xi,lambda,mtype);
%THRESH   Threshold (hard/soft)
%   Usage:  x=thresh(ttype,x,lambda);
%           [x,N]=thresh(ttype,x,lambda);
%           x=thresh(ttype,x,lambda,mtype);
%           [x,N]=thresh(ttype,x,lambda,mtype);
%
%   THRESH('hard',x,lambda) will perform hard thresholding on x, i.e. all
%   element with absolute value less than lambda will be set to zero.
%
%   THRESH('soft',x,lambda) will perform soft thresholding on x, i.e. lambda
%   will be subtracted from the absolute value of every element of x.
%
%   [x,N]=THRESH(ttype,x,lambda) additionally returns a number N specifying
%   how many numbers where kept.
%
%   THRESH(ttype,x,lambda,'full') returns the output as a full matrix. This
%   is the default.
%
%   THRESH(ttype,x,lambda,'sparse') returns the output as a sparse matrix.
%
%   The function WTHRESH in the Matlab Wavelet toolbox implements the same
%   functionality.
%
%   See also: largestr, largestn

%   AUTHOR : Peter Soendergaard and Bruno Torresani.  
%   TESTING: OK
%   REFERENCE: OK

error(nargchk(3,4,nargin));

if (prod(size(lambda))~=1 || ~isnumeric(lambda))
  error('lambda must be a scalar.');
end;

if ischar(ttype)
  ttype=lower(ttype);
end;

dosparse=0;
if nargin==4
  switch(lower(mtype))
    case {'full'}
    case {'sparse'}
      if ndims(xi)>2
	error('Sparse output is only supported for 1D/2D input. This is a limitation of Matlab/Octave.');
      end;
      dosparse=1;      
    otherwise
      error('The output type (last argument) must be either "full" or "sparse".');
  end;
end;

dohard=1;
switch(lower(ttype))
 case {'hard'}
 case {'soft'}
  dohard=0;
 otherwise
  error('Unknown thresholding type.');
end;


if dosparse
  xo=sparse(size(xi,1),size(xi,2));
else
  xo=zeros(size(xi));
end;

signifmap=find(abs(xi)>=lambda);

if dohard
  xo(signifmap)=xi(signifmap);
else
  %    xo(signifmap)=xi(signifmap) - sign(xi(signifmap))*lambda;
  xo(signifmap)=(abs(xi(signifmap)) - lambda) .* ...
      exp(i*angle(xi(signifmap)));
  % The line above produces very small imaginary values when the input
  % is real-valued. The next line fixes this
  if isreal(xi)
    xo=real(xo);
  end;
  
end;

if nargout==2
  N=numel(signifmap);
end;

