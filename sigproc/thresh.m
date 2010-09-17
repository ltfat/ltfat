function [xo,N]=thresh(xi,lambda,varargin);
%THRESH   Threshold (hard/soft)
%   Usage:  x=thresh(x,lambda,...);
%           [x,N]=thresh(x,lambda,...);
%
%   THRESH(x,lambda) will perform hard thresholding on x, i.e. all
%   element with absolute value less than lambda will be set to zero.
%
%   THRESH(lambda,'soft') will perform soft thresholding on x, i.e. lambda
%   will be subtracted from the absolute value of every element of x.
%
%   [x,N]=THRESH(x,lambda) additionally returns a number N specifying
%   how many numbers where kept.
%
%   THRESH takes the following flags at the end of the line of input
%   arguments:
%
%-     'hard'   - Perform hard thresholding. This is the default.
%
%-     'soft'   - Perform soft thresholding.  
%
%-     'full'   - Returns the output as a full matrix. This is the default.
%
%-     'sparse' - Returns the output as a sparse matrix.
%
%   The function WTHRESH in the Matlab Wavelet toolbox implements the same
%   functionality.
%
%   See also: largestr, largestn

%   AUTHOR : Peter Soendergaard and Bruno Torresani.  
%   TESTING: OK
%   REFERENCE: OK

if nargin<2
  error('Too few input parameters.');k
end;

if (prod(size(lambda))~=1 || ~isnumeric(lambda))
  error('lambda must be a scalar.');
end;

% Define initial value for flags and key/value pairs.
definput.flags.iofun={'hard','soft'};
definput.flags.outclass={'full','sparse'};

[flags,keyvals]=ltfatarghelper({},definput,varargin);

if flags.do_sparse
  if ndims(xi)>2
    error('Sparse output is only supported for 1D/2D input. This is a limitation of Matlab/Octave.');
  end;
end;

if flags.do_sparse
  xo=sparse(size(xi,1),size(xi,2));
    
  if flags.do_hard
    % Create a significance map pointing to the non-zero elements.
    signifmap=find(abs(xi)>=lambda);
    
    xo(signifmap)=xi(signifmap);
  else
    % Create a significance map pointing to the non-zero elements.
    signifmap=find(abs(xi)>lambda);
    
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
    
else
  xo=zeros(size(xi));
  
  % Create a mask with a value of 1 for non-zero element. For full
  % matrices, this is faster than the significance map.

  if flags.do_hard
    if nargout==2
      mask=abs(xi)>=lambda;
      N=sum(mask(:));
      xo=xi.*mask;
    else
      xo=xi.*(abs(xi)>=lambda);    
    end;
    
  else
    % In the following lines, the +0 is significant: It turns
    % -0 into +0, oh! the joy of numerics.
    
    if nargout==2
      xa=abs(xi)-lambda;    
      mask=xa>=0;
      xo=(mask.*xa+0).*sign(xi);
      N=sum(mask(:));
      
    else
      xa=abs(xi)-lambda;    
      xo=((xa>=0).*xa+0).*sign(xi);
    end;

  end;
end;

