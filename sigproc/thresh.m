function [xo,N]=thresh(xi,lambda,varargin);
%THRESH   Coefficient thresholding
%   Usage:  x=thresh(x,lambda,...);
%           [x,N]=thresh(x,lambda,...);
%
%   `thresh(x,lambda)` will perform hard thresholding on *x*, i.e. all
%   elements with absolute value less than *lambda* will be set to zero.
%
%   `thresh(x,lambda,'soft')` will perform soft thresholding on *x*,
%   i.e. *lambda* will be subtracted from the absolute value of every element
%   of *x*.
%
%   `[x,N]=thresh(x,lambda)` additionally returns a number *N* specifying
%   how many numbers where kept.
%
%   `thresh` takes the following flags at the end of the line of input
%   arguments:
%
%     'hard'    Perform hard thresholding. This is the default.
%
%     'wiener'  Perform empirical Wiener shrinkage. This is in between
%               soft and hard thresholding.
%
%     'soft'    Perform soft thresholding.  
%
%     'full'    Returns the output as a full matrix. This is the default.
%
%     'sparse'  Returns the output as a sparse matrix.
%
%   The function `wthresh` in the Matlab Wavelet toolbox implements some of
%   the same functionality.
%
%   The following code produces a plot to demonstrate the difference
%   between hard and soft thresholding for a simple linear input:::
%
%     t=linspace(-4,4,100);
%     plot(t,thresh(t,1,'soft'),'r',...
%          t,thresh(t,1,'hard'),'.b',...
%          t,thresh(t,1,'wiener'),'--g');
%     legend('Soft thresh.','Hard thresh.','Wiener thresh.','Location','NorthWest');
%
%   See also: largestr, largestn
%
%   References: lim1979enhancement ghael1997improved

%   AUTHOR : Kai Siedenburg, Bruno Torresani and Peter L. Søndergaard.
%   TESTING: OK
%   REFERENCE: OK

if nargin<2
  error('Too few input parameters.');k
end;

if (prod(size(lambda))~=1 || ~isnumeric(lambda))
  error('lambda must be a scalar.');
end;

% Define initial value for flags and key/value pairs.
definput.import={'thresh'};
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
  end;
  
  if flags.do_wiener
    signifmap=find(abs(xi)>lambda);
    xo(signifmap) = 1 - (lambda./abs(xi(signifmap))).^2;
    xo(signifmap) = xi(signifmap).*xo(signifmap); 
  end;
  
  if flags.do_soft
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
  
  % Create a mask with a value of 1 for non-zero elements. For full
  % matrices, this is faster than the significance map.

  if flags.do_hard
    if nargout==2
      mask=abs(xi)>=lambda;
      N=sum(mask(:));
      xo=xi.*mask;
    else
      xo=xi.*(abs(xi)>=lambda);    
    end;
    
  end;
  
  if flags.do_soft
      % In the following lines, the +0 is significant: It turns
      % -0 into +0, oh! the joy of numerics.
      
      if nargout==2
          xa=abs(xi)-lambda;    
          mask=xa>=0;
          xo=(mask.*xa+0).*sign(xi);
          N=sum(mask(:))-sum(xa(:)==0);      
      else
          xa=abs(xi)-lambda;    
          xo=((xa>=0).*xa+0).*sign(xi);
      end;
  end;
  
  if flags.do_wiener
      xa = lambda./abs(xi);
      xa(isinf(xa)) = 0;
      xa = 1 - xa.^2;
      
      if nargout==2
          mask = xa>0;
          xo = xi.*xa.*mask;
          N = sum(mask(:));
      else
          xo = xi.*xa.*(xa>0);
      end
      
  end;
end;

