function varargout=scalardistribute(varargin)
%SCALARDISTRIBUTE  Copy scalar to array shape for parameter handling
%   Usage:  [...] = scalardistribute(...);
%
%   `[...]=scalardistribute(...)` copies the input parameters to the
%   output parameters.
%
%   * If one of the input parameters is an array, all the output parameters
%     will be column vectors containing the same number of elements. If one
%     of the other input parameters is a scalar, it will be replicated to
%     the correct length. This allows a scalar value to be repeated for
%     all conditions.
%
%   * If two or more input parameters are arrays, the must have the exact
%     same size. They will be converted to vectors and returned in the
%     output parameters. This allows two arrays to co-vary at the same time.
%
%   This operator is usefull for sanitizing input parameters: The user is
%   allowed to enter scalars or arrays as input paremeters. These input
%   are in turn passed to `scalardistribute`, which makes sure that the
%   arrays have the same shape, and that scalars are replicated. The user
%   of `scalardistibute` can now generate conditions based on all the
%   parameters, and be sure the have the right sizes.
%
%   As an example, consider:::
%
%     [a,b,c]=scalardistribute(1,[2,3],[4,5])
%
  
  npars=length(varargin);
  
  scalartest=zeros(npars,1);
  for ii=1:npars
    if ~isnumeric(varargin{ii})
      error('%s: Input no. %i must be numerical.',upper(mfilename),ii);
    end;
    if isempty(varargin{ii})
      error('%s: Input no. %i is empty.',upper(mfilename),ii);
    end;

    scalartest(ii)=isscalar(varargin{ii});
  end;

  idx=find(scalartest==0);
  
  % Tjeck that all non-scalars have the same shape
  for jj=1:numel(idx)
    ref_size  = size(varargin{idx(1)});
    this_size = size(varargin{idx(jj)});
    if any((ref_size-this_size)~=0)
      error('%s: Input no. %i and no. %i must have the same shape.',upper(mfilename),idx(1),idx(jj));
    end;
  end;  

  if numel(idx)==0
    % All arguments are scalar, so this is just a dummy, but it must
    % still be defined for the code not to fail.
    shape=1;
  else
    shape=ones(numel(varargin{idx(1)}),1);
  end;
  
  varargout=cell(1,npars);
  for ii=1:npars
    if scalartest(ii)
      % Replicate scalar
      varargout{ii}=shape*varargin{ii};
    else
      % Copy input and turn it into a vector. This would be a one-liner
      % in Octave.
      tmp=varargin{ii};
      varargout{ii}=tmp(:);
    end;
  end;
  
  
