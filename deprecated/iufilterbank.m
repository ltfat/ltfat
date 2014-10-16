function f=iufilterbank(varargin);  
%IUFILTERBANK  Filter bank inversion, DEPRECATED
%   Usage:  f=iufilterbank(c,g,a);
%
%   `iufilterbank` has been deprecated by |ifilterbank|. Call |ifilterbank|
%   with the exact same parameters as the old call to `iufilterbank`.
%
%   See also: ifilterbank

warning(['LTFAT: IUFILTERBANK has been deprecated, please use IFILTERBANK ' ...
         'instead.']);
  

f=ifilterbank(varargin{:});

