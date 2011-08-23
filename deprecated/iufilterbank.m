function f=iufilterbank(varargin);  
%IUFILTERBANK  Filter bank inversion, DEPRECATED
%   Usage:  f=iufilterbank(c,g,a);
%
%   IUFILTERBANK has been deprecated by IFILTERANK. Call IFILTERBANK with
%   the exact same parameters as the old call to IUFILTERBANK.
%
%   See also: ifilterbank

warning(['LTFAT: IUFILTERBANK has been deprecated, used IFILTERBANK ' ...
         'instead.']);
  

f=ifilterbank(varargin{:});

