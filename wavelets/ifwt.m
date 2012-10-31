function f = ifwt(c,varargin)
%IFWT   Inverse Fast Wavelet Transform 
%   Usage:  f = ifwt(c,g)
%           f = ifwt(c,g,Ls,...)
%           f = ifwt(c,Lc,g,...)
%           f = ifwt(c,Lc,g,Ls,...)
%
%   Input parameters:
%         c     : Coefficients stored in J+1 cell-array or in packed format.
%         g     : Synthesis wavelet filters.
%         Ls    : Length of the reconstructed signal.
%         Lc    : Lengths of the wavelet coefficients in c.
%
%   The following flags are supported:
%
%         'dec','undec'
%                Type of the wavelet transform.
%
%         'per','zpd','sym','symw','asym','asymw','ppd','sp0'
%                Type of the boundary handling.
%
%   Output parameters:
%         f     : Reconstructed data.
%
%
%   See also:
%
%   Demos:
%
%   References:

%   AUTHOR : Zdenek Prusa.
%   TESTING: TEST_IFWT
%   REFERENCE: REF_IFWT

if(iscell(c))
    if nargin<2
      error('%s: Too few input parameters.',upper(mfilename));
    end;
else
    if nargin<3
      error('%s: Too few input parameters.',upper(mfilename));
    end; 
end
%% PARSE INPUT
definput.keyvals.Ls=[];    
definput.import = {'fwt'};

if(iscell(c))
    g = varargin{1};
    [flags,kv,Ls]=ltfatarghelper({'Ls'},definput,{varargin{2:end}});
elseif(isnumeric(c))
    Lc = varargin{1};
    g = varargin{2};
    c = pack2cell(c,Lc); % TO DO: adapt comp functions to work with the packed format
    [flags,kv,Ls]=ltfatarghelper({'Ls'},definput,{varargin{3:end}});
else
    error('%s: Unrecognized coefficient format.',upper(mfilename));
end
%% CHECK INPUT
do_definedfb = 0;
if(iscell(g))
    if( length(g)<2)
       error('%s: h is expected to be a cell array containing two wavelet filters.',upper(mfilename)); 
    end

    if(length(g{1})< 2 && length(g{2})< 2)
        error('%s: Wavelet filters should have at least two coefficients.',upper(mfilename)); 
    end

    if(length(g{1})~=length(g{2}))
        error('%s: Wavelet filters have to have equal length.',upper(mfilename));
    end
elseif(isstruct(g))
    do_definedfb = 1;
elseif(ischar(g))
    g = waveletfb(g);
    do_definedfb = 1;
else
   error('%s: Unrecognized Wavelet filters definition.',upper(mfilename)); 
end




if(do_definedfb)
    % if type ws not defined use 
    if(flags.do_type_null)
       flags.type = g.type; 
    end

    if(flags.do_ext_null)
       flags.ext = g.ext; 
    end
    
    g = g.g;
else
    % manually setting defaults
    if(flags.do_type_null)
       flags.type = 'dec'; 
    end

    if(flags.do_ext_null)
       flags.ext = 'per'; 
    end
end



% Determine J from number of elements of c
J = length(c)-1;

if isempty(Ls)
    % Estimate output signal length from the number of coefficients
   [sigHalfLen,W] = size(c{end});
   % for 'undec' the Ls can be estimated exactly
   if(strcmp(flags.type,'undec'))
       if(strcmp(flags.ext,'per'))
           Ls = sigHalfLen;  
       else
           Ls = sigHalfLen - (length(g{1})-1);
       end
   elseif(strcmp(flags.type,'dec'))
       if(strcmp(flags.ext,'per'))
            % estimated Ls can be one sample more, if the original input
            % signal length was odd
           Ls = 2*sigHalfLen; 
       else
            % estimated Ls can be one sample more, if the original input
            % signal length plus length(h{1})-1 was an even number
           Ls = 2*sigHalfLen - (length(g{1})-2);
       end
   end
end


f = comp_ifwt_all(c,g,J,Ls,flags.type,flags.ext);


