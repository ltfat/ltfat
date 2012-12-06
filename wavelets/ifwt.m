function f = ifwt(c,g,J,varargin)
%IFWT   Inverse Fast Wavelet Transform 
%   Usage:  f = ifwt(c,g)
%           f = ifwt(c,g,Ls,...)
%
%   Input parameters:
%         c     : Coefficients stored in J+1 cell-array.
%         g     : Synthesis wavelet filters.
%         Ls    : Length of the reconstructed signal.
%
%   Output parameters:
%         f     : Reconstructed data.
%
%   The following flags are supported:
%
%         'dec','undec'
%                Type of the wavelet transform.
%
%         'per','zpd','sym','symw','asym','asymw','ppd','sp0'
%                Type of the boundary handling.
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

    if nargin<2
      error('%s: Too few input parameters.',upper(mfilename));
    end;

%% PARSE INPUT
definput.keyvals.Ls=[];    
definput.import = {'fwt'};

if(iscell(c))
    [flags,kv,Ls]=ltfatarghelper({'Ls'},definput,varargin);
else
    error('%s: Unrecognized coefficient format.',upper(mfilename));
end
%% CHECK INPUT
do_definedfb = 0;
if(iscell(g))
    if( length(g)<2)
       error('%s: h is expected to be a cell array containing two or more wavelet filters.',upper(mfilename)); 
    end


    for ii=2:numel(g)
       if(length(g{1})~=length(g{ii}))
           error('%s: Wavelet filters have to have equal length.',upper(mfilename));
       end
    end
    
    if(length(g{1})< 2)
        error('%s: Wavelet filters should have at least two coefficients.',upper(mfilename)); 
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
    a = g.a;
    g = g.g;
else
    % manually setting defaults
    if(flags.do_type_null)
       flags.type = 'dec'; 
    end

    if(flags.do_ext_null)
       flags.ext = 'per'; 
    end
    % critical subsampling by default 
    a = length(g)*ones(length(g),1);
end



% TO REMOVE: Determine J from number of elements of c
[cR cC] = size(c);


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


f = comp_ifwt_all(c,g,J,a,Ls,flags.type,flags.ext);


