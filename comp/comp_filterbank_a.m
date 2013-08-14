function [a,info]=comp_filterbank_a(a,M,info)
%COMP_FILTERBANK_A  Return sanitized a
%   Usage:  [a,info]=comp_filterbank_a(a,M);
%   
%   `[a,info]=comp_filterbank_a(a,M)` returns a sanitized version of *a*
%   expand to a $M\times 2$ matrix, and update the information in *info*.

% FIXME: Not sufficiently safe in the case where there is only one
% channel, or when attempting to specify a uniform filterbank with
% fractional downsampling.
%


info.isfractional=0;
info.isuniform=0;

if M==1 && size(a,1)~=1
   error('%s: One channel, but more a.',upper(mfilename));
end

if M==1 && size(a,2)<2
   a = [a,1];
end

if isvector(a) && M~=1
    [a,~]=scalardistribute(a(:),ones(M,1));
    
    if all(a==a(1))
        info.isuniform=1;
    end;

    a=[a,ones(M,1)];
else
    % We need to check against the case where this routine has already
    % been run
    if isequal(a(:,2),ones(M,1))
        if all(a(:,1)==a(1))
            info.isuniform=1;
        end;        
    else
        info.isfractional=1;
    end;
    
    % If the filterbank uses fractional downsampling, it cannot be
    % treated by the uniform algorithms, even though the sampling rate is uniform.
    
    % FIXME: Fractional, uniform filterbanks are not handled, they are
    % not allowed.

end;

info.a=a;

