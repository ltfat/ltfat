function f = ifwt(c,g,J,varargin)
%IFWT   Inverse Fast Wavelet Transform 
%   Usage:  f = ifwt(c,g,J)
%
%   Input parameters:
%         c     : Coefficients stored in J+1 cell-array.
%         h     : Reconstruction wavelet filters.
%         J     : Number of filterbank iterations.
%   Output parameters:
%         f     : Output data.
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


% works with coefficients obtained from input signals of lengts a*2^J
definput.flags.type = {'dec','undec'};
[flags,kv]=ltfatarghelper({},definput,varargin);

if(J<1)
   error('%s: J must be a positive integer.',upper(callfun)); 
end


gpad = cell(J,length(g));
[sigHalfLen,W] = size(c{end});  
sigLen = ceil((sigHalfLen*2)/2^J)*2^J;
flen = length(g{1});

       for hidx=1:length(g)
            for j=0:J-1
               gpad{j+1,hidx}=zeros(2*length(c{j+2}(:,1)),1);
               gpad{j+1,hidx}(1:flen) = g{hidx}(:);
               gpad{j+1,hidx} = circshift(gpad{j+1,hidx},-ceil(flen/2));
           end
        end

f = zeros(sigLen,W);
  
  
for w=1:W  
  tempa = c{1}(:,w);
    for j=1:J
       tempa = pconv(ups(tempa,2,2), gpad{j,1}) + pconv(ups(tempa,2,2), gpad{j,2});
    end
    f(1:length(tempa),w) = tempa;
end