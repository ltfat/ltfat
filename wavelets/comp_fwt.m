function c=comp_fwt(f,h,J)
%COMP_FWT   Computes Fast Wavelet Transform 
%   Usage:  c=comp_fwt(f,h,J);
%  
%   Input parameters:
%         f     : Input data.
%         h     : Wavelet filters.
%         J     : Number of filterbank iterations.
%   Output parameters:
%         c     : Coefficients stored in J+1 cell-array.
%
%   See also:
%
%   Demos:
%
%   References:

%   AUTHOR : Zdenek Prusa.
%   TESTING: COMP_FWT
%   REFERENCE: COMP_FWT

[L,W] = size(f);
c = cell(J+1,1);

flen = length(h{1});
% prepare impulse responses
hpad = cell(J,length(h));
 for hidx=1:length(h)
            for j=0:J-1
               hpad{j+1,hidx}=zeros(L/2^j,1);
               hpad{j+1,hidx}(1:flen) = h{hidx}(:);
               hpad{j+1,hidx} = circshift(hpad{j+1,hidx},-ceil(flen/2)+1);
            end
 end
 
c{1} = zeros(L/2^J,W);
for j=J:-1:1
    c{end-j+1} = zeros(L/2^j,W);
end


     for w=1:W    
         tempa = f(:,w);
         for j=1:J
             c{end+1-j}(:,w) = downs(pconv(tempa,hpad{j,2}));
             tempa = downs(pconv(tempa,hpad{j,1}));
         end
         c{1}(:,w) = tempa;
     end