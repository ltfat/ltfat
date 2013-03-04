function f = comp_iufwt(c,g,J,a)
%COMP_IFWT_ALL Compute Inverse DWT
%   Usage:  f = comp_ifwt_all(c,g,J,Ls,type,ext);
%
%   Input parameters:
%         c     : Coefficients stored in a cell-array.
%         g     : Synthesis wavelet filters.
%
%   Output parameters:
%         f     : Reconstructed data.
%

% see comp_fwt for explanantion
assert(a(1)==a(2),'First two elements of a are not equal. Such wavelet filterbank is not suported.');

% Number of Filters
filtNo = numel(g);
% Determine number of output channels
if(length(size(c))>2)
   W = size(c,3);
else
   W = 1; 
end
% Determine output length
Ls = size(c,1);
% Allocate output data
f = zeros(Ls,W);
upFac= 1;
% Since no downsampling takes place, normalize impulse responses
tmpg = cell(filtNo,1);
for ii = 1:filtNo
   tmpg{ii} = g{ii}.h/sqrt(a(ii));
end

for ch=1:W
  tempca = c(:,1,ch);
  runPtr = 2;
  for jj=1:J
     filtUps = a(1)^(J-jj); 
     skip = floor(filtUps*g{1}.d) - filtUps; 
     tempca = comp_upconv({tempca},Ls,{tmpg{1}},upFac,skip,'per',filtUps); 
     for ff=2:filtNo
         skip = floor(filtUps*g{ff}.d) - filtUps; 
         tempca = tempca + comp_upconv({c(:,runPtr,ch)},Ls,{tmpg{ff}},upFac,skip,'per',filtUps);
         runPtr = runPtr + 1;
     end
  end
  f(:,ch) = tempca;
end

    
    