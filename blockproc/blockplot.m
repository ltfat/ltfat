function cola=blockplot(p,F,c,cola)
%BLOCKPLOT Plot block coefficients
%   Usage: blockplot(p,F,c);
%
%   Input parameters:
%         p     : JAVA object of the class net.sourceforge.ltfat.SpectFrame.
%         F     : Frame object.
%         c     : Block coefficients.
%         cola  : (Optional) overlap from previous block.
%
%   Output parameters:
%         cola  : Overlap to the next block.
%
%   `blockplot(p,F,c)` appends the block coefficients `c` to the running 
%   coefficient plot in `p`.
%
%   `cola=blockplot(p,F,c,cola)` doest the same, but adds `cola` to the 
%   first respective coefficients in `c` and returns last coefficients from
%   `c`. This is only relevant for the sliced window blocking approach.

if size(c,2)>1
   error('%s: Only one channel input is supported.',upper(mfilename));
end

ctf = framecoef2tf(F,c(:,1));

if strcmp(F.blockalg,'sliced')
   % DO the coefficient overlapping or cropping
   %ctf = ctf(:,floor(end*3/8):floor(end*5/8)+1);
   
   if nargin>3 
      olLen = ceil(size(ctf,2)/2);
      if isempty(cola)
         cola = zeros(size(ctf,1),olLen,class(ctf));
      end
         
      ctf(:,1:olLen) = ctf(:,1:olLen) + cola;
      cola = ctf(:,end+1-olLen:end);
      ctf = ctf(:,1:olLen);
   end
end


ctf = abs(ctf);

if isoctave
   ctf = cast(ctf,'double');
   javaMethod('append',p,ctf(:),size(ctf,1),size(ctf,2));
else
   ctf = cast(ctf,'single');
   javaMethod('append',p,ctf);
end


