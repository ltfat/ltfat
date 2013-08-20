function blockplot(p,F,c)
%BLOCKPLOT Plot block coefficients
%   Usage: blockplot(p,F,c);
%
%   Input parameters:
%         p     : JAVA object of the class net.sourceforge.ltfat.SpectFrame.
%         F     : Frame object.
%         c     : Block coefficients.
%
%   `blockplot(p,F,c)` appends the block coefficients `c` to the running 
%   coefficient plot in `p`.

if size(c,2)>1
   error('%s: Only one channel input is supported.',upper(mfilename));
end

ctf = framecoef2tf(F,c(:,1));

if strcmp(F.blockalg,'sliced')
   % DO the coefficient overlapping or cropping
end

ctf = cast(ctf,'single');
javaMethod('append',p,ctf);


