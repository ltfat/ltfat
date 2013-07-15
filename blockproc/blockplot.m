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



ctf = framecoef2tf(F,c);

javaMethod('append',p,ctf);


