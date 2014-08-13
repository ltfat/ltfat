function c = block_fwt( f, w, J)
%BLOCK_FWT FWT func. wrapper for a block processing
%   Usage: c = block_fwt( f, w, J);
%
%   Input parameters:
%         f     : Input data.
%         w     : Analysis Wavelet Filterbank. 
%         J     : Number of filterbank iterations.
%
%   Output parameters:
%         c      : Coefficient vector.
%
%   `c = block_fwt(f,w,J)` accepts suitably extended block of data *f*
%   and produces correct coefficients using the SegDWT algorithm (based on
%   overlap-save block convolution) with wavelet filters defined by *w* 
%   and *J* levels. *f* is expected to be a column vector or a matrix and 
%   the processing is done column-wise.
%
%   Do not call this function directly. The function is called from 
%   |blockana| when used with frame type 'fwt' and 'segola' block transform
%   handling see |blockframeaccel|.
%
%   Function should be independent of block_interface.
%
%   See also: block, block_ifwt
%
%   References: ltfatnote026

if nargin<3
  error('%s: Too few input parameters.',upper(mfilename));
end;

% Initialize the wavelet filters structure
%h = fwtinit(h,'ana');

if any(w.a~=w.a(1))
   error('%s: Non-equal subsampling factors are not supported.',upper(mfilename));
end

w = fwtinit(w);
% Extended block length 
Ls = size(f,1);
% Low-pass filter length
m = numel(w.h{1}.h);
% Low-pass subsampling factor
a = w.a(1);
% Extension length
rred = (a^J-1)/(a-1)*(m-a);
% Block boundaries
blocksize=w.a(1)^J;
% Input signal samples to be processed

% This is effectivelly the "negative" right extension described in chapter
% 4.1.4 in the reference.
L=rred+floor((Ls-rred)/blocksize)*blocksize;

levelLen = L;
filtNo = length(w.h);
subbNo = (filtNo-1)*J+1;
Lc = zeros(subbNo,1);
runPtr = 0; 
for jj=1:J
   for ff=filtNo:-1:2
      Lc(end-runPtr) = floor((levelLen-m-1)/w.a(ff));
      runPtr = runPtr + 1;
   end
   levelLen = floor((levelLen-m-1)/w.a(1));
end
Lc(1)=levelLen; 

% 
%[Lc, L] = fwtclength(Ls,h,J,'valid');

% Crop to the right length
if(Ls>L)
   f=postpad(f,L); 
end

if Ls<rred+a^J
   error('%s: Insufficient input signal length for the %s flag. Minimum is %i.',upper(mfilename),'''valid''',rred+a^J);
end

c = comp_fwt(f,w.h,w.a,J,'valid');

% Do the cropping 
runPtr = 0; 
for jj=1:J-1
   for ff=filtNo:-1:2
      cstart = (a^(J-jj)-1)/(a-1)*(m-a);
      c{end-runPtr} = c{end-runPtr}(cstart+1:end,:);
      runPtr = runPtr + 1;
   end
end

% To the pack format
c = cell2mat(c);
