function RGB = jpeg2rgb(YCbCr)
%JPEG2RGB  Coverts from RGB format to YCbCr format
%   Usage:  RGB = jpeg2rgb(YCbCr);
% 
%   Input parameters:
%         YCbCr : 3d data-cube, containing the YCbCr information of the
%                 image
% 
%   Output parameters:
%         RGB   : 3d data-cube, containing RGB information of the image
% 
%   'jpeg2rgb(YCbCr)' performs a transformation of the 3d data-cube *YCbCr* with
%   dimensions $N \times M \times 3$, which contains information of
%   "luminance", "chrominance blue" and "chrominance red".  The output
%   variable *RGB* is a 3d data-cube of the same size containing information
%   about the colours "red", "green" and "blue". The output will be of
%   the `uint8` type.
% 
%   For more information, see
%   `<http://en.wikipedia.org/wiki/YCbCr>`_ and `<http://de.wikipedia.org/wiki/JPEG>`_
%
%   See also:   rgb2jpeg

% AUTHOR:   Markus Faulhuber, February 2013

[s1,s2,s3] = size(YCbCr);
YCbCr = double(YCbCr);

if s3 ~= 3
    disp('Sorry, this routine is for YCbCr of dimension NxMx3 only')
    return;
end

RGB(:,:,1) = YCbCr(:,:,1)+1.402*(YCbCr(:,:,3)-128);
RGB(:,:,2) = YCbCr(:,:,1)-0.3441*(YCbCr(:,:,2)-128)-0.7141*(YCbCr(:,:,3)-128);
RGB(:,:,3) = YCbCr(:,:,1)+1.772*(YCbCr(:,:,2)-128);

RGB = uint8(RGB);
