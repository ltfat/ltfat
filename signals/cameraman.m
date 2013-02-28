function s=cameraman();
%CAMERAMAN  Load the 'cameraman' test image
%   Usage: s=cameraman;
% 
%   `cameraman` loads a $256 \times 256$ greyscale image of a cameraman.
% 
%   The returned matrix `s` consists of integers between 0 and 255,
%   which have been converted to double precision.
% 
%   To display the image, use `imagesc` with a gray colormap:::
% 
%     imagesc(cameraman); colormap(gray); axis('image');
% 
%   See `<ftp://nic.funet.fi/pub/graphics/misc/test-images/>`_ or
%   `<http://sipi.usc.edu/database/database.cgi?volume=misc>`_.

%   AUTHOR : Peter L. Søndergaard
%   TESTING: TEST_SIGNALS

  
if nargin>0
  error('This function does not take input arguments.')
end;

f=mfilename('fullpath');

s=double(imread([f,'.png']));
