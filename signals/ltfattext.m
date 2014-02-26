function s=ltfattext();
%LTFATTEXT  Load the 'ltfattext' test image
%   Usage: s=ltfattext;
% 
%   `ltfattext` loads a $401 \times 600$ black and white image of the word
%   'LTFAT'.
% 
%   The image is assumed to be used as a spectrogram with 800 channels
%   as produced by |dgtreal|.
% 
%   The returned matrix `s` consists of the integers 0 and 1, which have
%   been converted to double precision.
% 
%   To display the image, use `imagesc` with a gray colormap:::
% 
%     imagesc(ltfattext);
%     colormap(gray);
%     axis('xy');
% 
%   See also: ltfatlogo, dgtreal


%   AUTHOR : Peter L. Søndergaard
%   TESTING: TEST_SIGNALS
%   REFERENCE: OK
  
if nargin>0
  error('This function does not take input arguments.')
end;

f=mfilename('fullpath');

s=flipud(double(imread([f,'.png'])))/255;
