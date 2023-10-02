%this is a very basic demo for partitioned convolution, for helping me 
%understand how it works
clear all;

addpath([ltfatbasepath, '/blockproc/partconv']) %this is temporary


%% 1.) basic partitioned convolution with two different IRs

%pick one second of a test signal
[f,fs] = gspi;
ftmp = f(1:fs,:);
L = length(ftmp);

%choose the IR at the same length as the test signal
hl = fs;
%and choose an initial block size
B = 2048;
%channel number of the input signal
W = size(ftmp,2);

%TODO: load an IR via SOFAload
h = ones(hl,2);
h = bsxfun(@times,h,1:2);

%bring the input signal to a suitable length
Lpad = ceil((L + hl)/B)*B;
f = postpad(ftmp, Lpad);

%initialize
state = partconv_init( B, W, h);

%prepare the block matrix for each channel of the input signal
f_blocks = zeros( B, W, Lpad/B);
for w = 1:W
    f_blocks(:,w,:) = reshape( f(:,w), B, Lpad/B );
end

fout_blocks = zeros( B, W, size(h,2), Lpad/B );
block_count = Lpad/B;

for b = 1:block_count
    [ fout_blocks(:,:,:,b), state ] =  partconv_execute( f_blocks(:,:,b), state );
end

fout = reshape(permute(fout_blocks,[1,4,2,3]),Lpad,W,size(h,2));

%single channel signal convolved with two different IRs
figure
plot(squeeze(fout(:,:,1)))
hold on
plot(squeeze(fout(:,:,2)))


%% 2.) basic partitioned convolution with a stereo signal

%pick one second of a test signal
[f,fs] = gspi;
ftmp = f(1:fs,:);
L = length(ftmp);
%make the signal 2 channels
ftmp = repmat(ftmp, 1, 2);
ftmp(:,2) = 0.5*ftmp(:,2);

%choose the IR at the same length as the test signal
hl = fs;
%and choose an initial block size
B = 2048;
%channel number of the input signal
W = size(ftmp,2);

%TODO: load an IR via SOFAload
h = ones(hl,1);

%bring the input signal to a suitable length
Lpad = ceil((L + hl)/B)*B;
f = postpad(ftmp, Lpad);

%initialize
state = partconv_init( B, W, h);

%prepare the block matrix for each channel of the input signal
f_blocks = zeros( B, W, Lpad/B);
for w = 1:W
    f_blocks(:,w,:) = reshape( f(:,w), B, Lpad/B );
end

fout_blocks = zeros( B, W, size(h,2), Lpad/B );
block_count = Lpad/B;

for b = 1:block_count
    [ fout_blocks(:,:,:,b), state ] =  partconv_execute( f_blocks(:,:,b), state );
end

fout = reshape(permute(fout_blocks,[1,4,2,3]),Lpad,W,size(h,2));

%the stereo signal can be plotted immediately
figure
plot(fout)


%% 3.) partitioned convolution with two different block sizes on a dual
% frequency delay line

%pick one second of a test signal
[f,fs] = gspi;
ftmp = f(1:fs,:);
L = length(ftmp);
%make the signal 2 channels
ftmp = repmat(ftmp, 1, 2);
ftmp(:,2) = 0.5*ftmp(:,2);

%choose the IR at the same length as the test signal
hl = fs;

%now, there are two different block sizes
B_long  = 8192;
B_short = B_long/8;
%channel number of the input signal
W = size(ftmp,2);

%TODO: load an IR via SOFAload
h = ones(hl,1);

Lpad = ceil((L + hl)/B_long)*B_long;
f = postpad(ftmp,Lpad);

%the state has now many more parameters...
state = partconv_dualfdl_init( B_short, B_long, W, h);

%again, prepare the blocks (according to the shortest "unit of processing")
f_blocks = zeros( B_short, W, Lpad/B_short);
for w = 1:W
    f_blocks(:,w,:) = reshape( f(:,w), B_short, Lpad/B_short );
end

fout_blocks = zeros( B_short, W, size(h,2), Lpad/B_short );
block_count = Lpad/B_short;
for b = 1:block_count
    [ fout_blocks(:,:,:,b), state ] =  partconv_dualfdl_execute( f_blocks(:,:,b), state );
end

fout = reshape(permute(fout_blocks,[1,4,2,3]),Lpad,W,size(h,2));

%the stereo signal can be plotted immediately
figure
plot(fout)