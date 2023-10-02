function state = partconv_dualfdl_init( B_short, B_long, W, h)
%PARTCONV_DUALFDL_INIT Initialize dual FDL partitioned convolution
%   Usage: state = partconv_dualfdl_init( B_short, B_long, W, h)
%
%   Input parameters:
%       B_short : Initial buffer length
%       B_long  : Secondary buffer length
%       W       : Number of channels
%       h       : Impulse response
%
%   Output parameters:
%       state   : Convolution state
%
%   state = PARTCONV_DUALFDL_INIT( B_short, B_long, W, h) initializes a struct *state* for 
%   computing *W*-channel dual FDL partitioned convolution with impulse response *h*
%   divided into blocks of length *B_long* with the first block sub-divided into blocks of length *B_short*.
%   *B_long* must be divisible by *B-short*. 
%
%   The short FDL computes the convolution with the first *B_long* samples of *h*, the long one with the rest.
%
%   If *h* is a *hl x hw* matrix, its collumns are treated as *hw* impulse responses and
%   multiple convolutions are computed in parallel.

if 0 ~= rem(B_long,B_short)
    error('%s: Long buffer length must be divisible by short buffer length', mfilename);
end

[hl,hw] = size(h);
hlpad   = max( [ 2*B_long, ceil(hl/B_long)*B_long ] );
h       = postpad(h,hlpad); 

state.H_parts_long       = fftreal( reshape( h , B_long, hlpad/B_long, hw ), 2*B_long );
state.H_parts_long(:,1,:) = []; % The first block is replaced by shorter ones
state.FDL_long           = zeros([size(state.H_parts_long),W]);
state.B_long             = B_long;

% Double-buffering between the FDLs
state.FDL_buf_in         = zeros([B_long,W,2]);
state.FDL_buf_in_widx    = 1;
state.FDL_buf_in_ridx    = 2;
state.FDL_buf_out        = zeros([B_long,W,hw,2]);
state.FDL_buf_out_widx   = 1;
state.FDL_buf_out_ridx   = 2;
state.FDL_counter        = 0;

state.H_parts_short     = fftreal( reshape( postpad( h, B_long ), B_short, B_long/B_short, hw ), 2*B_short );
state.FDL_short         = zeros([size(state.H_parts_short),W]);
state.B_short           = B_short;

state.W                 = W;
state.prevblock_long    = zeros(B_long ,W);
state.prevblock_short   = zeros(B_short,W);

end
