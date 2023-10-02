function state = partconv_init( B, W, h)
%PARTCONV_INIT Initialize FDL partitioned convolution
%   Usage: state = partconv_init( B, W, h )
%
%   Input parameters:
%       B   : Block length
%       W   : Number of channels
%       h   : Impulse response
%
%   state = PARTCONV_INIT( B, W, h) initializes a state struct *state* for computing
%   *W*-channel partitioned convolution with impulse response *h* using frequency delay line 
%   with block length *B*.
%   If *h* is a *hl x hw* matrix, its collumns are treated as *hw* impulse responses and
%   multiple convolutions are computed in parallel.
%   

%   AUTHOR : Zdenek Prusa (2023)

[hl,hw] = size(h);
hlpad   = ceil(hl/B)*B ;

state.H_parts       = fftreal(reshape( postpad( h, hlpad ), B, hlpad/B, hw ),2*B);
state.FDL           = zeros([size(state.H_parts),W]);
state.B             = B;
state.W             = W;
state.prevblock     = zeros(B,W);

end % function


