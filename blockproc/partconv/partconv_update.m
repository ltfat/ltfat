function state = partconv_update( state, h)
%PARTCONV_INIT Initialize FDL partitioned convolution
%   Usage: state = partconv_init( B, W, h )
%
%   Input parameters:
%       h   : Impulse response
%
%   state = PARTCONV_INIT( B, W, h) initializes a state struct *state* for computing
%   *W*-channel partitioned convolution with impulse response *h* using frequency delay line 
%   with block length *B*.
%   If *h* is a *hl x hw* matrix, its collumns are treated as *hw* impulse responses and
%   multiple convolutions are computed in parallel.
%   

[hl,hw] = size(h);
hlpad   = ceil(hl/state.B)*state.B ;

state.H_parts = fftreal(reshape( postpad( h, hlpad ), state.B, hlpad/state.B, hw ),2*state.B);

end % function


