function state = partconv_update( state, h)
%PARTCONV_UPDATE Update FDL partitioned convolution
%   Usage: state = partconv_update( B, W, h )
%
%   Input parameters:
%       h   : Impulse response
%
%   state = PARTCONV_UPDATE( B, W, h) updates a state struct *state* for computing
%   *W*-channel partitioned convolution with impulse response *h* using frequency delay line 
%   with block length *B*.
%   If *h* is a *hl x hw* matrix, its collumns are treated as *hw* impulse responses and
%   multiple convolutions are computed in parallel.
%   

%   AUTHOR : Zdenek Prusa (2023)

[hl,hw] = size(h);
hlpad   = ceil(hl/state.B)*state.B ;

state.H_parts = fftreal(reshape( postpad( h, hlpad ), state.B, hlpad/state.B, hw ),2*state.B);

end % function


