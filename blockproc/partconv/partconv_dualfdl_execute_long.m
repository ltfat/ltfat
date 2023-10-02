function [y,state] = partconv_dualfdl_execute_long( x, state )
%PARTCONV_DUALFDL_EXECUTE_LONG Compute dual FDL partitioned convolution
%   Usage: [y, state] = partconv_dualfdl_execute_long( x, state )
%
%   Input parameters:
%       x       : Input signal
%       state   : Input convolution state
%
%   Output parameters:
%       y       : Output signal
%       state   : Output convolution state
%
%   [y,state] = PARTCONV_DUALFDL_EXECUTE_LONG( x, state ) computes dual-length FDL partitioned covolution. The size of
%   *x* must be *state.B_long * state.W*. This is a wrapper using the longer block size in order to avoid comp. spikes.
%   This erases the benefit of lower latency introduced by the shorter blocks.

%   AUTHOR : Zdenek Prusa (2023)

[B,W] = size(x);
[hw] = size(state.H_parts_long,3);

if B ~= state.B_long
    error('Wrong buffer length');
end

y = zeros([B,W,hw]);

for ii = 0:state.B_long/state.B_short - 1
    idx_slice = ( ii * state.B_short + 1 ) : ((ii+1) * state.B_short);
    [y(idx_slice,:,:),state] = partconv_dualfdl_execute( x(idx_slice,:), state );
end

end % function

