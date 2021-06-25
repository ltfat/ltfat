function [g,info] = gabwin(g,a,M,varargin);
%GABWIN  Compute a Gabor window from text or cell array
%   Usage: [g,info] = gabwin(g,a,M,L);
%
%   `[g,info]=gabwin(g,a,M,L)` computes a window that fits well with the
%   specified number of channels *M*, time shift *a* and transform length
%   *L*. The window itself is specified by a text description or a cell array
%   containing additional parameters.
%
%   The window can be specified directly as a vector of numerical
%   values. In this case, `gabwin` only checks assumptions about transform
%   sizes etc.
%
%   `[g,info]=gabwin(g,a,M)` does the same, but the window must be a FIR
%   window, as the transform length is unspecified.
%
%   `gabwin(g,a,M,L,lt)` or `gabwin(g,a,M,[],lt)` does as above but for a
%   non-separable lattice specified by *lt*. Please see the help of
%   |matrix2latticetype| for a precise description of the parameter *lt*.
%
%   The window can be specified as one of the following text strings:
%  
%     'gauss'      Gaussian window fitted to the lattice,
%                  i.e. $tfr=a\cdot M/L$.
%
%     'dualgauss'  Canonical dual of Gaussian window.
%
%     'tight'      Tight window generated from a Gaussian.
%
%   In these cases, a long window is generated with a length of *L*.
%
%   It is also possible to specify one of the window names from |firwin|. In
%   such a case, `gabwin` will generate the specified FIR window with a length
%   of *M*.
%
%   The window can also be specified as cell array. The possibilities are:
%
%     `{'gauss',...}`
%         Additional parameters are passed to |pgauss|. When no additional
%         parameters are passed, the window is generated according to the
%         defaults in |pgauss|.
%
%     `{'dual',...}`
%         Canonical dual window of whatever follows. See the examples below.
%
%     `{'tight',...}` 
%         Canonical tight window of whatever follows.
%
%   It is also possible to specify one of the window names from |firwin| as
%   the first field in the cell array. In this case, the remaining
%   entries of the cell array are passed directly to |firwin|.
%
%   Some examples: To compute a Gaussian window of length *L* fitted for a
%   system with time-shift *a* and *M* channels use::
%
%     g=gabwin('gauss',a,M,L);
%
%   To compute Gaussian window with equal time and frequency support
%   irrespective of *a* and *M*::
%
%     g=gabwin({'gauss'},a,M,L);
%
%   To compute the canonical dual of a Gaussian window fitted for a
%   system with time-shift *a* and *M* channels::
%
%     gd=gabwin('gaussdual',a,M,L);
%
%   To compute the canonical tight window of the Gaussian window fitted
%   for the system::
%
%     gd=gabwin({'tight','gauss'},a,M,L);
%
%   To compute the dual of a Hann window of length 20::  
% 
%     g=gabwin({'dual',{'hann',20}},a,M,L);
%
%   The structure *info* provides some information about the computed
%   window:
%
%     `info.gauss`
%        True if the window is a Gaussian.
%
%     `info.tfr`
%        Time/frequency support ratio of the window. Set whenever it makes sense.
%
%     `info.wasrow`
%        Input was a row window
%
%     `info.isfir`
%        Input is an FIR window
%
%     `info.isdual`
%        Output is the dual window of the auxiliary window.
%
%     `info.istight`
%        Output is known to be a tight window.
%
%     `info.auxinfo`
%        Info about auxiliary window.
%   
%     `info.gl`
%        Length of window.
%
%   See also: pgauss, firwin, wilwin
  
% Assert correct input.
if nargin<3
  error('%s: Too few input parameters.',upper(mfilename));
end;

definput.keyvals.L=[];
definput.keyvals.lt=[0 1];
definput.keyvals.callfun='GABWIN';
definput.flags.phase={'freqinv','timeinv'};
[~,kv,L,lt]=ltfatarghelper({'L','lt'},definput,varargin,'gabwin');

[g,info] = comp_window(g,a,M,L,lt,kv.callfun);

