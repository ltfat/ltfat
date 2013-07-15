function [par,varargout] = blockpanelget(p,varargin)
%BLOCKPANELGET Get parameters from GUI
%   Usage: [par,...] = blockpanelget(p,spar,...)
%
%   Input parameters:
%         p    : JAVA object
%         spar : String. Name of the parameter.
%
%   Output parameters:
%         par  : Single value or vector of parameters.
%
%   `par = blockpanelget(p,'spar')` gets a current value of the parameter
%   `'spar'` from the GUI specified in `p`. 
%   `par = blockpanelget(p,'spar1','spar2','spar3')` gets current values
%   of the parameters `'spar1'`, `'spar2'` and `'spar3'` from the GUI and
%   stores them in a vector `p`.
%   `[par1,par2,par3] = blockpanelget(p,'spar1','spar2','spar3')` gets 
%   current values of the parameters `'spar1'`, `'spar2'` and `'spar3'` 
%   from the GUI and stores them in the separate variables `par1`,`par2`
%   and `par3`. 
%

nout = nargout();

if nargin<2
   error('%s: Too few input arguments.',upper(mfilename));
end

if nout~=1 && nout ~= numel(varargin)
   error('%s: Number of inputs does not match with number of outputs.',upper(mfilename));
end


res = javaMethod('getParams',p,varargin);

if nout == 1
   par = res;
else
   par = res(1);
   varargout = mat2cell(res(2:end),nout-1,1);
end

