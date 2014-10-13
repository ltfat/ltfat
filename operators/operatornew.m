function Op=operatornew(otype,varargin)
%OPERATORNEW  Construct a new operator
%   Usage: F=operatornew(otype,...);
%
%   `Op=operatornew(otype,...)` constructs a new operator object *Op* of type
%   *otype*. Arguments following *otype* are specific to the type of operator
%   chosen.
%
%   Frame multipliers
%   -----------------
%
%   `operatornew('framemul',Fa,Fs,s)` constructs a frame multiplier with
%   analysis frame *Fa*, synthesis frame *Fs* and symbol *s*. See the help on
%   |framemul|.
%
%   Spreading operators
%   -------------------
%
%   `operatornew('spread',s)` constructs a spreading operator with symbol
%   *s*. See the help on |spreadop|.
%  
%   See also: operator, ioperator

  
if nargin<2
  error('%s: Too few input parameters.',upper(mfilename));
end;

if ~ischar(otype)
  error(['%s: First agument must be a string denoting the type of ' ...
         'frame.'],upper(mfilename));
end;

otype=lower(otype);

switch(otype)
  case 'framemul'
    Op.Fa=varargin{1};
    Op.Fs=varargin{2};
    Op.s =varargin{3};
    Op.L =framelengthcoef(Op.Fs,size(Op.s,1));
  case 'spread'
    Op.s =varargin{1};
    Op.L =length(Op.s);
  otherwise
    error('%s: Unknows operator type: %s',upper(mfilename),otype);  
end;

Op.type=otype;
