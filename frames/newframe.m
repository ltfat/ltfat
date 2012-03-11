function F=newframe(ftype,varargin);
  
ftype=lower(ftype);
switch(ftype)
 case 'dgt'
  F.g=varargin{1};
  F.a=varargin{2};
  F.M=varargin{3};
 case 'dgtreal'
  F.g=varargin{1};
  F.a=varargin{2};
  F.M=varargin{3};
 case {'dwilt','wmdct'}
  F.g=varargin{1};
  F.M=varargin{2};  
end;

F.type=ftype;
F.gd={'dual',F.g};