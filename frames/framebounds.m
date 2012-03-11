function varargout=framebounds(F);
  
switch(ftype)
 case {'dgt','dgtreal'}
  varargout=gabframebounds(F.g,F.a,F.M); 
 case {'dwilt','wmdct'}
  varargout=wilbounds(F.g,F.M); 
end;

  