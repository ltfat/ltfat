function F=frameaccel(F,L);  
  
switch(F.type)
 case {'dgt','dgtreal'}
  [F.g,F.g_info]  =gabwin(F.g,F.a,F.M,L);
  [F.gd,F.gd_info]=gabwin(F.g,F.a,F.M,L);
 case {'dwilt','wmdct'}
  [F.g,F.g_info]  =gabwin(F.gd,F.M,L);
  [F.gd,F.gd_info]=gabwin(F.gd,F.M,L);
end;

  