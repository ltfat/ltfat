function red=framered(F);
  
switch(F.type)
 case 'dgt'
  red=F.M/F.a;
 case 'dgtreal'
  red=F.M/F.a;
 case {'dwilt','wmdct'}
  red=1;
end;

  