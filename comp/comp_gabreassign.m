function sr=comp_gabreassign(s,tgrad,fgrad,a);
%COMP_GABREASSIGN  Reassign time-frequency distribution.
%   Usage:  sr = comp_gabreassign(s,tgrad,fgrad,a);
%
%   COMP_GABREASSIGN(s,tgrad,fgrad,a) will reassign the values of the positive
%   time-frequency distribution s using the instantaneous time and frequency
%   fgrad and ifdummy. The lattice is determined by the time shift _a and
%   the number of channels deduced from the size of s.
%
%   See also: gabreassign
%
%R  aufl95

%   AUTHOR : Peter Soendergaard.
%   TESTING: OK
%   REFERENCE: OK
  
[M,N,W]=size(s);
L=N*a;
b=L/M;

freqpos=fftindex(M);  
for w=1:W
  tgrad(:,:,w)=tgrad(:,:,w)/b+repmat(freqpos,1,N);
end;

timepos=fftindex(N);
for w=1:W
  fgrad(:,:,w)=fgrad(:,:,w)/a+repmat(timepos.',M,1);
end;

tgrad=round(tgrad);
fgrad=round(fgrad);

tgrad=mod(tgrad,M);
fgrad=mod(fgrad,N);  
  
sr=zeros(M,N,W);

fgrad=fgrad+1;
tgrad=tgrad+1;

for ii=1:M
  for jj=1:N      
    sr(tgrad(ii,jj),fgrad(ii,jj)) = sr(tgrad(ii,jj),fgrad(ii,jj))+s(ii,jj);
  end;
end;  


%OLDFORMAT
