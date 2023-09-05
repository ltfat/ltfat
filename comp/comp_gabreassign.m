function sr=comp_gabreassign(s,tgrad,fgrad,a);
%COMP_GABREASSIGN  Reassign time-frequency distribution.
%   Usage:  sr = comp_gabreassign(s,tgrad,fgrad,a);
%
%   `comp_gabreassign(s,tgrad,fgrad,a)` will reassign the values of the positive
%   time-frequency distribution *s* using the instantaneous time and frequency
%   *fgrad* and *tgrad*. The lattice is determined by the time shift *a* and
%   the number of channels deduced from the size of *s*.
%
%   See also: gabreassign
%
%   References: aufl95

%   AUTHOR : Peter L. Søndergaard.
%   TESTING: OK
%   REFERENCE: OK

[M,N,W]=size(s);
L=N*a;
b=L/M;

freqpos=fftindex(M);  
tgrad=bsxfun(@plus,tgrad/b,freqpos);

timepos=fftindex(N);
fgrad=bsxfun(@plus,fgrad/a,timepos.');

tgrad=round(tgrad);
fgrad=round(fgrad);

tgrad=mod(tgrad,M);
fgrad=mod(fgrad,N);  
  
sr=zeros(M,N,W,assert_classname(s,tgrad,fgrad));

fgrad=fgrad+1;
tgrad=tgrad+1;

for w=1:W
    for ii=1:M
        for jj=1:N      
            sr(tgrad(ii,jj),fgrad(ii,jj),w) = sr(tgrad(ii,jj),fgrad(ii,jj),w)+s(ii,jj,w);
        end;
    end;  
end;



