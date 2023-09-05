function sr=comp_gabreassignreal(s,tgrad,fgrad,a,M);
%COMP_GABREASSIGNREAL  Reassign time-frequency distribution.
%   Usage:  sr = comp_gabreassign(s,tgrad,fgrad,a);
%
%   `comp_gabreassignreal(s,tgrad,fgrad,a,M)` will reassign the values of the 
%   positive time-frequency distribution *s* using the instantaneous time and 
%   frequency *fgrad* and *tgrad*. The lattice is determined by the time shift 
%   *a* and the number of channels *M*.
%
%   See also: gabreassignreal
%
%   References: aufl95

%   AUTHOR : Peter L. Søndergaard.
%   TESTING: OK
%   REFERENCE: OK

[M2,N,W]=size(s);
L=N*a;
b=L/M;

freqpos=(0:M2-1).';
tgrad=bsxfun(@plus,tgrad/b,freqpos);

timepos=fftindex(N);
fgrad=bsxfun(@plus,fgrad/a,timepos.');

tgrad=round(tgrad);
fgrad=round(fgrad);

tgrad=mod(tgrad,M);
fgrad=mod(fgrad,N);  
  
sr=zeros(M2,N,W,assert_classname(s,tgrad,fgrad));

fgrad=fgrad+1;
tgrad=tgrad+1;

% In theory, this should not be necessary, but we make sure that fgrad is
% in the range 1 to M2
zeroIdxs = (tgrad > floor(3/4*M)+1);
tgrad(zeroIdxs) = 1;
tgrad(tgrad > M2 & ~ zeroIdxs) = M2;

for w=1:W
    for ii=1:M2
        for jj=1:N      
            sr(tgrad(ii,jj),fgrad(ii,jj),w) = sr(tgrad(ii,jj),fgrad(ii,jj),w)+s(ii,jj,w);
        end;
    end;  
end;



