function sr=gabreassign(s,tgrad,fgrad,a,p5);
%GABREASSIGN  Reassign time-frequency distribution
%   Usage:  sr = gabreassign(s,tgrad,fgrad,a);
%           sr = gabreassign(s,tgrad,fgrad,a,'aa');
%
%   `gabreassign(s,tgrad,fgrad,a)` reassigns the values of the positive
%   time-frequency distribution *s* using the phase gradient given by *fgrad*
%   and *tgrad*. The lattice is determined by the time shift *a* and the number
%   of channels deduced from the size of *s*.
%
%   *fgrad* and *tgrad* can be obtained by the routine |gabphasegrad|.
%
%   Examples:
%   ---------
%
%   The following example demonstrates how to manually create a
%   reassigned spectrogram. An easier way is to just call |resgram|:::
%
%     % Create reassigned vector field of the bat signal.
%     a=4;
%     M=100
%     [tgrad, fgrad, c] = gabphasegrad('dgt',bat,'gauss',a,M);
%
%     % Perform the actual reassignment
%     sr = gabreassign(abs(c).^2,tgrad,fgrad,a);
%
%     % Display it using plotdgt
%     plotdgt(sr,a,143000,50);
%  
%   See also: resgram, gabphasegrad
%
%   References: aufl95

% AUTHOR: Peter L. Søndergaard, 2008.
  
error(nargchk(4,5,nargin));


sr=comp_gabreassign(s,tgrad,fgrad,a);




% The following code is currently not actived. It calculates the
% reassigment using anti-aliasing, but it make very little visual
% difference, and it is slower.
  %   [M,N,W]=size(s);
  %   L=N*a;
  %   b=L/M;
    
  %   freqpos=fftindex(M);  
  %   tgrad=bsxfun(@plus,tgrad/b,freqpos);
        
  %   timepos=fftindex(N);
  %   fgrad=bsxfun(@plus,fgrad/a,timepos.');
    
  %   tgrad=round(tgrad);
  %   fgrad=round(fgrad);
    
  %   tgrad=mod(tgrad,M);
  %   fgrad=mod(fgrad,N);  
    
  %   sr=zeros(M,N,W);
    
  %   fk=mod(floor(tgrad),M)+1;
  %   ck=mod(ceil(tgrad),M)+1;
  %   fn=mod(floor(fgrad),N)+1;
  %   cn=mod(ceil(fgrad),N)+1;
    
  %   alpha = fgrad-floor(fgrad);
  %   beta  = tgrad-floor(tgrad);
  %   m1 =(1-alpha).*(1-beta).*s;
  %   m2 =(1-alpha).*beta.*s;
  %   m3 =alpha.*(1-beta).*s;
  %   m4 =alpha.*beta.*s;
  %   for ii=1:M
  %     for jj=1:N
  %       sr(fk(ii,jj),fn(ii,jj))=sr(fk(ii,jj),fn(ii,jj))+m1(ii,jj);
  %       sr(ck(ii,jj),fn(ii,jj))=sr(ck(ii,jj),fn(ii,jj))+m2(ii,jj);
  %       sr(fk(ii,jj),cn(ii,jj))=sr(fk(ii,jj),cn(ii,jj))+m3(ii,jj);
  %       sr(ck(ii,jj),cn(ii,jj))=sr(ck(ii,jj),cn(ii,jj))+m4(ii,jj);
        
  %     end;
  %   end;
  % end;
