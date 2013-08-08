function c=comp_filterbank(f,g,a);
%COMP_FILTERBANK  Compute filtering
%
%   If numel(g.h)<=crossover, the routine will use the time-side algorithm for
%   FIR filters, otherwise it will always do the multiplication in the
%   frequency domain.

[L,W]=size(f);
M=numel(g);
c = cell(M,1);

classname=assert_classname(f);

if size(a,2)==1
    N=L./a;
else
    N=L./a(:,1).*a(:,2);
    afrac=a(:,1)./a(:,2);
end;

% Divide filters to time domain and frequency domain groups
mFreq = 1:M;
mTime = mFreq(cellfun(@(gEl) isfield(gEl,'h') ,g)>0); 
mFreq(mTime) = [];

if ~isempty(mTime)
   % Pick imp. resp.
   gtime = cellfun(@(gEl) gEl.h, g(mTime),'UniformOutput',0);

   % Call the routine
   gskip = cellfun(@(gEl) gEl.offset ,g(mTime));
   c(mTime) = comp_filterbank_td(f,gtime,a(mTime),gskip,'per');
end

if ~isempty(mFreq)
   F=fft(f);
end

for mId=1:numel(mFreq)
    m = mFreq(mId);
        G = g{m}.H;
        Lg = numel(G);
        c{m}=zeros(N(m),W,classname);
        % Band-limited case. 
        if (g{m}.foff~=0 && Lg<L) && size(a,2)==1 && ~g{m}.realonly && a(m)>1

           fsuppStart = g{m}.foff;
           fsuppEnd = g{m}.foff+numel(g{m}.H);
           
           % Partition the support to intervals of length N(m)
           fsuppSidx = floor(fsuppStart/N(m));
           fsuppEidx = ceil(fsuppEnd/N(m));
           fsuppS = fsuppSidx*N(m);
           fsuppE = fsuppEidx*N(m);
           intNo = (fsuppE-fsuppS)/N(m);
           Gtmp = zeros(intNo*N(m),1);
           %Gtmp(fsuppStart-fsuppS+1:end-(fsuppE-fsuppEnd)) = G;
           
           fsuppRangeSmall = mod(fsuppStart:fsuppEnd-1,L)+1;
           %fsuppRange = mod(fsuppS:fsuppS+intNo*N(m)-1,L)+1;
           for w=1:W
              interstuff = G.*F(fsuppRangeSmall,w);
              Gtmp(fsuppStart-fsuppS+1:end-(fsuppE-fsuppEnd)) = interstuff;
              c{m}(:,w) = ifft(sum(reshape(Gtmp,N(m),intNo),2))/a(m);
           end
        else
           G = comp_transferfunction(g{m},L);
           if size(a,2)==1
             for w=1:W
                c{m}(:,w)=ifft(sum(reshape(F(:,w).*G,N(m),a(m)),2))/a(m);
             end;
           else
              % Fractional case
              Llarge=ceil(L/N(m)+1)*N(m);
              amod=Llarge/N(m);
            
              for w=1:W
                 c{m}(:,w)=ifft(sum(reshape(fir2long(F(:,w).*G,Llarge),N(m),amod),2))/afrac(m);                
              end;
           end;
        end;

end;




