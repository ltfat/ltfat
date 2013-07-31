function c=comp_filterbank(f,g,a,crossover);
%COMP_FILTERBANK  Compute filtering
%
%   If numel(g.h)<=crossover, the routine will use the time-side algorithm for
%   FIR filters, otherwise it will always do the multiplication in the
%   frequency domain.

[L,W]=size(f);
M=numel(g);
c = cell(M,1);

classname=assert_classname(f);

l=(0:L-1).'/L;

% Test if we will need to compute the fft of f
dofft=0;
for m=1:M
   if ~(isfield(g{m},'h') && numel(g{m}.h)<=crossover)
       dofft=1;
   end;
end;

if dofft
    F=fft(f);
end;

if size(a,2)==1
    N=L./a;
else
    N=L./a(:,1).*a(:,2);
    afrac=a(:,1)./a(:,2);
end;

% Divide filters to time domain and frequency domain groups
mFreq = 1:M;
mTime = mFreq(cellfun(@(gEl) isfield(gEl,'h') && numel(gEl.h)<=crossover,g)>0); 
mFreq(mTime) = [];

if ~isempty(mTime)
   % Pick imp. resp. and modulate
   gtime = cellfun(@(gEl)...
           gEl.h.*exp(2*pi*1i*round(gEl.fc*L/2)*(gEl.offset:gEl.offset+numel(gEl.h)-1).'),...
           g(mTime),'UniformOutput',0);
   % Handle realonly
   gtimeId = cellfun(@(gEl) gEl.realonly ,g(mTime));
   gtime(gtimeId>0) = real(gtimeId(gtimeId>0));
   
   % Call the routine
   gskip = cellfun(@(gEl) -gEl.offset ,g(mTime));
   c(mTime) = comp_filterbank_td(f,gtime,a(mTime),gskip,'per');
   
   % Handle realoutput
   realoutputId = mTime(cellfun(@(gEl) isreal(f) && isreal(gEl.h) && gEl.fc==0 ,g(mTime))>0);
   c(realoutputId) = cellfun(@(cEl) real(cEl),c(realoutputId),'UniformOutput',0);

% The original code
%    for mId=1:numel(mTime)
%          m = mTime(mId);
%          realoutput = isreal(f) && isfield(g{m},'h') && isreal(g{m}.h) && g{m}.fc==0;
%          % Use a direct algorithm
%          g_time=circshift(postpad(g{m}.h,L),g{m}.offset).*...
%                 exp(2*pi*1i*round(g{m}.fc*L/2)*l);            
%          if g{m}.realonly
%              g_time=real(g_time);
%          end;
%          g_time=conj(involute(g_time));
%          for n=0:N(m)-1
%              c{m}(n+1,:)=sum(bsxfun(@times,f,circshift(g_time,a(m)*n)));
%          end;
%          if realoutput
%             c{m}=real(c{m});
%          end;
%   end
   
end

for mId=1:numel(mFreq)
    m = mFreq(mId);
      
        % Zero-extend and use a full length fft algorithm. This case can be further
        % optimized by not calling comp_transferfunction, but getting the
        % support correctly.
        
        G=comp_transferfunction(g{m},L);
       
        c{m}=zeros(N(m),W,classname);
        
        if size(a,2)==1
            for w=1:W
                c{m}(:,w)=ifft(sum(reshape(F(:,w).*G,N(m),a(m)),2))/a(m);
            end;
        else
            Llarge=ceil(L/N(m)+1)*N(m);
            amod=Llarge/N(m);
            
            for w=1:W
                c{m}(:,w)=ifft(sum(reshape(fir2long(F(:,w).*G,Llarge),N(m),amod),2))/afrac(m);                
            end;
        end;

end;




