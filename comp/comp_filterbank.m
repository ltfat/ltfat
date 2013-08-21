function c=comp_filterbank(f,g,a);
%COMP_FILTERBANK  Compute filtering
%
%   Function groups filters in g according to a presence of .h and .H
%   fields. If .H is present, it is further decided whether it is a full
%   frequency response or a band-limited freq. resp.
%

[L,W]=size(f);
M=numel(g);
c = cell(M,1);

% Divide filters into time domain and frequency domain groups
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
   % Filtering in the frequency domain
   F=fft(f);
   % Pick frequency domain filters
   gfreq = g(mFreq);
   % Divide filters into the full-length and band-limited groups
   mFreqFullL = 1:numel(gfreq);
   amFreqCell = mat2cell(a(mFreq,:).',size(a,2),ones(1,numel(mFreq)));
   mFreqBL = mFreqFullL(cellfun(@(gEl,aEl) numel(gEl.H)~=L || (numel(aEl)>1 && aEl(2) ~=1), gfreq(:),amFreqCell(:))>0);
   mFreqFullL(mFreqBL) = [];
   
   mFreqFullL = mFreq(mFreqFullL);
   mFreqBL = mFreq(mFreqBL);
   
   if ~isempty(mFreqFullL)
      G = cellfun(@(gEl) gEl.H, g(mFreqFullL),'UniformOutput',0);
      c(mFreqFullL) = comp_filterbank_fft(F,G,a(mFreqFullL));
   end
   
   if ~isempty(mFreqBL)
      G = cellfun(@(gEl) gEl.H, g(mFreqBL),'UniformOutput',0);
      foff = cellfun(@(gEl) gEl.foff, g(mFreqBL));
      % Cast from logical to double.
      realonly = cellfun(@(gEl) cast(isfield(gEl,'realonly') && gEl.realonly,'double'), g(mFreqBL));
      c(mFreqBL) = comp_filterbank_fftbl(F,G,foff,a(mFreqBL,:),realonly);
   end
end

% for mId=1:numel(mFreq)
%     m = mFreq(mId);
%     G = g{m}.H;
%     Lg = numel(G);
%     c{m}=zeros(N(m),W,classname);
%         % Band-limited case. 
%         if (g{m}.foff~=0 && Lg<L) && size(a,2)==1 && ~g{m}.realonly && a(m)>1
% 
%            fsuppStart = g{m}.foff;
%            fsuppEnd = g{m}.foff+numel(g{m}.H);
%            
%            % Partition the support to intervals of length N(m)
%            fsuppSidx = floor(fsuppStart/N(m));
%            fsuppEidx = ceil(fsuppEnd/N(m));
%            fsuppS = fsuppSidx*N(m);
%            fsuppE = fsuppEidx*N(m);
%            intNo = (fsuppE-fsuppS)/N(m);
%            Gtmp = zeros(intNo*N(m),1);
%            %Gtmp(fsuppStart-fsuppS+1:end-(fsuppE-fsuppEnd)) = G;
%            
%            fsuppRangeSmall = mod(fsuppStart:fsuppEnd-1,L)+1;
%            %fsuppRange = mod(fsuppS:fsuppS+intNo*N(m)-1,L)+1;
%            for w=1:W
%               interstuff = G.*F(fsuppRangeSmall,w);
%               Gtmp(fsuppStart-fsuppS+1:end-(fsuppE-fsuppEnd)) = interstuff;
%               c{m}(:,w) = ifft(sum(reshape(Gtmp,N(m),intNo),2))/a(m);
%            end
%         else
%            G = comp_transferfunction(g{m},L);
%            if size(a,2)==1
% 
%              for w=1:W
%                 c{m}(:,w)=ifft(sum(reshape(F(:,w).*G,N(m),a(m)),2))/a(m);
%              end;
%            else
%               % Fractional case
%               %Llarge=ceil(L/N(m)+1)*N(m);
%               
%               postpadL = ceil(max([N(m),numel(g{m}.H)])/N(m))*N(m);
%               foff = g{m}.foff + floor(numel(g{m}.H)/2);
%               amod=postpadL/N(m);
%               
%               for w=1:W
%                  c{m}(:,w)=ifft(sum(reshape(circshift(middlepad(circshift(F(:,w).*G,-foff),postpadL),foff),N(m),amod),2))/afrac(m);                
%               end;
%            end;
%         end;
% end
%end




