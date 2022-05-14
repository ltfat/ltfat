function f = comp_ifilterbank(c,g,a,L)
%COMP_IFILTERBANK Compute inverse filterbank

%   called by comp_ifilterbank, performs frequency domain
%   filtering for filters with finite impulse responses

M = numel(g);
classname = assert_classname(c{1});


% Divide filters into time domain and frequency domain groups
mFreq = 1:M;
mTime = mFreq(cellfun(@(gEl) isfield(gEl,'h') ,g)>0); 
mFreq(mTime) = [];

f = [];

if ~isempty(mTime)
   % Pick imp. resp.
   gtime = cellfun(@(gEl) gEl.h, g(mTime),'UniformOutput',0);

   % Call the routine
   gskip = cellfun(@(gEl) gEl.offset ,g(mTime));
   f = comp_ifilterbank_td(c(mTime),gtime,a(mTime),L,gskip,'per');
end

if ~isempty(mFreq)
   % Pick frequency domain filters
   gfreq = g(mFreq);
   % Divide filters into the full-length and band-limited groups
   mFreqFullL = 1:numel(gfreq);
   amFreqCell = mat2cell(a(mFreq,:).',size(a,2),ones(1,numel(mFreq)));
   mFreqBL = mFreqFullL(cellfun(@(gEl,aEl) numel(gEl.H)~=L || (numel(aEl)>1 && aEl(2) ~=1), gfreq(:),amFreqCell(:))>0);
   mFreqFullL(mFreqBL) = [];
   
   mFreqFullL = mFreq(mFreqFullL);
   mFreqBL = mFreq(mFreqBL);
   
   F = [];
   if ~isempty(mFreqBL)
      conjG = cellfun(@(gEl) cast(gEl.H,classname), g(mFreqBL),'UniformOutput',0);
      foff = cellfun(@(gEl) gEl.foff, g(mFreqBL));
      % Cast from logical to double.
      realonly = cellfun(@(gEl) cast(isfield(gEl,'realonly') && gEl.realonly,'double'), g(mFreqBL));
      F = comp_ifilterbank_fftbl(c(mFreqBL),conjG,foff,a(mFreqBL,:),realonly);
   end   
   
   if ~isempty(mFreqFullL)
      conjG = cellfun(@(gEl) cast(gEl.H,classname), g(mFreqFullL),'UniformOutput',0);
      
      % In case some of the filters were BL
      if isempty(F)
         F = comp_ifilterbank_fft(c(mFreqFullL),conjG,a(mFreqFullL));
      else
         F = F + comp_ifilterbank_fft(c(mFreqFullL),conjG,a(mFreqFullL));
      end
   end
   
   % In case some of the filters were TD
   if isempty(f)
      f = ifft(F);
   else
      f = f + ifft(F);
   end
end




% W = size(c{1},2);
% M = numel(g);
% classname = assert_classname(c{1});
% 
% f=zeros(L,W,classname);
% 
% % This routine must handle the following cases
% %
% %   * Time-side or frequency-side filters (test for  isfield(g,'H'))
% %
% %   * Cell array or matrix input (test for iscell(c))
% %
% %   * Regular or fractional subsampling (test for info.isfractional)
% 
% 
% for m=1:M
%     conjG=conj(comp_transferfunction(g{m},L));
%         
%     % For Octave 3.6 compatibility
%     conjG=cast(conjG,classname);
%     
%     % Handle fractional subsampling (this implies frequency side filters)
%     if isfield(g{m},'H') && numel(g{m}.H)~=L
%         N=size(c{m},1);
%         Llarge=ceil(L/N)*N;
%         amod=Llarge/N;
%         
%         for w=1:W                        
%             % This repmat cannot be replaced by bsxfun
%             innerstuff=middlepad(circshift(repmat(fft(c{m}(:,w)),amod,1),-g{m}.foff),L);
%             innerstuff(numel(g{m}.H)+1:end) = 0;
%             f(:,w)=f(:,w)+(circshift(innerstuff.*circshift(conjG,-g{m}.foff),g{m}.foff));
%         end;                
%     else
%         if iscell(c)
%             for w=1:W
%                 % This repmat cannot be replaced by bsxfun
%                 f(:,w)=f(:,w)+(repmat(fft(c{m}(:,w)),a(m),1).*conjG);
%             end;
%         else
%             for w=1:W
%                 % This repmat cannot be replaced by bsxfun
%                 f(:,w)=f(:,w)+(repmat(fft(c(:,m,w)),a(m),1).*conjG);
%             end;            
%         end;
%     end;
% end;
% 
% f = ifft(f);
