function c=comp_filterbank_fftbl(F,G,foff,a,realonly)
%COMP_FILTERBANK_FFTBL  Compute filtering in FD
%
%   does the same as comp_filterbank_fft, but for filters
%   with bandlimited frequency responses

M = numel(G);
[L,W] = size(F);
c = cell(M,1);
if size(a,2)>1
   afrac=a(:,1)./a(:,2);
else
   afrac = a(:,1);
end

N = L./afrac;
assert(all(N-round(N)<1e-6),'%s: Bad output length. \n',upper(mfilename));
N = round(N);

fsuppRangeSmall = cellfun(@(fEl,GEl) mod([fEl:fEl+numel(GEl)-1].',L)+1 ,num2cell(foff),G,'UniformOutput',0);


for m=1:M
    c{m}=zeros(N(m),W,assert_classname(F,G{m}));

    for w=1:W
        Ftmp = F(fsuppRangeSmall{m},w).*G{m};
        postpadL = ceil(max([N(m),numel(G{m})])/N(m))*N(m);
        Ftmp = postpad(Ftmp,postpadL);
        
        Ftmp = sum(reshape(Ftmp,N(m),numel(Ftmp)/N(m)),2);
        
        Ftmp = circshift(Ftmp,foff(m));
        
        c{m}(:,w)=ifft(Ftmp)/afrac(m);
    end;
end


% Handle the real only as a separate filter using recursion
realonlyRange = 1:M;
realonlyRange = realonlyRange(realonly>0);

if ~isempty(realonlyRange)
   Gconj = cellfun(@(gEl) conj(gEl(end:-1:1)),G(realonlyRange),'UniformOutput',0);
   LG = cellfun(@(gEl) numel(gEl),Gconj);
   foffconj = -L+mod(L-foff(realonlyRange)-LG,L)+1;
   aconj = a(realonlyRange,:);

   cconj = comp_filterbank_fftbl(F,Gconj,foffconj,aconj,0);
   for ii=1:numel(cconj)
      c{realonlyRange(ii)} = (c{realonlyRange(ii)} + cconj{ii})/2;
   end
end
