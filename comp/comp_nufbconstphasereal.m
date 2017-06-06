function [newphase, usedmask] = comp_nufbconstphasereal(s,tgrad,fgrad,N,NEIGH,posInfo,cfreq,a,M,tol,mask,usephase)

%N = cellfun(@(sEl) size(sEl,1),s);
chanStart = [0;cumsum(N)];

%s = cell2mat(s);
%tgrad = cell2mat(tgrad);
%fgrad = cell2mat(fgrad);

W = size(s,2);

sMax = max(s(:));
absthr = max(sMax)*tol;

if isempty(mask)
    usedmask= zeros(chanStart(end),W);
else
    usedmask = cell2mat(mask);
end

phasetype = 0;
% % Build the phase
% if isempty(mask)
%     % s is real and positive
%     newphase=comp_heapintreal_ufb(s,tgrad,fgrad,cfreq,a,M,tol(1),1);
% 
%     % Find all small coefficients and set the mask
%     bigenoughidx = s>absthr(1);
%     usedmask(bigenoughidx) = 1;
% else
%     newphase=comp_maskedheapintreal_ufb(s,tgrad,fgrad,cfreq,mask,a,M,tol(1),...
%                                     1,usephase);
%     % Find all small coefficients in the unknown phase area
%     missingidx = find(usedmask==0);
%     bigenoughidx = s(missingidx)>absthr(1);
%     usedmask(missingidx(bigenoughidx)) = 1;
% end
% 
% % Do further tol
% for ii=2:numel(tol)
%     newphase=comp_maskedheapintreal_ufb(s,tgrad,fgrad,cfreq,usedmask,a,M,tol(ii),...
%                                     1,newphase);
%     missingidx = find(usedmask==0);
%     bigenoughidx = s(missingidx)>absthr(ii);
%     usedmask(missingidx(bigenoughidx)) = 1;
% end
newphase = comp_heapintreal_fb(s,tgrad,fgrad,NEIGH,posInfo,cfreq,a,M,N,chanStart,tol,phasetype);


% % Convert the mask so it can be used directly for indexing
% usedmask = logical(usedmask);
% % Assign random values to coefficients below tolerance
% zerono = numel(find(~usedmask));
% newphase(~usedmask) = rand(zerono,1)*2*pi;

% Convert to cell array again
%usedmask = mat2cell(usedmask,N,W);
%newphase = mat2cell(newphase,N,W);