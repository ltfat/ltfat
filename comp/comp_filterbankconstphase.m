function [newphase,usedmask] = comp_filterbankconstphase(abss,tgrad,fgrad,NEIGH,posInfo,fc,mask,usephase,a,M,N,tol,phasetype)

%sets the mask for the heap integration and acts as a wrapper for
%the corresponding C-function

%% The actual phase construction
sMax = max(abss(:));
absthr = max(sMax)*tol;

if isempty(mask)
    usedmask= zeros(size(abss));
else
    usedmask = cell2mat(mask);
end

%phasetype = 0;
% % Build the phase
if isempty(mask)
     % s is real and positive
     newphase = comp_filterbankheapint(abss,tgrad,fgrad,NEIGH,posInfo,fc,a,M,N,tol(1),phasetype);
 
     % Find all small coefficients and set the mask
     bigenoughidx = abss>absthr(1);
     usedmask(bigenoughidx) = 1;
 else
     newphase=comp_filterbankmaskedheapint(abss,tgrad,fgrad,NEIGH,posInfo,fc,mask,a,M,N,tol(1),...
                                     phasetype,usephase);
     % Find all small coefficients in the unknown phase area
     missingidx = find(usedmask==0);
     bigenoughidx = abss(missingidx)>absthr(1);
     usedmask(missingidx(bigenoughidx)) = 1;
end
% 
% % Do further tol
for ii=2:numel(tol)
     newphase=comp_filterbankmaskedheapint(abss,tgrad,fgrad,NEIGH,posInfo,fc,usedmask,a,M,N,tol(ii),...
                                     phasetype,newphase);
     missingidx = find(usedmask==0);
     bigenoughidx = abss(missingidx)>absthr(ii);
     usedmask(missingidx(bigenoughidx)) = 1;
end

% Convert the mask so it can be used directly for indexing
usedmask = logical(usedmask);
% Assign random values to coefficients below tolerance
zerono = numel(find(~usedmask));
newphase(~usedmask) = rand(zerono,1)*2*pi;