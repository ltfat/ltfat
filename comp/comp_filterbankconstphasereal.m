function [newphase, usedmask] = comp_filterbankconstphasereal(s,tgrad,fgrad,cfreq,a,M,tol,mask,usephase)

absthr = max(s(:))*tol;

if isempty(mask)
    usedmask = zeros(size(s));
else
    usedmask = mask;
end

% Build the phase
if isempty(mask)
    % s is real and positive
    newphase=comp_heapintreal_ufb(s,tgrad,fgrad,cfreq,a,M,tol(1),1);

    % Find all small coefficients and set the mask
    bigenoughidx = s>absthr(1);
    usedmask(bigenoughidx) = 1;
else
    newphase=comp_maskedheapintreal_ufb(s,tgrad,fgrad,cfreq,mask,a,M,tol(1),...
                                    1,usephase);
    % Find all small coefficients in the unknown phase area
    missingidx = find(usedmask==0);
    bigenoughidx = s(missingidx)>absthr(1);
    usedmask(missingidx(bigenoughidx)) = 1;
end

% Do further tol
for ii=2:numel(tol)
    newphase=comp_maskedheapintreal_ufb(s,tgrad,fgrad,cfreq,usedmask,a,M,tol(ii),...
                                    1,newphase);
    missingidx = find(usedmask==0);
    bigenoughidx = s(missingidx)>absthr(ii);
    usedmask(missingidx(bigenoughidx)) = 1;
end


% Convert the mask so it can be used directly for indexing
usedmask = logical(usedmask);
% Assign random values to coefficients below tolerance
zerono = numel(find(~usedmask));
newphase(~usedmask) = rand(zerono,1)*2*pi;