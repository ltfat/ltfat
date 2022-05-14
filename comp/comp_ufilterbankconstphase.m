function [newphase,usedmask] = comp_ufilterbankconstphase(abss,tgrad,fgrad,fc,mask,usephase,a,tol,phasetype,do_real)

%   part of the heap integration for reconstructing phase from the
%   magnitude of filterbank coefficients

%% DO the heap integration
absthr = max(abss(:))*tol;
if isempty(mask)
    usedmask = zeros(size(abss));
else
    usedmask = mask;
end

%phasetype = 0;
if isempty(mask)
    % Build the phase (calling a MEX file)
    newphase=comp_ufilterbankheapint(abss,tgrad,fgrad,fc,a,do_real,tol(1),phasetype);
    % Set phase of the coefficients below tol to random values
    bigenoughidx = abss>absthr(1);
    usedmask(bigenoughidx) = 1;
else
    newphase=comp_ufilterbankmaskedheapint(abss,tgrad,fgrad,fc,mask,a,do_real,tol(1),phasetype,...
                                usephase);
    % Set phase of small coefficient to random values
    % but just in the missing part
    % Find all small coefficients in the unknown phase area
    missingidx = find(usedmask==0);
    bigenoughidx = abss(missingidx)>absthr(1);
    usedmask(missingidx(bigenoughidx)) = 1;
end

% Do further tol
for ii=2:numel(tol)
    newphase=comp_ufilterbankmaskedheapint(abss,tgrad,fgrad,fc,usedmask,a,do_real,tol(ii),phasetype,...
                                newphase);
    missingidx = find(usedmask==0);
    bigenoughidx = abss(missingidx)>absthr(ii);
    usedmask(missingidx(bigenoughidx)) = 1;                  
end

% Convert the mask so it can be used directly for indexing
usedmask = logical(usedmask);
% Assign random values to coefficients below tolerance
zerono = numel(find(~usedmask));
newphase(~usedmask) = rand(zerono,1)*2*pi;