function [xo,N]=ref_adaptlasso(ttype,xi,lambda,group);
%TF_ADAPTLASSO   adaptive lasso estimate (hard/soft) in time-frequency domain
%   Usage:  xo=tf_adaptlasso(ttype,x,lambda,group);
%           [xo,N]=tf_grouplasso(ttype,x,lambda,group));
%
%   TF_ADAPTLASSO('hard',x,lambda,'time') will perform
%   time hard adaptive thresholding on x, i.e. coefficients
%   below a column dependent threshold are set to zero. The
%   threshold is computed from the l-1 norm of its
%   time-frequency column.
%
%   TF_ADAPTLASSO('soft',x,lambda,'time') will perform
%   time soft adaptive thresholding on x, i.e. all coefficients
%   below a column dependent threshold are set to zero.
%   The threshold value is substracted from coefficients above
%   the threshold
%
%   TF_ADAPTLASSO(ttype,x,lambda,'frequency') will perform
%   frequency adaptive thresholding on x, i.e. coefficients
%   below a row dependent threshold are set to zero. The
%   threshold is computed from the l-1 norm of its
%   time-frequency row.
%
%   [xo,N]=TF_ADAPTLASSO(ttype,x,lambda,group) additionally returns
%   a number N specifying how many numbers where kept.

complainif_notenoughargs(nargin,4,mfilename);
  
group=lower(group);
ttype=lower(ttype);

tmp = size(xi);
NbFreqBands = tmp(1);
NbTimeSteps = tmp(2);

xo = zeros(size(xi));

if (strcmp(group,'time')),
    for t=1:NbTimeSteps,
        threshold = norm(xi(:,t),1);
        threshold = lambda*threshold /(1+NbFreqBands*lambda);
        mask = abs(xi(:,t)) >= threshold;
        if(strcmp(ttype,'soft'))
            xo(mask,t) = sign(xi(mask,t)) .* (abs(xi(mask,t))-threshold);
        elseif(strcmp(ttype,'hard'))
            xo(mask,t) = sign(xi(mask,t)) .* abs(xi(mask,t));
        end
    end
elseif (strcmp(group,'frequency')),
    for f=1:NbFreqBands,
        threshold = norm(xi(f,:),1);
        threshold = lambda*threshold /(1+NbTimeSteps*lambda);
        mask = abs(xi(f,:)) >= threshold;
        if(strcmp(ttype,'soft'))
            xo(f,mask) = sign(xi(f,mask)) .* (abs(xi(f,mask))-threshold);
        elseif(strcmp(ttype,'hard'))
            xo(f,mask) = sign(xi(f,mask)) .* abs(xi(f,mask));
        end
    end
end

if nargout==2
    signif_map = (abs(xo)>0);
    N = sum(signif_map(:));
end
    

