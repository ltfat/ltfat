function newphase=comp_maskedheapint(abss,tgrad,fgrad,a,tol,phasetype,usephase)
%
%  This function is a wrapper for the comp_heapint mex function.

if exist('comp_maskedheapint', 'file') == 3
    newphase=comp_maskedheapint(abss,tgrad,fgrad,a,tol,phasetype,usephase);
else
    error('%s: This function needs to be compiled. Please run LTFATMEX.',upper(mfilename));
end