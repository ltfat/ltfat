function newphase=comp_heapint(abss,tgrad,fgrad,a,tol,phasetype)
%
%  This function is a wrapper for the comp_heapint mex function.

if exist('comp_heapint', 'file') == 3
    newphase=comp_heapint(abss,tgrad,fgrad,a,tol,phasetype);
else
    error('%s: This function needs to be compiled. Please run LTFATMEX.',upper(mfilename));
end