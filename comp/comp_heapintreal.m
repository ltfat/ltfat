function newphase=comp_heapintreal(abss,tgrad,fgrad,a,M, tol,phasetype)
%
%  This function is a wrapper for the comp_heapintreal mex function.

if exist('comp_heapintreal', 'file') == 3
    newphase=comp_heapintreal(abss,tgrad,fgrad,a,M, tol,phasetype);
else
    error('%s: This function needs to be compiled. Please run LTFATMEX.',upper(mfilename));
end