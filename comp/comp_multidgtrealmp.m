function [c,atoms,iter,status] = comp_multidgtrealmp(fpad,g,a,M,phasetype,...
                            kernthr,errdb,iter1,iter2, pedanticsearch, algorithm )
                        
                        
if exist('comp_multidgtrealmp', 'file') == 3
    [c,atoms,iter,status] = comp_multidgtrealmp(fpad,g,a,M,phasetype,...
                            kernthr,errdb,iter1,iter2, pedanticsearch, algorithm );
else
    error('%s: This function needs to be compiled. Please run LTFATMEX.',upper(mfilename));
end