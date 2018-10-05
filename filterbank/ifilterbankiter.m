function [f,relres,iter]=ifilterbankiter(c,g,a,varargin)
%IFILTERBANKITER  Filter bank iterative inversion
%   Usage:  f=ifilterbankiter(c,g,a);
%
%   `ifilterbankiter(c,g,a)` iteratively synthesizes a signal *f* from the
%   coefficients *c* which were obtained using the filters stored in *g* for
%   a channel subsampling rate of *a* (the hop-size).
%
%   The filter bank *g* and the subsampling rate *a* must be the same
%   as used in |filterbank| or |ufilterbank|.
%
%   This function is useful if there is no way how to explicitly compute
%   a dual system using |filterbankdual| or |filterbankrealdual|.
%
%   Additional parameters
%   ---------------------
%
%   The function calls |frsyniter| and passes all the optional arguments to it.
%   Please refer to help of |frsyniter| for further details.
%
%   Please note that by default, the function expects filterbank *g*
%   to be created for real signals i.e. *g* cover only the positive frequencies.
%   Additional flag 'complex' is required if the filterbank is defined for 
%   positive and negative frequencies.
%
%   Examples:
%   ---------
%
%   The following example compares convergence rates of CG and PCG for a
%   filterbank which forms a frame, but it is neither uniform or painless:::
%
%       [f,fs] = greasy; L = size(f,1);
%       [g,a,fc]=erbfilters(fs,L,'fractional','bwmul',0.6,'redmul',4/5,'complex');
%       filterbankfreqz(g,a,L,'plot','linabs');
%       % Filterbankdual does not work
%       try
%           gd=filterbankdual(g,a,L);
%       catch
%           disp('FILTERBANKDUAL exited with error.');
%       end
%
%       c = filterbank(f,g,a);
%       [fpcg,~,iterpcg] = ifilterbankiter(c,g,a,'complex','pcg');
%       [fcg,~,itercg] = ifilterbankiter(c,g,a,'complex','cg');
%
%       fprintf('CG achieved error %e in %d iterations.\n',norm(f-fcg), itercg);
%       fprintf('PCG achieved error %e in %d iterations.\n',norm(f-fpcg), iterpcg);
%
%   Similar example with real filterbank:::
%
%       [f,fs] = greasy; L = size(f,1);
%       [g,a,fc]=erbfilters(fs,L,'fractional','bwmul',0.6,'redmul',4/5);
%       filterbankfreqz(g,a,L,'plot','linabs');
%       % Filterbankrealdual does not work
%       try
%           gd=filterbankrealdual(g,a,L);
%       catch
%           disp('FILTERBANKREALDUAL exited with error.');
%       end
%
%       c = filterbank(f,g,a);
%       [fpcg,~,iterpcg] = ifilterbankiter(c,g,a,'pcg');
%       [fcg,~,itercg] = ifilterbankiter(c,g,a,'cg');
%
%       fprintf('CG achieved error %e in %d iterations.\n',norm(f-fcg), itercg);
%       fprintf('PCG achieved error %e in %d iterations.\n',norm(f-fpcg), iterpcg);
%
%   See also: filterbank, ufilterbank, ifilterbank, filterbankdual
%
%   References: ltfatnote027 nehobaprpide18

complainif_notenoughargs(nargin,3,'IFILTERBANKITER');

tolchooser.double=1e-9;
tolchooser.single=1e-5;

definput.keyvals.Ls=[];
definput.keyvals.tol=tolchooser.(classofc(c));
definput.keyvals.maxit=100;
definput.flags.alg={'pcg','cg'};
definput.flags.real={'real','complex'};
[flags,kv,Ls]=ltfatarghelper({'Ls'},definput,varargin);


if isnumeric(c)
    if flags.do_real
        F = frame('ufilterbankreal',g,a,numel(g));
    else
        F = frame('ufilterbank',g,a,numel(g));
    end
else
    if flags.do_real
        F = frame('filterbankreal',g,a,numel(g));
    else
        F = frame('filterbank',g,a,numel(g));
    end
end

[f,relres,iter] = frsyniter(F,framenative2coef(F,c),'Ls',Ls,flags.alg,'maxit',kv.maxit,'tol',kv.tol);

function cl = classofc(c)
if iscell(c)
    cl = class(c{1});
else
    cl = class(c);
end

