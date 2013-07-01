function [g,a,fc]=erbfilters(fs,varargin)

definput.keyvals.L=[];
definput.keyvals.N=[];
definput.flags.uniform  = {'nonuniform','uniform'};
definput.flags.real     = {'real','complex'};
definput.flags.sampling = {'regsampling','fractional'};

[flags,kv,L]=ltfatarghelper({'L'},definput,varargin);

% Construct the Erb filterbank

N=kv.N;
if isempty(N)
    N=ceil(freqtoerb(fs/2))+1;
end;

fc=erbspace(0,fs/2,N);
% "*3" is just a heuristic, no justification
fsupp=round(audfiltbw(fc)*4);

% Improve the scaling of the first and last channel
scal=ones(1,N);
scal(1)=scal(1)/sqrt(2);
scal(N)=scal(N)/sqrt(2);


if flags.do_nonuniform
    % Do the non-uniform case
    % Energy scaling works best
    g=blfilter('hanning',fsupp,fc,'fs',fs,'scal',scal,'2');
    
    % Find suitable channel subsampling rates
    aprecise=round(fs./fsupp/2); % "/2" is just a heuristic, no justification
    aprecise=aprecise(:);
    
    if flags.do_fractional
        Nfilts=round(L./aprecise);
        a=[repmat(L,N,1),Nfilts];                
        
    else
        a=ceil23(aprecise); % Grow "a" to the next composite number
    
        % Determine the minimal transform length
        L=filterbanklength(1,a);
    end;

else
    % Do the uniform case
    % Peak-frequency scaling works best
    g=blfilter('hanning',fsupp,fc,'fs',fs,'scal',scal,'inf');
    a=3;

end;


