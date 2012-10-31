function fout = extendBoundary(f,extLen,ext)
% 

fout = zeros(length(f) + 2*extLen,1);
fout(extLen+1:end-extLen) = f;

legalExtLen = min([length(f),extLen]);

% zero padding by default
% ext: 'per','zpd','sym','symw','asym','asymw','ppd','sp0'
if(strcmp(ext,'perdec')) % possible last sample replication
    if(mod(length(f),2)==1)
        fout(end-extLen+1) = f(end);
        fout(end-extLen+2:end-extLen+1+legalExtLen-1) = f(1:legalExtLen-1);
        fout(extLen) = f(end);
        fout(1:legalExtLen-1)= f(end-legalExtLen+2:end);
    else
        fout(1:legalExtLen) = f(end-legalExtLen+1:end);
        fout(end-extLen+1:end-extLen+legalExtLen) = f(1:legalExtLen);
    end
elseif(strcmp(ext,'per') || strcmp(ext,'ppd'))
    fout(1:legalExtLen) = f(end-legalExtLen+1:end);
    fout(end-extLen+1:end-extLen+legalExtLen) = f(1:legalExtLen);
elseif(strcmp(ext,'sym'))
    fout(1:legalExtLen) = f(legalExtLen:-1:1);
    fout(end-extLen+1:end-extLen+legalExtLen) = f(end:-1:end-legalExtLen+1);
elseif(strcmp(ext,'symw'))
    legalExtLen = min([length(f)-1,extLen]);
    fout(1:legalExtLen) = f(legalExtLen+1:-1:2);
    fout(end-extLen+1:end-extLen+legalExtLen) = f(end-1:-1:end-legalExtLen);
elseif(strcmp(ext,'asym'))
    fout(1:legalExtLen) = -f(legalExtLen:-1:1);
    fout(end-extLen+1:end-extLen+legalExtLen) = -f(end:-1:end-legalExtLen+1);
elseif(strcmp(ext,'asymw'))
    legalExtLen = min([length(f)-1,extLen]);
    fout(1:legalExtLen) = -f(legalExtLen+1:-1:2);
    fout(end-extLen+1:end-extLen+legalExtLen) = -f(end-1:-1:end-legalExtLen);
elseif(strcmp(ext,'sp0'))
    fout(1:extLen) = f(1);
    fout(end-extLen+1:end) = f(end);
end