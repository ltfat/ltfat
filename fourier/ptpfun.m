function g = ptpfun(L,w,varargin)
%PTPFUN Sampled, periodized totally positive function of finite type
%   Usage: g=ptpfun(L,w)
%          g=ptpfun(L,w,width)
%
%   Input parameters:
%         L    : Window length.
%         w    : Vector of reciprocals $w_j=1/\delta_j$ in Fourier representation of *g*
%         width: Integer stretching factor for the essential support of g 
%
%   Output parameters:
%         g    : The periodized totally positive function.
%
%   `ptpfun(L,w)` computes samples of a periodized totally positive
%   function of finite type >=2 with weights *w* for a system of length *L*.
%   The Fourier representation of the continuous TP function is given as:
%
%   ..             m
%      ghat(f) = prod (1+2pijf/w(i))^(-1),
%                 i=1
%
%   .. math:: \hat{g}(\xi)=\prod_{i=1}^{m}\left(1+2\pi i j\xi /w(i)\right)^{-1},
%
%   where $m$=`numel(w)`$\geq 2$. The samples are obtained by discretizing
%   the Zak transform of the function.
%
%   *w* controls the function decay in the time domain. More specifically
%   the function decays as $exp(max(w)x)$ for $x->\infty$ and $exp(min(w)x)$
%   for $x->-\infty$ assuming *w* contains both positive and negative
%   numbers.
%
%   `ptpfun(L,w,width)` additionally stretches the function by a factor of 
%   *width*.
%
%   See also: dgt, ptpfundual, gabdualnorm, setnorm
%
%   References: grst13 kl12 bagrst14 klst14
%

%   AUTHORS: Joachim Stoeckler, Tobias Kloos  2012, 2014

complainif_notenoughargs(nargin,2,upper(mfilename));
complainif_notposint(L,'L',upper(mfilename));

if isempty(w) || ~isnumeric(w) || numel(w)<2
    error(['%s: w must be a nonempty numeric vector with at least',...
           ' 2 elements.'], upper(mfilename));
end

if any(w==0)
    error('%s: All weights w must be nonzero.', upper(mfilename));
end

% Define initial value for flags and key/value pairs.
%definput.import={'setnorm'};
definput.keyvals.width=floor(sqrt(L));
[flags,~,width]=ltfatarghelper({'width'},definput,varargin);
complainif_notposint(width,'width',upper(mfilename));

w = -sort(w(:))*L/width;
x = [0:L-1]'/L;
n = length(w);
x = repmat(x,1,2*n-1) + repmat([-n+1:n-1],L,1);
x = x(:)';
Y = zeros(n-1,length(x));

for k = 1:n-1
    if w(k) == w(k+1)
        if w(k) < 0
            Y(k,(n-1)*L+1:(n+1)*L) = w(k)^2/(1-exp(w(k)))^2*[x((n-1)*L+1:n*L).*exp(w(k)*x((n-1)*L+1:n*L)) , ...
                (2-x(n*L+1:(n+1)*L)).*exp(w(k)*x(n*L+1:(n+1)*L))];
        else
            Y(k,(n-1)*L+1:(n+1)*L) = w(k)^2/(exp(-w(k))-1)^2*[x((n-1)*L+1:n*L).*exp(w(k)*(x((n-1)*L+1:n*L)-2)) , ...
                (2-x(n*L+1:(n+1)*L)).*exp(w(k)*(x(n*L+1:(n+1)*L)-2))];
        end
    else
        if w(k) < 0 && w(k+1) < 0
            Y(k,(n-1)*L+1:(n+1)*L) = w(k)/(1-exp(w(k)))*w(k+1)/(1-exp(w(k+1)))*[(exp(w(k)*x((n-1)*L+1:n*L))-exp(w(k+1)*x((n-1)*L+1:n*L)))/(w(k)-w(k+1)) , ...
                (exp(w(k))*exp(w(k+1)*(x(n*L+1:(n+1)*L)-1))-exp(w(k+1))*exp(w(k)*(x(n*L+1:(n+1)*L)-1)))/(w(k)-w(k+1))];
        elseif w(k) > 0 && w(k+1) < 0
            Y(k,(n-1)*L+1:(n+1)*L) = w(k)/(exp(-w(k))-1)*w(k+1)/(1-exp(w(k+1)))*[(exp(w(k)*(x((n-1)*L+1:n*L)-1))-exp(-w(k))*exp(w(k+1)*x((n-1)*L+1:n*L)))/(w(k)-w(k+1)) , ...
                (exp(w(k+1)*(x(n*L+1:(n+1)*L)-1))-exp(w(k+1))*exp(w(k)*(x(n*L+1:(n+1)*L)-2)))/(w(k)-w(k+1))];
        elseif w(k) < 0 && w(k+1) > 0
            Y(k,(n-1)*L+1:(n+1)*L) = w(k+1)/(exp(-w(k+1))-1)*w(k)/(1-exp(w(k)))*[(exp(w(k+1)*(x((n-1)*L+1:n*L)-1))-exp(-w(k+1))*exp(w(k)*x((n-1)*L+1:n*L)))/(w(k+1)-w(k)) , ...
                (exp(w(k)*(x(n*L+1:(n+1)*L)-1))-exp(w(k))*exp(w(k+1)*(x(n*L+1:(n+1)*L)-2)))/(w(k+1)-w(k))];
        elseif w(k) > 0 && w(k+1) > 0
            Y(k,(n-1)*L+1:(n+1)*L) = w(k)/(exp(-w(k))-1)*w(k+1)/(exp(-w(k+1))-1)*[(exp(-w(k+1))*exp(w(k)*(x((n-1)*L+1:n*L)-1))-exp(-w(k))*exp(w(k+1)*(x((n-1)*L+1:n*L)-1)))/(w(k)-w(k+1)) , ...
                (exp(w(k+1)*(x(n*L+1:(n+1)*L)-2))-exp(w(k)*(x(n*L+1:(n+1)*L)-2)))/(w(k)-w(k+1))];
        end
    end
end

for k = 2:n-1
    for j = 1:n-k
        if w(j) == w(j+k)
            if w(j) < 0
                Y(j,(k-1)*L+1:end) = -w(j)/(1-exp(w(j)))*(x((k-1)*L+1:end)/k .* Y(j,(k-1)*L+1:end) + ...
                    exp(w(j))*(k+1-x((k-1)*L+1:end))/k .* Y(j,(k-2)*L+1:end-L));
            else
                Y(j,(k-1)*L+1:end) = -w(j)/(exp(-w(j))-1)*(exp(-w(j))*x((k-1)*L+1:end)/k .* Y(j,(k-1)*L+1:end) + ...
                    (k+1-x((k-1)*L+1:end))/k .* Y(j,(k-2)*L+1:end-L));
            end
        else
            if w(j) < 0 && w(j+k) < 0
                Y(j,(k-1)*L+1:end) = ( (-w(j)*(Y(j+1,(k-1)*L+1:end) - exp(w(j))*Y(j+1,(k-2)*L+1:end-L)))/(1-exp(w(j))) - ...
                    (-w(j+k)*(Y(j,(k-1)*L+1:end) - exp(w(j+k))*Y(j,(k-2)*L+1:end-L)))/(1-exp(w(j+k))) )/ ...
                    (w(j+k)-w(j));
            elseif w(j) > 0 && w(j+k) < 0
                Y(j,(k-1)*L+1:end) = ( (-w(j)*(exp(-w(j))*Y(j+1,(k-1)*L+1:end) - Y(j+1,(k-2)*L+1:end-L)))/(exp(-w(j))-1) - ...
                    (-w(j+k)*(Y(j,(k-1)*L+1:end) - exp(w(j+k))*Y(j,(k-2)*L+1:end-L)))/(1-exp(w(j+k))) )/ ...
                    (w(j+k)-w(j));
            elseif w(j) < 0 && w(j+k) > 0
                Y(j,(k-1)*L+1:end) = ( (-w(j)*(Y(j+1,(k-1)*L+1:end) - exp(w(j))*Y(j+1,(k-2)*L+1:end-L)))/(1-exp(w(j))) - ...
                    (-w(j+k)*(exp(-w(j+k))*Y(j,(k-1)*L+1:end) - Y(j,(k-2)*L+1:end-L)))/(exp(-w(j+k))-1) )/ ...
                    (w(j+k)-w(j));
            elseif w(j) > 0 && w(j+k) > 0
                Y(j,(k-1)*L+1:end) = ( (-w(j)*(exp(-w(j))*Y(j+1,(k-1)*L+1:end) - Y(j+1,(k-2)*L+1:end-L)))/(exp(-w(j))-1) - ...
                    (-w(j+k)*(exp(-w(j+k))*Y(j,(k-1)*L+1:end) - Y(j,(k-2)*L+1:end-L)))/(exp(-w(j+k))-1) )/ ...
                    (w(j+k)-w(j));
            end
        end
    end
end

if n == 1
    if w < 0
        g = -w/(1-exp(w))*exp(w*x(1:L)) * sqrt(width) / L;
    else
        g = -w/(exp(-w)-1)*exp(w*(x(1:L)-1)) * sqrt(width) / L;
    end
else
    g = sum(reshape(Y(1,end-n*L+1:end),L,n),2) * sqrt(width) / L;
end

%g = setnorm(g(:),flags.norm);
g = g(:);

end
