function h=ref_pfilt(f,g,a)
%REF_PFILT  Reference pfilt handling structs
%   Usage:  h=pfilt(f,g);
%           h=pfilt(f,g,a,dim);
%
%   `pfilt(f,g)` applies the filter *g* to the input *f*. If *f* is a
%   matrix, the filter is applied along each column.
%
%   `pfilt(f,g,a)` does the same, but downsamples the output keeping only
%   every a'th sample (starting with the first one).
%
%   `pfilt(f,g,a,dim)` filters along dimension dim. The default value of
%   [] means to filter along the first non-singleton dimension.
%
%   The filter *g* can be a vector, in which case the vector is treated
%   as a zero-delay FIR filter.
%
%   The filter *g* can be a cell array. The following options are
%   possible:
%
%     * If the first element of the cell array is the name of one of the
%       windows from |firwin|, the whole cell array is passed onto
%       |firfilter|.
%
%     * If the first element of the cell array is `'bl'`, the rest of the
%     cell array is passed onto |blfilter|.
%
%     * If the first element of the cell array is `'pgauss'`, `'psech'`,
%       the rest of the parameters is passed onto the respective
%       function. Note that you do not need to specify the length *L*.
%
%   The coefficients obtained from filtering a signal *f* by a filter *g* are
%   defined by
%
%   ..          L-1
%      c(n+1) = sum f(l+1) * g(an-l+1)
%               l=0
%
%   .. math:: c\left(n+1\right)=\sum_{l=0}^{L-1}f\left(l+1\right)g\left(an-l+1\right)
%
%   where $an-l$ is computed modulo $L$.
%
%   See also: pconv

if nargin<3
    a=1;
end;

[L,W]=size(f);

l=(0:L-1).'/L;
if isstruct(g)
    if isfield(g,'h')
        g_time=circshift(postpad(g.h,L),g.offset).*exp(2*pi*1i*(round(g.fc*L/2))*l);
        G=fft(g_time);
    elseif isfield(g,'H')
        G=circshift(postpad(g.H(L),L),g.foff(L)).*exp(-2*pi*1i*round(g.delay)*l);  
    else
       error('%s: Unknown filter definition.',upper(mfilename));
    end;
    
    if g.realonly
        G=(G+involute(G))/2;
    end;

else
    G=fft(fir2long(g,L));
end;

N=L/a;
h=zeros(N,W,assert_classname(f));
for w=1:W
    F=fft(f(:,w));
    h(:,w)=ifft(sum(reshape(F.*G,N,a),2))/a;
end;

