function delays = lowdiscrepancy(name, varargin)
%LOWDISCREPANCY  Returns a low discrepancy sequence
%   Usage: delays=lowdiscrepancy(name)
%
%   Input parameters:
%         name  : Name of the low discrepancy sequence
%   Output parameters:
%         delays     : Anonymous function specifying the sequence
%
%   `lowdiscrepancy(name)` returns a low discrepancy sequence for the usage
%   as a delay generating function in conjunction with `waveletfilters`.
%   Currently, a kronecker sequence and a digital net are implemented.
%
%   See also: waveletfilters

%   Authors: Nicki Holighaus, Clara Hollomey, Guenther Koliander

definput.keyvals.s = ceil(log2(4096));

[~, kv] = ltfatarghelper({}, definput, varargin);


switch name
    case 'digital'
        input = (0:2^kv.s-1);
        bin_vecs = [fliplr(dec2binary(input))';zeros(1,2^kv.s)];

        temp = tril(ones(kv.s+1));
        temp(3:end,1:end-2) = temp(3:end,1:end-2) - temp(1:end-2,1:end-2);
        C1 = temp;

        out = mod(C1*bin_vecs,2);
        ord = zeros(1,size(out,2));
        for kk = 1:size(out,1)
            ord = ord + out(kk,:).*2^(-kk);
        end    
        delays = @(n,a) a*(mod(ord(n+1)+.5,1)-.5);
    
    case 'kronecker'
        alpha = 1-2/(1+sqrt(5)); % 1-1/(goldenratio) delay sequence
        delays = @(n,a) a*(mod(n*alpha+.5,1)-.5);
    otherwise
        disp('Sequence not yet implemented.');
end
end

function out = dec2binary(int)
    ll = floor(log2(max(int)))+1;
    out = rem(floor(int(:)*pow2(1-ll:0)),2);
end