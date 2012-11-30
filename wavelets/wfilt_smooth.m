function [H, G, a] = wfilt_smooth(N,varargin)

H = cell(2,1);



H{1}= real(poly(hroots));
% normalize to norm(H{1})=sqrt(2)
H{1}= sqrt(2)*H{1}/sum(H{1});
% QMF modulation low-pass -> highpass
H{2}= (-1).^(1:flen).*H{1}(end:-1:1);


if(nargout>1)
   % Building reconstruction filterbank
   G = cell(2,1); 
   % flip
   G{1} = H{1}(end:-1:1);
   % modulation
   G{2} = -(-1).^(1:flen).*H{1};

   if(nargout>2)
       a = [2;2];
   end
   if(J>1)
        % make J level multirate identity filterbank
        G = makeMultirateIdentity(G{1},G{2},J);
       if(nargout>2)
         a = zeros(J+1,1);
         a(1) = 2^J;
         for j=1:J a(end+1-j)= 2^j; end;
       end
   end
end

if(J>1)
    % make J level multirate identity filterbank
    H = makeMultirateIdentity(H{1},H{2},J);
end





