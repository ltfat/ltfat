function f=comp_ifwt_undec(c,g,J)

flen = length(g{1});
gpad = cell(J,length(g));
[L,W] = size(c{end});
upfLen = @(j) max([flen, 2^j*flen-(2^j-1)]);


       for hidx=1:length(g)
            for j=0:J-1
                flenTemp = upfLen(j);
                gpad{J-j,hidx}=zeros(L,1);
                range = 1:flenTemp;
                gpad{J-j,hidx}(range) = ups(g{hidx}(:),2^j,1);
 
                % circshift to compensate for the delay
                gpad{J-j,hidx} = circshift(gpad{J-j,hidx},-(2^j*flen)/2 +2^j);
           end
        end

f = zeros(L,W);
  
  
for w=1:W  
  tempa = c{1}(:,w);
    for j=1:J
       tempa = ( pconv(tempa, gpad{j,1}) + pconv(c{j+1}(:,w), gpad{j,2}) )/2;
    end
    f(:,w) = tempa; 
end