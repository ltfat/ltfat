function f=comp_ifwt(c,g,J)

flen = length(g{1});
gpad = cell(J,length(g));
[Lhalf,W] = size(c{end});
L = 2*Lhalf;

       for hidx=1:length(g)
            for j=0:J-1
               gpad{j+1,hidx}=zeros(2*length(c{j+2}(:,1)),1);
               gpad{j+1,hidx}(1:flen) = g{hidx}(:);
               gpad{j+1,hidx} = circshift(gpad{j+1,hidx},-ceil(flen/2));
           end
        end

f = zeros(L,W);
  
  
for w=1:W  
  tempa = c{1}(:,w);
    for j=1:J
       tempa = pconv(ups(tempa,2,2), gpad{j,1}) + pconv(ups(c{j+1}(:,w),2,2), gpad{j,2});
    end
    f(:,w) = tempa;
end