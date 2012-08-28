function [s0,s1,X] = shearfind(a,b,s,N)

    if nargin < 4 
        error('Too few input arguments');
    end

    if a>N || b>N || N/a ~= round(N/a) || N/b ~= round(N/b)
        error('a and b must be divisors of N');
    end

    if s >= b || s < 0
        error('Choose an s < b');
    elseif s == 0
        s0 = 0; s1 = 0; 
        X = b;
        return;
    end

    [Nabfac,sfac,lenNabfac] = lattfac(a,b,s,N);

    %lenNabfac = size(Nabfac,2);

    if s/a == round(s/a)
        if s/a <= b/2
            s1 = -s/a;
        else 
            s1 = b-s/a;
        end
        s0 = 0;
        %alpha = 0;
        X = b;
    elseif ones(1,lenNabfac) == (Nabfac(3,:) < sfac(2,1:end-1))
        s1 = 0;
        [X,alpha,temp] = gcd(s,b);
        if alpha < 0
            alpha = b/X + alpha;
        end
        s0 = alpha*a/X;    
    else
        s1fac = (Nabfac(3,:) == sfac(2,1:end-1)).*(Nabfac(3,:) < Nabfac(4,:));
        s1 = prod(Nabfac(1,:).^s1fac);

        if s1*a/b == round(s1*a/b) 
            s1 = 0; 
        else 
           B = prod(Nabfac(1,:).^max(Nabfac(4,:)-Nabfac(3,:),0));
           if s1 > B/2
                s1 = s1-B;
           end   
        end

        [X,alpha,temp] = gcd(s1*a+s,b);
        if alpha < 0
            alpha = b/X + alpha;
        end
        tempX = factor(X);
        tempalph = factor(alpha);

        Xfac = zeros(1,length(lenNabfac));
        alphfac = zeros(1,length(lenNabfac)+1);

        for kk = 1:lenNabfac
            Xfac(kk) = sum(tempX == Nabfac(1,kk));
            tempX = tempX(tempX ~= Nabfac(1,kk));
            alphfac(kk) = sum(tempalph == Nabfac(1,kk));
            tempalph = tempalph(tempalph ~= Nabfac(1,kk));
        end

        alphfac(lenNabfac+1) = prod(tempalph);

        s0fac = [Nabfac(3,:)+min(alphfac(1:end-1),Nabfac(4,:)-Xfac)-Xfac,0];
        pwrs = max(Nabfac(4,:)-Xfac-alphfac(1:end-1),0);
        pwrs2 = max(-Nabfac(4,:)+Xfac+alphfac(1:end-1),0);

        K = ceil(alphfac(end).*prod(Nabfac(1,:).^(pwrs2-pwrs))-.5);

        s0fac(end) = K*prod(Nabfac(1,:).^pwrs) - alphfac(end).*prod(Nabfac(1,:).^pwrs2);

        s0 = prod(Nabfac(1,:).^s0fac(1:end-1))*s0fac(end);

        if s0*X^2/(a*b) == round(s0*X^2/(a*b)) 
            s0 = 0; 
        end
    end
    
end

function [Nabfac,sfac,lenNabfac] = lattfac(a,b,s,N)    

    tempN = factor(N);
    tempa = factor(a);
    
    if tempa == 1
        tempa = [];
    end
    tempb = factor(b);
    if tempb == 1
        tempb = [];
    end

    Nabfac = unique(tempN);
    lenNabfac = length(Nabfac);
    Nabfac = [Nabfac;zeros(3,lenNabfac)];

    for kk = 1:lenNabfac
       Nabfac(2,kk) = sum(tempN == Nabfac(1,kk));
       tempN = tempN(tempN ~= Nabfac(1,kk));
       Nabfac(3,kk) = sum(tempa == Nabfac(1,kk));
       tempa = tempa(tempa ~= Nabfac(1,kk));
       Nabfac(4,kk) = sum(tempb == Nabfac(1,kk));
       tempb = tempb(tempb ~= Nabfac(1,kk));
    end

    if isempty(tempa) == 0 || isempty(tempb) == 0
        error('a and b must be divisors of N');
    end

    if s*N/(a*b) ~= round(s*N/(a*b));
        error('s must be a multiple of a*b/N');
    end

    temps = factor(s);

    sfac = [Nabfac(1,:),0;zeros(1,lenNabfac+1)];

    for kk = 1:lenNabfac
       sfac(2,kk) = sum(temps == sfac(1,kk));
       temps = temps(temps ~= sfac(1,kk));
    end

    sfac(:,lenNabfac+1) = [prod(temps);1];

end
