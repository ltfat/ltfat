function c = comp_dct(f,type)

[L,W] = size(f);


switch type
   case 1
      c=zeros(L,W,assert_classname(f));
      f2=[f;flipud(f(2:L-1,:))]/sqrt(2);
      f2(1,:)=f2(1,:)*sqrt(2);
      f2(L,:)=f2(L,:)*sqrt(2);
  
      % Do DFT.
      s1=fft(f2)/sqrt(2*L-2);

      % This could be done by a repmat instead.
      for w=1:W
         c(:,w)=s1(1:L,w)+[0;s1(2*L-2:-1:L+1,w);0];
      end;

      c(2:L-1,:)=c(2:L-1,:)/sqrt(2);
   case 2
      c=zeros(L,W,assert_classname(f));
      m1=1/sqrt(2)*exp(-(0:L-1)*pi*i/(2*L)).';
      m1(1)=1;

      m2=1/sqrt(2)*exp((1:L-1)*pi*i/(2*L)).';

      s1=fft([f;flipud(f)]);

      % This could be done by a repmat instead.
      for w=1:W
         c(:,w)=s1(1:L,w).*m1+[0;s1(2*L:-1:L+2,w).*m2];
      end;

      c=c/sqrt(L)/2;
   case 3
      c=zeros(2*L,W,assert_classname(f));

      m1=1/sqrt(2)*exp(-(0:L-1)*pi*i/(2*L)).';
      m1(1)=1;
  
      m2=1/sqrt(2)*exp((L-1:-1:1)*pi*i/(2*L)).';

      for w=1:W
         c(:,w)=[m1.*f(:,w);0;m2.*f(L:-1:2,w)];
      end;

      c=fft(c)/sqrt(L);

      c=c(1:L,:);
   case 4
      s1=zeros(2*L,W,assert_classname(f));
      c=zeros(L,W,assert_classname(f));

      m1=1/sqrt(2)*exp(-(0:L-1)*pi*i/(2*L)).';
      m2=1/sqrt(2)*exp((1:L)*pi*i/(2*L)).';

      for w=1:W
        s1(:,w)=[m1.*f(:,w);flipud(m2).*f(L:-1:1,w)];
      end;
  
      s1=exp(-pi*i/(4*L))*fft(s1)/sqrt(2*L);

      % This could be done by a repmat instead.
      for w=1:W
        c(:,w)=s1(1:L,w).*m1+s1(2*L:-1:L+1,w).*m2;
      end;
   otherwise
      error('%s: Not supported type.',upper(mfilename));
end


if isreal(f)
   c=real(c);
end;