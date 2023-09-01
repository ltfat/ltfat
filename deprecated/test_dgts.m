% This is currently the tester for the
% new TF-transforms.
%
% It does not do multiwindow.

%Lr=[24,144,108,144,24,135,35,77];
%ar=[ 4,  9,  9, 12, 6,  9, 5, 7];
%Mr=[ 6, 16, 12, 24, 8,  9, 7,11];

Lr=[24,20,15,12];
ar=[ 4, 4, 3, 3];
Mr=[ 6, 5, 5, 4];

Wmax=2;

ref_funs={{'ref_rdgt','ref_rdgt_1'},...
	  {'ref_rdgt3','ref_rdgt3_1'},...
	  {'ref_edgtii','ref_edgtii_1'},...
	  };

refinv_funs={{'ref_irdgt','ref_irdgt_1'},...
	     {'ref_irdgt3','ref_irdgt3_1'},...
	     {'ref_iedgtii','ref_iedgtii_1'},...
	     };

inv_funs={{'ref_rdgt','ref_irdgt'},...
	  {'ref_rdgt3','ref_irdgt3'},...
	  {'ref_edgtii','ref_iedgtii'},...
	  };

% ---------- reference testing --------------------
% Test that the transforms agree on values.

for funpair=ref_funs
    
    for ii=1:length(Lr);
      
      for W=1:Wmax

	L=Lr(ii);
	
	M=Mr(ii);
	a=ar(ii);
	
	%b=L/M;
	N=L/a;
    
	f=tester_rand(L,W);
	g=ref_win(funpair{1}{1},'test',L,a,M);
    
	c1=feval(funpair{1}{1},f,g,a,M);
	c2=feval(funpair{1}{2},f,g,a,M);
	
	res=norm(c1(:)-c2(:));
    
	s=sprintf('REF %10s L:%3i W:%2i a:%3i M:%3i %0.5g',funpair{1}{2},L,W,a,M,res);
	disp(s)
    
	%if res>1e-8
	  %	disp('Reference comparion failed:')
	  %	[L, M, N, a, b, c, d]
	  %      end;
	
      end;
      
    end;
end;

% ---------- reference inverse function testing ---------------
% Test that the transforms agree on values.

for funpair=refinv_funs
    
    for ii=1:length(Lr);
      
      for W=1:Wmax

	L=Lr(ii);
	
	M=Mr(ii);
	a=ar(ii);
	
	%b=L/M;
	N=L/a;
    
	c=tester_rand(M*N,W);
	g=ref_win(funpair{1}{1},'test',L,a,M);
    
	f1=feval(funpair{1}{1},c,g,a,M);
	f2=feval(funpair{1}{2},c,g,a,M);
	
	res=norm(f1(:)-f2(:));
    
	s=sprintf('REF %10s L:%3i W:%2i a:%3i M:%3i %0.5g',funpair{1}{2},L,W,a,M,res);
	disp(s)
    
	%if res>1e-8
	  %	disp('Reference comparion failed:')
	  %	[L, M, N, a, b, c, d]
	  %      end;
	
      end;
      
    end;
end;

%------------ inversion testing -----------------
% Test that the transforms are invertible

for funpair=inv_funs
    
    for ii=1:length(Lr);
      
      for W=1:Wmax

	L=Lr(ii);
	
	M=Mr(ii);
	a=ar(ii);
	
	%b=L/M;
	N=L/a;
    
	f=tester_rand(L,W);
	g=ref_win(funpair{1}{1},'test',L,a,M);
	gamma=ref_tgabdual(funpair{1}{1},g,L,a,M);
    
	c=feval(funpair{1}{1},f,g,a,M);
	fr=feval(funpair{1}{2},c,gamma,a,M);
	
	res=norm(f(:)-fr(:));
    
	s=sprintf('INV %10s L:%3i W:%2i a:%3i M:%3i %0.5g',funpair{1}{1},L,W,a,M,res);
	disp(s)
    
	%if res>1e-8
	  %	disp('Reference comparion failed:')
	  %	[L, M, N, a, b, c, d]
	  %      end;
	
      end;
      
    end;

end;




