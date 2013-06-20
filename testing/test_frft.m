function test_failed=test_frft



disp(' ===============  TEST_FRFT ===========');

Lr=[9,10,11,12];

test_failed=0;

% Test the hermite functions and discrete frft
for ii=1:length(Lr)
	L=Lr(ii);
	F=fft(eye(L))/sqrt(L);

	% check if hermite functions are eigenfunctions of F
	V=hermbasis(L,4);
	res=norm(abs(F*V)-abs(V));
	[test_failed,fail]=ltfatdiditfail(res,test_failed);          
        s=fprintf('HERMBASIS L:%3i %0.5g %s\n',L,res,fail);

	% Frft of order 1 becomes ordinary DFT
	f1=tester_crand(L,1);
	f2=tester_crand(1,L);

	p=4;
	frf1=dfracft(f1,1,[],p);
	frf2=dfracft(f2,1,2,p);
	res=norm(F*f1-frf1);
	[test_failed,fail]=ltfatdiditfail(res,test_failed);          
        s=fprintf('DFRACFT  L:%3i, %0.5g %s\n',L,res,fail);
	res=norm(f2*F-frf2);
	[test_failed,fail]=ltfatdiditfail(res,test_failed);          
        s=fprintf('DFRACFT  L:%3i, %0.5g %s\n',L,res,fail);

	frf1=dfracft(f1,1);
	frf2=dfracft(f2,1,2);
	res=norm(F*f1-frf1);
	[test_failed,fail]=ltfatdiditfail(res,test_failed);          
        s=fprintf('DFRACFT  L:%3i %0.5g %s\n',L,res,fail);
	res=norm(f2*F-frf2);
	[test_failed,fail]=ltfatdiditfail(res,test_failed);          
        s=fprintf('DFRACFT  L:%3i %0.5g %s\n',L,res,fail);


end

