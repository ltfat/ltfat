function test_failed = test_nonu2ufilterbank
% This tests equality of coefficients of a non uniform and
% an identical uniform filterbanks
test_failed=0;

disp('-------------TEST_NONU2UFILTERBANK-------------');

for ftypeCell = {'time','freq'}
    ftype = ftypeCell{1};
    disp(sprintf('--------------- %s ------------',ftype));
for M=[1,6];

for aCell = {(1:M)',(M:-1:1)',ones(M,1), M*ones(M,1)}
a = aCell{1};
N=10;

L=filterbanklength(1,a)*N;
for W = 1:3

g=cell(1,M);
if strcmp(ftype,'time')
for ii=1:M
  g{ii}=tester_crand(L/2,1);
end;
elseif strcmp(ftype,'freq')
 for ii=1:M
  g{ii}=struct('H',tester_crand(L,1),'foff',0,'L',L);
end;   
end

f = tester_crand(L,W);

% Target coefficients
c_non = filterbank(f,g,a);

% Do a uniform filterbank
[gu,au,pk] = nonu2ufilterbank(g,a);

% Identical cell-array
c_u = filterbank(f,gu,au);

% Identical matrix
c_uu = ufilterbank(f,gu,au);

res=0;
for m=1:M
    res=res+norm(c_u{m}-squeeze(c_uu(:,m,:)));  
end;
[test_failed,fail]=ltfatdiditfail(res,test_failed);
fprintf('COEFEQ L:%3i,W:%3i,M:%3i, %0.5g %s\n',L,W,M,res,fail);

% Convert back to nonuniform format
c = u2nonucfmt(c_u,pk);

res = norm(cell2mat(c_non) - cell2mat(c));
[test_failed,fail]=ltfatdiditfail(res,test_failed);
fprintf('FORMAT CELL L:%3i,W:%3i,M:%3i, %0.5g %s\n',L,W,M,res,fail);

% Convert back to nonuniform format
c = u2nonucfmt(c_uu,pk);

res = norm(cell2mat(c_non) - cell2mat(c));
[test_failed,fail]=ltfatdiditfail(res,test_failed);
fprintf('FORMAT MAT L:%3i,W:%3i,M:%3i, %0.5g %s\n',L,W,M,res,fail);


% Convert 
c_uuu = nonu2ucfmt(c_non,pk);

res = norm(cell2mat(c_uuu) - cell2mat(c_u));
[test_failed,fail]=ltfatdiditfail(res,test_failed);
fprintf('FORMATBACK L:%3i,W:%3i,M:%3i, %0.5g %s\n',L,W,M,res,fail);


end
end
end
end
