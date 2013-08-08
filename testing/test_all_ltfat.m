global LTFAT_TEST_TYPE;

tests_todo={
    'dgt','dwilt','wmdct',...
    'dgt_fb','multiwin',...
    'purefreq','zak',...
    'gabmulappr',...
    'dgt2','dwilt2','wmdct2',...
    'firwin',...
    'spread','thresh',...
    'pconv','involute',...
    'signals','realout','windrivers',...
    'nsdgt','filterbank',...
    'pgauss','pfilt',...
    'rangecompress',...
    'gabmuleigs','phaselock',...
    'fwtpr','ufwtpr','wfbtpr','uwfbtpr','wpfbtpr','uwpfbtpr','fwt2',...
    'wfbt2filterbank','gga',...
    'frames', 'frft'
           };



% Testing of pbspline has been removed, as it causes too much trouble.

total_tests_failed=0;
list_of_failed_tests={};

precarray={'double','single'};
for precidx=1:numel(precarray)
    prec=precarray{precidx};
    LTFAT_TEST_TYPE=prec;

    for ii=1:length(tests_todo)
        test_failed=feval(['test_',tests_todo{ii}]);
        total_tests_failed=total_tests_failed+test_failed;
        if test_failed>0
            list_of_failed_tests{end+1}=['test_',tests_todo{ii},' ',prec];
        end;
    end;
end;
LTFAT_TEST_TYPE=precarray{1};

disp(' ');
if total_tests_failed==0
  disp('ALL TESTS PASSED');
else
  s=sprintf('%i TESTS FAILED',total_tests_failed);
  disp(s);
  disp('The following test scripts contained failed tests');
  for ii=1:length(list_of_failed_tests)
    disp(['   ',list_of_failed_tests{ii}]);
  end;
end;



