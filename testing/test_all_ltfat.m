tests_todo={'dgt','dwilt','wmdct',...
            'dgt_fb',...
            'purefreq','zak',...
            'dgt2','dwilt2','wmdct2',...
            'firwin',...
            'gabmulappr',...
            'spread',...
            'pconv','involute',...
            'signals','realout','windrivers',...
            'nsdgt','filterbank',...
            'pgauss','pfilt',...
            'rangecompress'};

% Testing of pbspline has been removed, as it causes too much trouble.

total_tests_failed=0;
list_of_failed_tests={};

for ii=1:length(tests_todo)
  test_failed=feval(['test_',tests_todo{ii}]);
  total_tests_failed=total_tests_failed+test_failed;
  if test_failed>0
    list_of_failed_tests{end+1}=['test_',tests_todo{ii}];
  end;
end;

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
