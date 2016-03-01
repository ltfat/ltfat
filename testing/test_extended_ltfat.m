function test_extended_ltfat()
% This test suite runs extended tests which either:
%   
%    Take long time to finish
%   
%    Their failure can only be checked by inspecting plots
%
%    We do not care whether they work with single precision or not.
%
%    Any that should be run at least before doing the release
%    but do not fit to be included in test_all_ltfat
%

tests_todo = {
    'erbfilters',...
    'fbreassign',...
    'fbwarped_framebounds',...
    'wfilt',...
    'argfirwin',...
    'gabphasederiv',...
    'audfilters',...
    'demos'
};

total_tests_failed=0;
list_of_failed_tests={};

for name = tests_todo
       tmpfailed = feval(['test_',name{1}]);
       if tmpfailed>0
           list_of_failed_tests{end+1} = ['test_',name{1}];
           total_tests_failed = total_tests_failed + tmpfailed;
       end
       
end


