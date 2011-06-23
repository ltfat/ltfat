function [test_failed,fail]=ltfatdiditfail(res,test_failed);
%LTFATDIDITFAIL  Did a test fail
%
%  [test_fail,fail]=LTFATDIDITFAIL(res,test_fail) updates test_fail if
%  res is above threshhold and outputs the word FAIL in the variable
%  fail. Use only in testing scripts.
  
fail='';
if (abs(res)>10e-10) || isnan(res)
  fail='FAILED';
  test_failed=test_failed+1;
end;