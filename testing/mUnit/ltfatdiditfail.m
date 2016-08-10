function [test_failed,fail]=ltfatdiditfail(res,test_failed,tolerance);
%LTFATDIDITFAIL  Did a test fail
%
%  [test_fail,fail]=LTFATDIDITFAIL(res,test_fail) updates test_fail if
%  res is above threshhold and outputs the word FAIL in the variable
%  fail. Use only in testing scripts.
if nargin<3
  tolerance=1e-10;
  if isa(res,'single')
      tolerance=5e-4;
  end
end;
  
fail='';
if (abs(res)>tolerance) || isnan(res)
  fail='FAILED';
  test_failed=test_failed+1;
end;

