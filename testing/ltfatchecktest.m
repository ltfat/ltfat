function [test_failed]=ltfatchecktest(res,outstr,test_failed,mode,tolerance);
%LTFATCHECKTEST  Did a test fail, new method
%
%  `[test_fail,fail]=ltfatchecktest(res,test_fail)` updates *test_fail* if
%  *res* is above threshhold and outputs the word FAIL in the variable
%  fail. Use only in testing scripts.
%
%  mode = 0 prints all results
%
%  mode = 1 is quiet mode to spot failures

if nargin<5
  tolerance=1e-10;
end;
  
fail=0;
if (abs(res)>tolerance) || isnan(res)
  fail=1;
end;

if (mode==0) || (fail==1)
    if (fail==1)
        outstr=[outstr,' FAILED'];
    end;
    disp(outstr);
end;

test_failed=test_failed+fail;

