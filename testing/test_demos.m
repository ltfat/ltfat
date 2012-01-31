function test_failed=test_demos
%TEST_DEMOS  Test if all the demos runs without errors.
  
  test_failed=0;
  
  s=dir([ltfatbasepath,filesep,'demos',filesep,'demo*.m']);

  for ii=1:numel(s)
     filename = s(ii).name;
     
     disp(filename);
     
     eval(filename(1:end-2));
     
  end;

