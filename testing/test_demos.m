function test_failed=test_demos
%TEST_DEMOS  Test if all the demos runs without errors.
  
  test_failed=0;
  
  s=dir([ltfatbasepath,filesep,'demos',filesep,'demo_*.m']);

  for ii=1:numel(s)
     filename = s(ii).name;
     
     disp(filename);
     
     % The demo is run in separate function to avoid 
     % variable name clash
     rundemo(filename(1:end-2));
     
  end
  
  
function rundemo(demoname)
   close all;
   eval(demoname);
   clear all; 

