function test_failed=test_purefreq
  
disp(' ===============  TEST_PUREFREQ ================');

disp('--- Used subroutines ---');

which comp_fftreal
which comp_ifftreal

% This script test all the pure frequency routines.

% This is a list or pairs. The two functions in each pair
% will be evaluated with the same parameters, and the output
% should be the same.
ref_funs={{'ref_dcti','dcti'},...
          {'ref_dcti','ref_dcti_1'},...
	  {'ref_dctii','dctii'},...
	  {'ref_dctiii','dctiii'},...
	  {'ref_dctii','ref_dctii_1'},...
	  {'ref_dctiii','ref_dctiii_1'},...	  	  	  
	  {'ref_dctiv','dctiv'},...
	  {'ref_dsti','dsti'},...
	  {'ref_dsti','ref_dsti_1'},...
	  {'ref_dstii','dstii'},...
	  {'ref_dstiii','dstiii'},...
	  {'ref_dstii','ref_dstii_1'},...
	  {'ref_dstiii','ref_dstiii_1'},...
	  {'ref_dstiv','dstiv'},...
	  {'ref_dft','ref_dft_1'},...
	  {'ref_rdft','ref_rdft_1'},...
	  {'ref_irdft','ref_irdft_1'},...
	  };

% This is a list or pairs. The two functions in each pair will be evaluated
% on purely real input with the same parameters, and the output should be
% the same.
ref_realfuns={{'ref_dcti','dcti'},...
          {'ref_dcti','ref_dcti_1'},...
	  {'ref_dctii','dctii'},...
	  {'ref_dctiii','dctiii'},...
	  {'ref_dctii','ref_dctii_1'},...
	  {'ref_dctiii','ref_dctiii_1'},...	  	  	  
	  {'ref_dctiv','dctiv'},...
	  {'ref_dsti','dsti'},...
	  {'ref_dsti','ref_dsti_1'},...
	  {'ref_dstii','dstii'},...
	  {'ref_dstiii','dstiii'},...
	  {'ref_dstii','ref_dstii_1'},...
	  {'ref_dstiii','ref_dstiii_1'},...
	  {'ref_dstiv','dstiv'},...
	  {'ref_dft','ref_dft_1'},...
          {'ref_fftreal','fftreal'},...          
	  {'ref_rdft','ref_rdft_1'},...
	  {'ref_irdft','ref_irdft_1'},...
	  };


  %	  {'ref_dftii','ref_dftii_1'},...
%	  {'ref_dftii','ref_dftii_2'},...
%	  {'ref_idftii','ref_idftii_1'},...
%	  {'ref_dftiv','ref_dftiv_1'},...
%	  {'ref_rdftiii','ref_rdftiii_1'},...
%	  {'ref_irdftiii','ref_irdftiii_1'},...

% Each row is the size of a matrix. All the functions
% mentioned above will be executed with input matrixes
% of sizes mentioned below.
ref_sizes=[1 1;...
	   2 1;...
	   3 2;...
	   4 3;...
	   5 3];
	   
% As ref_funs, these functions should be inverses of each other.
inv_funs={{'dcti','dcti'},...
	  {'dctii','dctiii'},...
	  {'dctiv','dctiv'},...
	  {'ref_dsti','ref_dsti'},...
	  {'dstii','dstiii'},...
	  {'dstiv','dstiv'},...
	  {'dft','idft'},...
	  {'ref_rdft','ref_irdft'},...
	  };

% As ref_funs, these functions should be inverses of each other.
realinv_funs={{'fftreal','ifftreal'}}
  
  %	  {'ref_rdftiii','ref_irdftiii'},...
%	  {'ref_dftii','ref_idftii'},...
%	  {'ref_dftiii','ref_idftiii'},...
%	  {'ref_rdftiii','ref_irdftiii'},...
%	  {'ref_dftiv','ref_idftiv'},...

% Each function in the list should be unitary. They will be tested on
% matrix of sizes correspondin to the first column in ref_sizes 

nrm_funs={'dctii','dctiii','dctiv','ref_dft','ref_rdft'};

% ---------- reference testing --------------------
% Test that the transforms agree on values.

ref_failed=0;

for funpair=ref_funs
  for ii=1:size(ref_sizes,1);
    %a=rand(ref_sizes(ii,:));
    a=crand(ref_sizes(ii,1),ref_sizes(ii,2));

    c1=feval(funpair{1}{1},a);
    c2=feval(funpair{1}{2},a);

    res=norm(c1(:)-c2(:));

    s=sprintf('REF %7s L:%2i W:%2i %0.5g',funpair{1}{2},ref_sizes(ii,1),ref_sizes(ii,2),res);
    disp(s)

    if res>10e-10
      disp('FAILED');
      ref_failed=ref_failed+1;
    end;

  end;
end;


% ---------- real valued reference testing --------------------
% Test that the transforms agree on values.

ref_failed=0;

for funpair=ref_realfuns
  for ii=1:size(ref_sizes,1);
    a=rand(ref_sizes(ii,1),ref_sizes(ii,2));

    c1=feval(funpair{1}{1},a);
    c2=feval(funpair{1}{2},a);

    res=norm(c1(:)-c2(:));

    s=sprintf('REA %7s L:%2i W:%2i %0.5g',funpair{1}{2},ref_sizes(ii,1),ref_sizes(ii,2),res);
    disp(s)

    if res>10e-10
      disp('FAILED');
      ref_failed=ref_failed+1;
    end;

  end;
end;

%------------ inversion testing -----------------
% Test that the transforms are invertible

inv_failed=0;

for funpair=inv_funs
  for ii=1:size(ref_sizes,1);
    %a=rand(ref_sizes(ii,:));
    a=crand(ref_sizes(ii,1),ref_sizes(ii,2));

    ar=feval(funpair{1}{2},feval(funpair{1}{1},a));

    res=norm(a(:)-ar(:));

    s=sprintf('INV %7s L:%2i W:%2i %0.5g',funpair{1}{1},ref_sizes(ii,1),ref_sizes(ii,2),res);
    disp(s)

    if res>10e-10
      disp('FAILED');
      inv_failed=inv_failed+1;
    end;

  end;
end;


%----------- normalization test ----------------
% Test that the transforms are orthonormal

nrm_failed=0;

for funname=nrm_funs

  for ii=1:size(ref_sizes,1);

    L=ref_sizes(ii,1);

    F=feval(funname{1},eye(L));

    res=norm(F*F'-eye(L));

    s=sprintf('NRM %7s L:%2i %0.5g',funname{1},ref_sizes(ii,1),res);
    disp(s)

    if res>10e-10
      disp('FAILED');
      nrm_failed=nrm_failed+1;
    end;
  end;
end;

%------------ test fftreal inverseion ----------------
% Test that the transforms are invertible

realinv_failed=0;

for funpair=realinv_funs
  for ii=1:size(ref_sizes,1);
    a=rand(ref_sizes(ii,1),ref_sizes(ii,2));

    ar=ifftreal(fftreal(a),ref_sizes(ii,1));

    res=norm(a(:)-ar(:));

    s=sprintf('RIN %7s L:%2i W:%2i %0.5g',funpair{1}{1},ref_sizes(ii,1),ref_sizes(ii,2),res);
    disp(s)

    if res>10e-10
      disp('FAILED');
      inv_failed=inv_failed+1;
    end;

  end;
end;


test_failed=ref_failed+inv_failed+nrm_failed;


%OLDFORMAT
