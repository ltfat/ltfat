if 1
  % This is the setup used in the paper
  a=40;
  M=60; 
  L=a*M;
  W=4;
  nrep=50;
else
  a=50;
  M=200; 
  L=a*M; 
  W=4;
  nrep=50;  
  
end;

system('rm crossover.log');
%for gl=M:M:20*M
for gl=M:M:16*M
  s=sprintf('./time_dgt_fb %i %i %i %i %i %i >> crossover.log\n',a,M,L,W,gl,nrep);
  
  disp(s);
  system(s);
end;

s=sprintf('./time_dgt_fac %i %i %i %i %i > crossover.ref\n',a,M,L,W,nrep);

disp(s);
system(s);

system('rm crossover_real.log');
%for gl=M:M:20*M
for gl=M:M:16*M
  s=sprintf('./time_dgtreal_fb %i %i %i %i %i %i >> crossover_real.log\n',a,M,L,W,gl,nrep);
  
  disp(s);
  system(s);
end;

s=sprintf('./time_dgtreal_fac %i %i %i %i %i > crossover_real.ref\n',a,M,L,W,nrep);

disp(s);
system(s);

