function test_failed = test_blocprocoffline()


% Scanario 1) Block reading from a wav and block writing to a wav
%             
%             
inName = 'test_in.wav';
outName = 'test_out.wav';
f = 2*rand(44100,1)*0.9-1;
wavwsave(f,44100,inName);

fs=block(inName,'offline','outfile',outName);

flag = 1;
while flag
   [fb,flag] = blockread();
   blockwrite(fb/2);
end




% Scanario 2) blockwrite from vector to wav
%             
%             


f = gspi;
f2 = 2*rand(numel(f),1)*0.9-1;

fs = block([f,f2],'fs',44100,'offline','outfile',outName);

flag = 1;
while flag
   [fb,flag] = blockread(44100);
   blockwrite(fb/2);
end




delete(inName);
delete(outName);
