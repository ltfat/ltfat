function test_failed = test_argfirwin
test_failed = 0;

disp('---------Testing arg_firwin---------');
% This test checks wheter all options from arg_firwin are actually
% treated in firwin.
% This is necessary as firwin itself does not use arg_firwin
% Also tests for duplicates


% Get window types
g = getfield(getfield(arg_firwin,'flags'),'wintype');

for ii=1:numel(g)
    
        tryfailed = 0;
        fprintf('Testing %12s:',g{ii});
        try
            h = firwin(g{ii},64);
        catch
            test_failed = test_failed + 1;
            tryfailed = 1;
            fprintf('FAILED');
        end
        if ~tryfailed
            fprintf('SUCCESS');
        end
        
        if(sum(strcmp(g{ii},g))>1)
            fprintf(' has DUPLICATES');
        end
        fprintf('\n');
end
