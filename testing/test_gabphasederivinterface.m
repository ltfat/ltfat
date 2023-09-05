function test_failed = test_gabphasederivinterface
test_failed = 0;

global LTFAT_TEST_TYPE;
tolerance = 1e-7;
if strcmpi(LTFAT_TEST_TYPE,'single')
   tolerance = 1e-3;
end

%Test whether the batch interface does the same thing as indivisual calls

% TODO: test also wrong number of input arguments

disp('---------------test_gabphasederivinterface---------------------');

f = tester_rand(128,1);
a = 4;
M = 8;
g = 'gauss';

c = dgt(f,g,a,M);
algCell = {'dgt','abs','phase','cross'};
algArgs = {{'dgt',f,g,a,M}
    {'abs',abs(c),g,a}
    {'phase',angle(c),a}
    {'cross',f,g,a,M}
    };

phaseconvCell = {'freqinv','timeinv','symphase','relative'};

dtypecombCell = {% {'t','t'},...
%    {'t','f'},...
    {'f','t'},...
    {'t','f','t'},...
    {'tt','t'},...
    {'ff','f'},...
    {'tf','f'},...
    {'tf','f','tt'},...
    {'tf','f','tt','tf'},...
    {'t','f','tt','ff','ft','tf'},...
    {'t','f','t','ff','tt','tt','tt'},...
    };

for algId=1:numel(algArgs)
    alg = algArgs{algId};
    
    for dtypecombId=1:numel(dtypecombCell)
        dtypecomb = dtypecombCell{dtypecombId};
        for phaseconvId=1:numel(phaseconvCell)
            phaseconv=phaseconvCell{phaseconvId};
            
            
            individual = cell(1,numel(dtypecomb));
            for dtypecompId = 1:numel(dtypecomb)
                individual{dtypecompId} = ...
                    gabphasederiv(dtypecomb{dtypecompId},alg{:},phaseconv);
            end
            batch = gabphasederiv(dtypecomb,alg{:},phaseconv);
            
            
            res = sum(cellfun(@(iEl,bEl) norm(iEl(:)-bEl(:)),individual,batch))/numel(batch);
            [test_failed,failstring]=ltfatdiditfail( res,...
                           test_failed,tolerance);

            flagcompString = sprintf('%s, ',dtypecomb{:});
            fprintf('GABPHASEDERIV BATCH EQ alg:%s order:{%s} phaseconv:%s %d %s\n',alg{1},flagcompString(1:end-2),phaseconv,res,failstring);
        end
    end
end





