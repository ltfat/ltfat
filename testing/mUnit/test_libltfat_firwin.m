function test_failed = test_libltfat_firwin(varargin)
test_failed = 0;

fprintf(' ===============  %s ================ \n',upper(mfilename));

definput.flags.complexity={'double','single'};
[flags]=ltfatarghelper({},definput,varargin);
dataPtr = [flags.complexity, 'Ptr'];

[~,~,enuminfo]=libltfatprotofile;
Cenumnorms = enuminfo.LTFAT_FIRWIN;

d = arg_firwin;
wins = d.flags.wintype;
wins(strcmp('truncgauss',wins)) = [];

wins = [wins 'truncgauss01'];

names =fieldnames(Cenumnorms);

libwins = {};
for nameId = 1:numel(wins)
libwins{end+1} = Cenumnorms.(names{strcmpi(['ltfat_', wins{nameId}],names)});
end

Larr = [1,9,10,11,110,111];

for do_complex = 0:1
    complexstring = '';
    if do_complex, complexstring = 'complex'; end
    funname = makelibraryname('firwin',flags.complexity,do_complex);
    for L = Larr
        for nId = 1:numel(wins)
            win = wins{nId};
            lwin = libwins{nId};

            if do_complex
                z = cast((1:L)' + i*(1:L)',flags.complexity);
                zi = complex2interleaved(z);
                ziPtr = libpointer(dataPtr,zi);
            else
                z = cast((1:L)',flags.complexity);
                zi = z;
                ziPtr = libpointer(dataPtr,zi);
            end

            status = calllib('libltfat',funname,lwin,L,ziPtr);

            startswith = 'truncgauss';
            if regexpi(win,['^',startswith])
                percent = 1;
                if numel(win) > numel(startswith)
                    percent = str2double(win(numel(startswith)+1:end));
                end

                trueres = long2fir(pgauss(10*L,'width',L,'atheight',percent/100,'inf'),L);
                if rem(L,2) == 0
                    trueres(end/2+1) = 0; %
                end
            else
                trueres = firwin(win,L);
            end
            
            if do_complex
                res = norm(trueres - interleaved2complex(ziPtr.Value));
            else
                res = norm(trueres - ziPtr.Value);
            end

            [test_failed,fail]=ltfatdiditfail(res+status,test_failed);
            fprintf(['FIRWIN L:%3i, %s %s %s %s %s\n'],L,win,flags.complexity,complexstring,ltfatstatusstring(status),fail);
        end
     end
end


