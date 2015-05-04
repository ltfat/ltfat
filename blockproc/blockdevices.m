function devs = blockdevices()
%BLOCKDEVICES Lists audio devices
%   Usage: devs = blockdevices();
%
%   `blockdevices` lists the available audio input and output devices. The
%   ID can be used in the |block| function to specify which device should
%   be used.
%
%   See also: block

clear playrec;

devs = playrec('getDevices');

fprintf('\nAvailable output devices:\n');

for k=1:length(devs)
    if(devs(k).outputChans)
        fs = sprintf('%d, ',devs(k).supportedSampleRates);
        fs = ['[',fs(1:end-2),']' ];
        fprintf(['ID =%2d: %s (%s) %d chan., latency %d--%d ms,'...
                 ' fs %s\n'], ...
            devs(k).deviceID, devs(k).name, ...
            devs(k).hostAPI, devs(k).outputChans,...
            floor(1000*devs(k).defaultLowOutputLatency),...
            floor(1000*devs(k).defaultHighOutputLatency),...
            fs);

    end
end

fprintf('\nAvailable input devices:\n');

for k=1:length(devs)
    if(devs(k).inputChans)
        fs = sprintf('%d, ',devs(k).supportedSampleRates);
        fs = ['[',fs(1:end-2),']' ];
        fprintf(['ID =%2d: %s (%s) %d chan., latency %d--%d ms,'...
                 ' fs %s\n'], ...
            devs(k).deviceID, devs(k).name, ...
            devs(k).hostAPI, devs(k).inputChans,...
            floor(1000*devs(k).defaultLowInputLatency),...
            floor(1000*devs(k).defaultHighInputLatency),fs);

    end
end
