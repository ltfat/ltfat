function devs = blockdevices()
%BLOCKDEVICES Lists audio devices
%   Usage: devs = blockdevices();
%
%   `blockdevices` lists the available audio input and output devices. The
%   ID can be used in the |block| function to specify which device should
%   be used.
%
%   See also: block

devs = playrec('getDevices');

fprintf('\nAvailable output devices:\n');

for k=1:length(devs)
    if(devs(k).outputChans)
        fprintf('ID =%2d: %s (%s) %d channels, default latency %d--%d ms\n', ...
            devs(k).deviceID, devs(k).name, ...
            devs(k).hostAPI, devs(k).outputChans,...
            floor(1000*devs(k).defaultLowOutputLatency),...
            floor(1000*devs(k).defaultHighOutputLatency));

    end
end

fprintf('\nAvailable input devices:\n');

for k=1:length(devs)
    if(devs(k).inputChans)
        fprintf('ID =%2d: %s (%s) %d channels,  default latency %d--%d ms\n', ...
            devs(k).deviceID, devs(k).name, ...
            devs(k).hostAPI, devs(k).inputChans,...
            floor(1000*devs(k).defaultLowInputLatency),...
            floor(1000*devs(k).defaultHighInputLatency));

    end
end