function devs = blockdevices()
%BLOCKDEVICES Lists audio devices
%
% Function lists available audio input and output devices. ID can be used
% in the |block| function to specify which device should be used.

devs = playrec('getDevices');

fprintf('\nAvailable output devices:\n');

for k=1:length(devs)
    if(devs(k).outputChans)
        fprintf('ID =%2d: %s (%s) %d channels\n', ...
            devs(k).deviceID, devs(k).name, ...
            devs(k).hostAPI, devs(k).outputChans);

    end
end

fprintf('\nAvailable input devices:\n');

for k=1:length(devs)
    if(devs(k).inputChans)
        fprintf('ID =%2d: %s (%s) %d channels\n', ...
            devs(k).deviceID, devs(k).name, ...
            devs(k).hostAPI, devs(k).inputChans);

    end
end