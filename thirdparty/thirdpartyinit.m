status=2;

p=[bp,'thirdparty',filesep];

% PolygonClip is not distributed with the Octave-forge package
s=[p,'PolygonClip'];
if exist(s)
    addpath(s);
end;

addpath([p,'Playrec']);


