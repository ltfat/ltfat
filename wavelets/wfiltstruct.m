function wf = wfiltstruct(type)

% The single filter structure definition

% Filter type
% Possible values: Enum style: 'FIR','FREQ' 
wf.type = type;

% Filter delay (position of the zero index in the implulse response)
wf.d = [];

% Impulse response
wf.h = [];


