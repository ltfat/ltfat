function [s,fs]=bat()
%BAT  Load the 'bat' test signal
%   Usage:  s=bat;
%
%   `bat` loads the 'bat' signal. It is a 400 samples long recording
%   of a bat chirp sampled with a sampling period of 7 microseconds.
%   This gives a sampling rate of 143 kHz.
%
%   `[sig,fs]=bat` additionally returns the sampling frequency *fs*.
%
%   The signal can be obtained from
%   `<http://dsp.rice.edu/software/bat-echolocation-chirp>`_
%
%   Please acknowledge use of this data in publications as follows:
%
%     The author wishes to thank Curtis Condon, Ken White, and Al Feng of
%     the Beckman Institute of the University of Illinois for the bat data
%     and for permission to use it in this paper.
%
%   Examples:
%   ---------
%
%   Plot of 'bat' in the time-domain:::
%
%     plot((1:400)/143000,bat);
%     xlabel('Time (seconds)');
%     ylabel('Amplitude');
%
%   Plot of 'bat' in the frequency-domain:::
%
%     plotfftreal(fftreal(bat),143000,90);
%
%   Plot of 'bat' in the time-frequency-domain:::
%
%     sgram(bat,143000,90);
%
%   See also:  batmask

%   AUTHOR : Peter L. Søndergaard
%   TESTING: TEST_SIGNALS
%   REFERENCE: OK
  
if nargin>0
  error('This function does not take input arguments.')
end;

f=mfilename('fullpath');

s=load('-ascii',[f,'.asc']);
fs=143000;
