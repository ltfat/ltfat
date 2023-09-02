% LTFAT - Filterbanks
%
%  Peter L. Søndergaard, 2011 - 2023.
%
%  Transforms and basic routines
%    FILTERBANK             - Filter bank
%    UFILTERBANK            - Uniform Filter bank
%    IFILTERBANK            - Inverse normal/uniform filter bank
%    IFILTERBANKITER        - Iteratively inverse filter bank 
%    FILTERBANKWIN          - Evaluate filter bank window
%    FILTERBANKLENGTH       - Length of filter bank to expand signal
%    FILTERBANKLENGTHCOEF   - Length of filter bank to expand coefficients
%
%  Filter generators
%    CQTFILTERS             - Logarithmically spaced filters
%    ERBFILTERS             - ERB-spaced filters
%    WARPEDFILTERS          - Frequency-warped filters 
%    AUDFILTERS             - Filters based on auditory scales
%    GABFILTERS             - Linearly spaced Gabor filters
%    WAVELETFILTERS         - Waveletfilters
%  
%  Window construction and bounds
%    FILTERBANKDUAL         - Canonical dual filters
%    FILTERBANKTIGHT        - Canonical tight filters
%    FILTERBANKREALDUAL     - Canonical dual filters for real-valued signals
%    FILTERBANKREALTIGHT    - Canonical tight filters for real-valued signals
%    FILTERBANKBOUNDS       - Frame bounds of filter bank
%    FILTERBANKREALBOUNDS   - Frame bounds of filter bank for real-valued signals
%    FILTERBANKRESPONSE     - Total frequency response (a frame property)
%
%  Auxilary
%    FILTERBANKFREQZ        - Frequency responses of filters
%    FILTERBANKSCALE        - Scaling and normalization of filters
%    NONU2UFILTERBANK       - Non-uni. to uniform filter bank transformation
%    U2NONUCFMT             - Change format of coefficients
%    NONU2UCFMT             - Change format of coefficients back
%
%  Plots
%    PLOTFILTERBANK         - Plot normal/uniform filter bank coefficients
%
%  Reassignment and phase gradient
%    FILTERBANKPHASEGRAD      - Instantaneous time/frequency from signal
%    FILTERBANKREASSIGN       - Reassign filterbank spectrogram
%    FILTERBANKSYNCHROSQUEEZE - Synchrosqueeze filterbank spectrogram  
%
%  Phase reconstruction
%    FILTERBANKCONSTPHASE     - Construct suitable phase from the coefficient magnitude
%    
%
%  For help, bug reports, suggestions etc. please visit 
%  http://github.com/ltfat/ltfat/issues

