% LTFAT - Filterbanks
%
%  Peter L. Søndergaard, 2011 - 2013
%
%  Transforms and basic routines
%    FILTERBANK             - Filter bank
%    UFILTERBANK            - Uniform Filter bank
%    IFILTERBANK            - Inverse normal/uniform filter bank
%    FILTERBANKWIN          - Evaluate filterbank window
%    FILTERBANKLENGTH       - Length of filterbank to expand signal
%    FILTERBANKLENGTHCOEF   - Length of filterbank to expand coefficients
%
%  Auditory inspired filterbanks
%    CQT                    - Constant Q transform
%    ICQT                   - Inverse constant Q transform
%    ERBLETT                - Erb-let transform
%    IERBLETT               - Inverse Erb-let transform
%    INSDGFB                - Say what?
%
%  Filter generators
%    ERBFILTERS             - ERB-spaced filters
%  
%  Window construction and bounds
%    FILTERBANKDUAL         - Canonical dual filters
%    FILTERBANKTIGHT        - Canonical tight filters
%    FILTERBANKREALDUAL     - Canonical dual filters for real signals
%    FILTERBANKREALTIGHT    - Canonical tight filters for real signals
%    FILTERBANKBOUNDS       - Frame bounds of filter bank
%    FILTERBANKREALBOUNDS   - Frame bounds of filter bank for real signals
%    FILTERBANKRESPONSE     - Total frequency response
%
%  Plots
%    PLOTFILTERBANK         - Plot normal/uniform filter bank coefficients
%
%  For help, bug reports, suggestions etc. please send email to
%  ltfat-help@lists.sourceforge.net

